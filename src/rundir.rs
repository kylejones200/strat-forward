//! Output directories that are replaced whole, or not at all.
//!
//! A run writes into a fresh staging directory beside the target. On success the staging directory takes the
//! target's place; on failure it is kept as `<target>.failed-<run id>` and the target is not touched. So a
//! failed run never destroys the last good outputs or leaves old and new files mixed.
//!
//! Replacing an existing target is one atomic exchange where the filesystem supports it (`renamex_np` with
//! `RENAME_SWAP` on macOS, `renameat2` with `RENAME_EXCHANGE` on Linux): readers see the old directory or the
//! new one, never neither. Elsewhere it falls back to two renames (old aside, staging in), and the target is
//! briefly missing between them.
//!
//! The replaced outputs are kept as `<target>.previous` until the next successful run, not deleted at once.
//! Looking up a path and reading it are separate steps, so a reader can resolve `<target>` to the old
//! directory just before the swap; deleting that directory straight after made such a reader find the
//! target missing (ENOENT on the target itself: 47 to 76 in 20,000 replacements in a stress probe under CPU
//! load, none with the deletion left out). Kept, it stays readable for a whole run. `<target>.previous` is
//! removed before the next swap only if this tool wrote it. It also leaves one generation to roll back to,
//! at the cost of a second copy on disk.
//!
//! A run holds an exclusive lock on `<target>.lock` from `prepare` until it commits or abandons, so two runs
//! cannot write the same target at once; the second is refused. The operating system releases the lock if
//! the process dies, so a crash never leaves a stale lock. The lock file itself stays on disk (removing it
//! would let a waiting run lock a file that is no longer the one others open).
//!
//! A target that exists must be one this tool wrote (it holds the marker file, or a `summary.json` from this
//! engine) or be empty. Anything else is refused, so `--out` can never replace an unrelated directory.

use std::fs::File;
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};

/// Marker file written into every committed output directory.
pub const MARKER: &str = ".strat-forward";

pub struct RunDir {
    target: PathBuf,
    staging: PathBuf,
    run_id: String,
    // Held for the life of the run; dropping it releases the lock.
    _lock: File,
}

/// What `commit` left behind.
#[derive(Debug)]
pub struct Committed {
    pub path: PathBuf,
    /// True when the target appeared in one step: a fresh rename, or an atomic exchange with the old outputs.
    pub atomic: bool,
    /// Where the outputs it replaced are kept (`<target>.previous`), if there were any.
    pub previous: Option<PathBuf>,
}

fn sibling(target: &Path, suffix: &str) -> PathBuf {
    let name = target
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "out".into());
    target.with_file_name(format!("{name}.{suffix}"))
}

/// The lock file guarding `target`.
pub fn lock_path(target: &Path) -> PathBuf {
    sibling(target, "lock")
}

/// Where the outputs replaced by the last successful run are kept.
pub fn previous_path(target: &Path) -> PathBuf {
    sibling(target, "previous")
}

/// True when `dir` may be replaced: missing, empty, or written by this tool.
pub fn replaceable(dir: &Path) -> Result<bool> {
    if !dir.exists() {
        return Ok(true);
    }
    if !dir.is_dir() {
        return Ok(false);
    }
    if dir.join(MARKER).is_file() {
        return Ok(true);
    }
    // Outputs written before the marker existed.
    if let Ok(raw) = std::fs::read(dir.join("summary.json")) {
        if let Ok(v) = serde_json::from_slice::<serde_json::Value>(&raw) {
            if v.get("engine").and_then(|e| e.as_str()) == Some(env!("CARGO_PKG_NAME")) {
                return Ok(true);
            }
        }
    }
    let mut entries = std::fs::read_dir(dir).map_err(|e| Error::io(dir, e))?;
    Ok(entries.next().is_none())
}

/// Exchange two existing paths in one step. `Ok(false)` when this platform or filesystem cannot.
fn exchange(a: &Path, b: &Path) -> std::io::Result<bool> {
    #[cfg(any(target_os = "macos", target_os = "linux"))]
    {
        use std::ffi::CString;
        use std::os::unix::ffi::OsStrExt;
        let ca = CString::new(a.as_os_str().as_bytes())?;
        let cb = CString::new(b.as_os_str().as_bytes())?;
        // SAFETY: both arguments are valid NUL-terminated strings that outlive the call.
        #[cfg(target_os = "macos")]
        let rc = unsafe { libc::renamex_np(ca.as_ptr(), cb.as_ptr(), libc::RENAME_SWAP) };
        // SAFETY: as above; AT_FDCWD resolves relative paths against the working directory.
        #[cfg(target_os = "linux")]
        let rc = unsafe {
            libc::renameat2(
                libc::AT_FDCWD,
                ca.as_ptr(),
                libc::AT_FDCWD,
                cb.as_ptr(),
                libc::RENAME_EXCHANGE,
            )
        };
        if rc == 0 {
            return Ok(true);
        }
        let err = std::io::Error::last_os_error();
        // Not supported by this kernel or filesystem: let the caller fall back.
        match err.raw_os_error() {
            Some(libc::ENOTSUP) | Some(libc::EINVAL) | Some(libc::ENOSYS) => Ok(false),
            _ => Err(err),
        }
    }
    #[cfg(not(any(target_os = "macos", target_os = "linux")))]
    {
        let _ = (a, b);
        Ok(false)
    }
}

impl RunDir {
    /// Lock the target, check it, and create the staging directory.
    pub fn prepare(target: &Path, run_id: &str) -> Result<Self> {
        let target = if target.is_absolute() {
            target.to_path_buf()
        } else {
            std::env::current_dir()
                .map_err(|e| Error::io(".", e))?
                .join(target)
        };
        if let Some(parent) = target.parent() {
            std::fs::create_dir_all(parent).map_err(|e| Error::io(parent, e))?;
        }
        let lock_file = lock_path(&target);
        let lock = File::options()
            .create(true)
            .truncate(false)
            .write(true)
            .open(&lock_file)
            .map_err(|e| Error::io(&lock_file, e))?;
        match lock.try_lock() {
            Ok(()) => {}
            Err(std::fs::TryLockError::WouldBlock) => {
                return Err(Error::State(format!(
                    "another strat-forward run is writing {} (it holds {})",
                    target.display(),
                    lock_file.display()
                )))
            }
            Err(std::fs::TryLockError::Error(e)) => return Err(Error::io(&lock_file, e)),
        }
        if !replaceable(&target)? {
            return Err(Error::State(format!(
                "refusing to replace {}: it exists, is not empty, and was not written by strat-forward",
                target.display()
            )));
        }
        let staging = sibling(&target, &format!("staging-{run_id}"));
        // create_dir (not _all) fails if a staging directory for this run id already exists.
        std::fs::create_dir(&staging).map_err(|e| Error::io(&staging, e))?;
        Ok(RunDir {
            target,
            staging,
            run_id: run_id.to_string(),
            _lock: lock,
        })
    }

    /// Where the run writes.
    pub fn path(&self) -> &Path {
        &self.staging
    }

    /// Where the outputs will be once committed.
    pub fn target(&self) -> &Path {
        &self.target
    }

    /// Put the staged outputs in the target's place. The outputs it replaces are kept as `<target>.previous`
    /// until the next successful run replaces them in turn.
    pub fn commit(self) -> Result<Committed> {
        let marker = self.staging.join(MARKER);
        std::fs::write(&marker, format!("{{\"run_id\":\"{}\"}}\n", self.run_id))
            .map_err(|e| Error::io(&marker, e))?;
        if !self.target.exists() {
            std::fs::rename(&self.staging, &self.target)
                .map_err(|e| Error::io(&self.staging, e))?;
            return Ok(Committed {
                path: self.target.clone(),
                atomic: true,
                previous: None,
            });
        }
        // Re-check: the lock keeps other strat-forward runs out, not other programs.
        if !replaceable(&self.target)? {
            return Err(Error::State(format!(
                "{} changed during the run and is no longer replaceable; outputs left in {}",
                self.target.display(),
                self.staging.display()
            )));
        }
        // Make room for the outputs being replaced. Only a directory this tool wrote is removed.
        let previous = previous_path(&self.target);
        if previous.exists() {
            if !replaceable(&previous)? {
                return Err(Error::State(format!(
                    "refusing to remove {}: it was not written by strat-forward; outputs left in {}",
                    previous.display(),
                    self.staging.display()
                )));
            }
            std::fs::remove_dir_all(&previous).map_err(|e| Error::io(&previous, e))?;
        }
        let atomic =
            exchange(&self.staging, &self.target).map_err(|e| Error::io(&self.staging, e))?;
        if atomic {
            // The staging path now holds the replaced outputs.
            std::fs::rename(&self.staging, &previous).map_err(|e| Error::io(&self.staging, e))?;
        } else {
            std::fs::rename(&self.target, &previous).map_err(|e| Error::io(&self.target, e))?;
            if let Err(e) = std::fs::rename(&self.staging, &self.target) {
                // Put the previous outputs back.
                let _ = std::fs::rename(&previous, &self.target);
                return Err(Error::io(&self.staging, e));
            }
        }
        Ok(Committed {
            path: self.target.clone(),
            atomic,
            previous: Some(previous),
        })
    }

    /// Keep what the failed run wrote as `<target>.failed-<run id>`; the target is not touched.
    pub fn abandon(self) -> Result<PathBuf> {
        let failed = sibling(&self.target, &format!("failed-{}", self.run_id));
        std::fs::rename(&self.staging, &failed).map_err(|e| Error::io(&self.staging, e))?;
        Ok(failed)
    }
}
