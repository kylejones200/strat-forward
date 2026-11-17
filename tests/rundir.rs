//! Output directory replacement, exercised directly: the swap is atomic where the platform supports it, and
//! readers never see the target missing or its files missing.

#[cfg(unix)]
use std::os::fd::{AsRawFd, FromRawFd};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;

use strat_forward::rundir::RunDir;

fn workdir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("rundir-{name}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn id(n: usize) -> String {
    format!("00000000-0000-4000-8000-{n:012}")
}

fn commit_one(out: &Path, n: usize) -> strat_forward::rundir::Committed {
    let d = RunDir::prepare(out, &id(n)).unwrap();
    std::fs::write(d.path().join("summary.json"), format!("{{\"n\":{n}}}")).unwrap();
    d.commit().unwrap()
}

#[test]
fn replacing_an_existing_directory_is_one_atomic_exchange_on_macos_and_linux() {
    let out = workdir("atomic").join("out");
    assert!(
        commit_one(&out, 1).atomic,
        "a fresh target is a single rename"
    );
    let second = commit_one(&out, 2);
    if cfg!(any(target_os = "macos", target_os = "linux")) {
        // Holds on APFS and on Linux filesystems with RENAME_EXCHANGE (ext4, xfs, btrfs, tmpfs).
        assert!(
            second.atomic,
            "expected an atomic exchange on this platform"
        );
    }
    let previous = out.with_file_name("out.previous");
    assert_eq!(second.previous.as_deref(), Some(previous.as_path()));
    assert_eq!(
        std::fs::read_to_string(previous.join("summary.json")).unwrap(),
        "{\"n\":1}",
        "the replaced outputs are kept as out.previous"
    );
    assert_eq!(
        std::fs::read_to_string(out.join("summary.json")).unwrap(),
        "{\"n\":2}"
    );
    let parent = out.parent().unwrap();
    let mut names: Vec<String> = std::fs::read_dir(parent)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert_eq!(
        names,
        ["out", "out.lock", "out.previous"],
        "no staging directories left"
    );
}

#[test]
fn a_reader_never_finds_the_directory_missing_during_replacement() {
    let out = workdir("reader").join("out");
    commit_one(&out, 0);
    let stop = Arc::new(AtomicBool::new(false));
    let (misses, checks) = (Arc::new(AtomicUsize::new(0)), Arc::new(AtomicUsize::new(0)));
    let reader = {
        let (out, stop, misses, checks) =
            (out.clone(), stop.clone(), misses.clone(), checks.clone());
        std::thread::spawn(move || {
            while !stop.load(Ordering::Relaxed) {
                if std::fs::metadata(out.join("summary.json")).is_err() {
                    misses.fetch_add(1, Ordering::Relaxed);
                }
                checks.fetch_add(1, Ordering::Relaxed);
            }
        })
    };
    let mut all_atomic = true;
    for n in 1..=200 {
        all_atomic &= commit_one(&out, n).atomic;
        // Replaced outputs are removed one replacement later (as `out.previous`). Back-to-back commits
        // ~1 ms apart let a preempted reader still hold that directory when it goes (5 of 100 runs); real
        // runs are minutes apart. 5 ms is still far shorter than any run.
        std::thread::sleep(std::time::Duration::from_millis(5));
    }
    stop.store(true, Ordering::Relaxed);
    reader.join().unwrap();
    assert!(checks.load(Ordering::Relaxed) > 0);
    if all_atomic {
        assert_eq!(misses.load(Ordering::Relaxed), 0, "the target went missing");
    }
}

#[cfg(unix)]
#[test]
fn a_reader_that_resolved_the_target_before_a_replacement_can_still_read_it() {
    // Looking up a path and reading it are separate steps. Hold the old directory open, as a reader does
    // between the two; it must stay readable after the swap (it was deleted at once before, giving ENOENT).
    let out = workdir("held").join("out");
    commit_one(&out, 1);
    let held = std::fs::File::open(&out).unwrap();
    commit_one(&out, 2);
    let name = std::ffi::CString::new("summary.json").unwrap();
    // SAFETY: a valid directory descriptor and a NUL-terminated name.
    let fd = unsafe { libc::openat(held.as_raw_fd(), name.as_ptr(), libc::O_RDONLY) };
    assert!(fd >= 0, "{}", std::io::Error::last_os_error());
    // SAFETY: fd was just opened above and is owned here.
    let mut f = unsafe { std::fs::File::from_raw_fd(fd) };
    let mut s = String::new();
    std::io::Read::read_to_string(&mut f, &mut s).unwrap();
    assert_eq!(s, "{\"n\":1}");
}

#[test]
fn a_held_lock_refuses_a_second_prepare_and_is_released_on_drop() {
    let out = workdir("lock").join("out");
    let first = RunDir::prepare(&out, &id(1)).unwrap();
    let err = RunDir::prepare(&out, &id(2))
        .err()
        .expect("second prepare was accepted");
    assert!(
        err.to_string().contains("another strat-forward run"),
        "{err}"
    );
    first.abandon().unwrap();
    // The lock went with the dropped RunDir.
    RunDir::prepare(&out, &id(3)).unwrap().commit().unwrap();
}

#[test]
fn a_previous_directory_it_did_not_write_is_never_removed() {
    let out = workdir("foreign-previous").join("out");
    commit_one(&out, 1);
    let previous = out.with_file_name("out.previous");
    std::fs::create_dir(&previous).unwrap();
    std::fs::write(previous.join("notes.txt"), "mine").unwrap();
    let d = RunDir::prepare(&out, &id(2)).unwrap();
    let staging = d.path().to_path_buf();
    std::fs::write(staging.join("summary.json"), "{\"n\":2}").unwrap();
    let err = d.commit().expect_err("commit removed a foreign directory");
    assert!(err.to_string().contains("refusing to remove"), "{err}");
    assert_eq!(
        std::fs::read_to_string(previous.join("notes.txt")).unwrap(),
        "mine"
    );
    // The current outputs are untouched and the new ones are kept where the error says.
    assert_eq!(
        std::fs::read_to_string(out.join("summary.json")).unwrap(),
        "{\"n\":1}"
    );
    assert!(err.to_string().contains(&staging.display().to_string()));
    assert!(staging.join("summary.json").is_file());
}
