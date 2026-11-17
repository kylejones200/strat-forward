//! The crate's error type.

use std::path::PathBuf;

#[derive(Debug)]
pub enum Error {
    /// The configuration is invalid. Nothing was computed.
    Config(String),
    /// A file or directory could not be read or written.
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    /// The model's own consistency check failed: a bug, not bad input.
    Invariant(String),
    /// The call is not allowed in the model's current state (for example, running a model twice).
    State(String),
    /// The requested product cannot be built from this run (for example, no reservoir window qualifies).
    Unsupported(String),
}

pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    pub(crate) fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Error::Io {
            path: path.into(),
            source,
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Config(m) => write!(f, "invalid config: {m}"),
            Error::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Error::Invariant(m) => write!(f, "internal consistency check failed: {m}"),
            Error::State(m) => f.write_str(m),
            Error::Unsupported(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}
