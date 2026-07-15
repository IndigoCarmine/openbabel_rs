//! Error type for the safe OpenBabel API.

use std::fmt;

/// Something went wrong talking to OpenBabel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The input data could not be parsed in the requested format.
    Parse {
        /// The format id that was attempted (e.g. `"smi"`).
        format: String,
    },
    /// OpenBabel does not recognize the requested format id.
    UnknownFormat {
        /// The unrecognized format id.
        format: String,
    },
    /// A SMARTS pattern failed to compile.
    InvalidSmarts {
        /// The offending pattern.
        pattern: String,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Parse { format } => {
                write!(f, "failed to parse input as format {format:?}")
            }
            Error::UnknownFormat { format } => {
                write!(f, "unknown OpenBabel format {format:?}")
            }
            Error::InvalidSmarts { pattern } => {
                write!(f, "invalid SMARTS pattern {pattern:?}")
            }
        }
    }
}

impl std::error::Error for Error {}
