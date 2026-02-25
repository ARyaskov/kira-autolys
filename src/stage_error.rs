use std::fmt;

#[derive(Debug)]
pub enum StageError {
    Io(std::io::Error),
    Format(String),
    Validation(String),
    FeatureDisabled(String),
}

impl fmt::Display for StageError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "io error: {err}"),
            Self::Format(msg) => write!(f, "format error: {msg}"),
            Self::Validation(msg) => write!(f, "validation error: {msg}"),
            Self::FeatureDisabled(msg) => write!(f, "feature disabled: {msg}"),
        }
    }
}

impl std::error::Error for StageError {}

impl From<std::io::Error> for StageError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}
