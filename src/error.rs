/// Possible RTM errors.
#[derive(Debug)]
pub enum RtmError {
    /// The inputs don't have the expected shape(s)
    InconsistentInputs,
    /// Couldn't find the surface index
    NoSurface,
    /// An array is not contiguous when it was assumed to be
    NotContiguous,
    /// The operation was aborted early
    Cancelled,
}

impl std::fmt::Display for RtmError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RtmError::InconsistentInputs => {
                write!(f, "inputs to RTM have the wrong shape")
            }
            RtmError::NoSurface => {
                write!(f, "couldn't find the surface index")
            }
            RtmError::NotContiguous => write!(f, "array slice not contiguous in memory"),
            RtmError::Cancelled => write!(f, "operation cancelled early"),
        }
    }
}

impl std::error::Error for RtmError {}
