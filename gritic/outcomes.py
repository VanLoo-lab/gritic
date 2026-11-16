"""Structured outcomes for valid samples outside GRITIC's analysis domain."""

EXIT_INELIGIBLE = 3


class IneligibleSampleError(ValueError):
    """Valid biological input cannot support the requested analysis.

    ``reason`` is a stable machine-readable identifier. Malformed input and
    numerical or operational failures must use their ordinary exception types.
    """

    def __init__(self, reason: str, message: str):
        super().__init__(message)
        self.reason = reason

    def as_dict(self, sample_id: str) -> dict:
        return {
            'format': 'gritic-outcome', 'status': 'ineligible',
            'sample_id': sample_id, 'reason': self.reason, 'message': str(self),
        }
