"""Keep reference provenance separate from numerical agreement and calibration."""

from pathlib import Path

UNVERIFIED_REFERENCE = 'REFERENCE COMPARISON / PROVENANCE UNVERIFIED'
SYNTHETIC_REFERENCE = 'SYNTHETIC REFERENCE / NOT EXPERIMENTAL'
MODEL_BASED = 'MODEL-BASED / NOT EXPERIMENTALLY VALIDATED'

# These origins are stated in the files' headers; no original measurement
# records or digitization projects have been independently checked yet.
SYNTHETIC_FILES = {
    'ingaas_pn_experimental_iv.csv',
    'si_pn_experimental_cv.csv',
    'si_pn_experimental_iv.csv',  # Header defines a Shockley model, not a measurement source.
}


def reference_status(paths):
    names = [Path(str(path)).name for path in paths if path]
    if not names:
        return MODEL_BASED
    if all(name in SYNTHETIC_FILES for name in names):
        return SYNTHETIC_REFERENCE
    return UNVERIFIED_REFERENCE


def reviewed_status(claim):
    """Legacy claims are not proof that a model passed independent validation."""
    if claim in ('EXPERIMENTALLY VALIDATED', 'EXPERIMENTAL', 'PARTIALLY VALIDATED'):
        return UNVERIFIED_REFERENCE
    if claim == 'FITTED TO EXPERIMENT':
        return 'CALIBRATED MODEL / PROVENANCE UNVERIFIED'
    return claim or MODEL_BASED
