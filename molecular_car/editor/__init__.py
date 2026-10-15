"""Independent molecule-workshop file/geometry service.

Only ``dispatch(path, payload)`` is needed by the host HTTP server. None of these
operations loads, stops, or mutates a race. Editor identity always follows array
indices; SMILES and graph matching are never used to recover atom identity.
"""

from .service import EditorError, dispatch

__all__ = ["EditorError", "dispatch"]
