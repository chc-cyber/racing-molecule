"""Shared race geometry in Angstrom; the playable rectangle is 20 nm square.

Longitudinal soft walls sit 1 nm beyond the nominal track ends, retaining
the original run-up and finish cushion. Geometry changes do not rescale the
crystal lattice, the drive field, or the integrator's time and length units.
Race markers use the post-relaxation COM origin plus this nominal track length;
that scoring offset does not translate the physical substrate or fields.
"""

TRACK_START_X = 0.0
TRACK_FINISH_X = 200.0
LANE_HALF_WIDTH = 100.0
TRACK_X_MIN = TRACK_START_X - 10.0
TRACK_X_MAX = TRACK_FINISH_X + 10.0

TRACK = {"start_x": TRACK_START_X, "finish_x": TRACK_FINISH_X,
         "half_width": LANE_HALF_WIDTH, "surface_z": 0.0}
