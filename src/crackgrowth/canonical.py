"""The canonical Milestone 1 wing-skin sanity case.

An idealized aircraft wing-skin through crack under constant-amplitude loading.
Every number here is an explicit engineering assumption, not a design input:

* ``Y = 1.0``            -- infinite-plate through crack; no finite-width effect.
* ``sigma_max = 120 MPa``, ``sigma_min = 20 MPa`` -- a plausible tension-dominated
  ground-air-ground-like cycle, giving ``delta_sigma = 100 MPa`` and ``R = 1/6``.
* ``a0 = 1.0 mm``        -- an assumed initial flaw size, NOT a detectability limit.
* ``a_target = 10.0 mm`` -- an IMPOSED target crack length. In Milestone 1 this is
  not derived from fracture toughness and is NOT a critical crack size.
* Paris curve           -- illustrative; see
  :data:`crackgrowth.paris.ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS`.
* ``K_IC = 30 MPa*sqrt(m)`` -- illustrative (Milestone 2); gives a critical
  crack length of 19.89 mm for this cycle. See
  :data:`crackgrowth.fracture.ILLUSTRATIVE_ALUMINIUM_LIKE_TOUGHNESS`.

``CANONICAL_TARGET_CRACK_LENGTH`` remains the Milestone 1 IMPOSED target and
is retained unchanged for regression. Milestone 2 derives its own endpoint
from toughness instead; the two are deliberately different.
"""

from __future__ import annotations

from .finite_width import FiniteWidthCenterCrack
from .fracture import ILLUSTRATIVE_ALUMINIUM_LIKE_TOUGHNESS
from .threshold import ILLUSTRATIVE_ALUMINIUM_LIKE_THRESHOLD
from .geometry import INFINITE_PLATE_THROUGH_CRACK
from .loading import StressCycle
from .paris import ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS

__all__ = [
    "CANONICAL_GEOMETRY",
    "CANONICAL_CYCLE",
    "CANONICAL_PARIS_LAW",
    "CANONICAL_INITIAL_CRACK_LENGTH",
    "CANONICAL_TARGET_CRACK_LENGTH",
    "CANONICAL_FRACTURE_TOUGHNESS",
    "CANONICAL_PLATE_WIDTH",
    "CANONICAL_FINITE_WIDTH_GEOMETRY",
    "CANONICAL_GROWTH_THRESHOLD",
]

CANONICAL_GEOMETRY = INFINITE_PLATE_THROUGH_CRACK
CANONICAL_CYCLE = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)
CANONICAL_PARIS_LAW = ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS

#: Assumed initial crack length [m].
CANONICAL_INITIAL_CRACK_LENGTH = 1.0e-3

#: Imposed target crack length [m]. NOT a fracture-toughness-derived critical size.
#: Retained from Milestone 1 for regression; Milestone 2 supersedes it with
#: a toughness-derived critical size.
CANONICAL_TARGET_CRACK_LENGTH = 10.0e-3

#: Canonical illustrative mode-I fracture toughness (Milestone 2).
CANONICAL_FRACTURE_TOUGHNESS = ILLUSTRATIVE_ALUMINIUM_LIKE_TOUGHNESS

#: Canonical panel width [m] (Milestone 3).
#:
#: Chosen after auditing W = 40, 50, 75, 100, 150, 200 and 500 mm against the
#: canonical cycle and flaw. At W = 100 mm the initial flaw is safely small
#: (a0/W = 0.01, Y(a0) = 1.00025), the critical crack sits at 17.09 mm with a
#: total crack length of only 34 % of the panel width, and the finite-width
#: effect is noticeable but not pathological: 14 % smaller critical crack and
#: 4.4 % less life than the infinite plate. Narrower panels give a more
#: dramatic penalty but push the total crack past half the width, where the
#: unmodelled net-section behaviour would realistically govern.
CANONICAL_PLATE_WIDTH = 100.0e-3

#: Canonical finite-width centre-cracked panel (Milestone 3).
#: ``a`` is the HALF crack length throughout.
CANONICAL_FINITE_WIDTH_GEOMETRY = FiniteWidthCenterCrack(
    plate_width=CANONICAL_PLATE_WIDTH
)

#: Canonical illustrative crack-growth threshold (Milestone 4).
#:
#: Chosen after auditing delta_K_th = 2, 3, 4, 5, 6, 8, 10 MPa*sqrt(m) against
#: the canonical finite-width case, where delta_K(a0) = 5.6064 MPa*sqrt(m):
#: values up to 5 leave the initial crack active, 6 and above arrest it.
#: 4 MPa*sqrt(m) gives a threshold ratio of 1.40 at a0 and a threshold crack
#: size of 0.51 mm, so a genuine active-growth interval exists.
CANONICAL_GROWTH_THRESHOLD = ILLUSTRATIVE_ALUMINIUM_LIKE_THRESHOLD
