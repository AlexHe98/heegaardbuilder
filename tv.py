"""
Experiment with computing Turaev-Viro invariants.
"""
from timeit import default_timer
from regina import *


def timeTV( tri, r ):
    """
    Returns the number of seconds that it takes to compute the Turaev-Viro
    invariant TV_r of the given triangulation tri.
    """
    start = default_timer()
    tri.turaevViro( r, True, ALG_TREEWIDTH )
    return default_timer() - start
