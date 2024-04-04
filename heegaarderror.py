"""
Custom exceptions for the HeegaardBuilder class.
"""
from regina import *


class HeegaardError(Exception):
    """
    Raised when a HeegaardBuilder object encounters an error.
    """
    pass


class EmptyHeegaardBuilder(HeegaardError):
    """
    Raised when attempting to perform an operation on a HeegaardBuilder
    object that is empty.
    """
    def __init__(self):
        super().__init__()


class BadTriangulation(HeegaardError):
    """
    Raised when a HeegaardBuilder object is given a bad triangulation.
    """
    def __init__( self, badProperty ):
        msg = "Triangulation is {}.".format(badProperty)
        super().__init__(msg)


class BadBouquet(HeegaardError):
    """
    Raised when a HeegaardBuilder object detects a bad filling bouquet.
    """
    pass


class WrongNumberOfWeights(BadBouquet):
    """
    Raised when a HeegaardBuilder object is given a triangulation and a list
    of edge weights such that the number of edges in the triangulation
    doesn't match the number of weights.
    """
    def __init__( self, numWeights, numEdges ):
        msg = ( "Given {} edge weights, but the ".format( numWeights ) +
                "triangulation has {} edges.".format( numEdges ) )
        super().__init__(msg)


class NegativeEdgeWeight(BadBouquet):
    """
    Raised when a HeegaardBuilder object is given a negative edge weight.
    """
    def __init__( self, edgeIndex, weight ):
        msg = ( "Edge weights must be non-negative, but edge " +
                "{} has weight {}.".format( edgeIndex, weight ) )
        super().__init__(msg)


class WeightOnInternalEdge(BadBouquet):
    """
    Raised when a HeegaardBuilder object is given a nonzero edge weight on an
    internal edge.
    """
    def __init__( self, edgeIndex, weight ):
        msg = ( "Edge {} is internal, so its weight ".format(edgeIndex) +
                "must be 0, not {}.".format(weight) )
        super().__init__(msg)


class ResolvedInternalEdge(BadBouquet):
    """
    Raised when a HeegaardBuilder object is asked to assign an internal edge
    as a resolved edge.
    """
    def __init__( self, edgeIndex ):
        msg = ( "Edge {} is internal, so it cannot be ".format(edgeIndex) +
                "assigned as a resolved edge." )
        super().__init__(msg)


class WeightOnResolvedEdge(BadBouquet):
    """
    Raised when a HeegaardBuilder object is given a nonzero edge weight on an
    edge that is supposed to be resolved.
    """
    def __init__( self, edgeIndex, weight ):
        msg = ( "Edge {} is resolved, so its weight ".format(edgeIndex) +
                "must be 0, not {}.".format(weight) )
        super().__init__(msg)


class FailedMatchingConstraints(BadBouquet):
    """
    Raised when a HeegaardBuilder object is given edge weights that fail to
    satisfy the matching constraints for a filling bouquet.
    """
    def __init__( self, faceIndex ):
        msg = ( "Edge weights fail to satisfy the matching constraints " +
                "in triangle {}.".format(faceIndex) )
        super().__init__(msg)


class NotCombinatoriallyAdmissible(BadBouquet):
    """
    Raised when a HeegaardBuilder object detects that the stored filling
    bouquet is not combinatorially admissible.
    """
    def __init__(self):
        msg = ( "The filling bouquet is not combinatorially admissible." )
        super().__init__(msg)


class NotTopologicallyAdmissible(BadBouquet):
    """
    Raised when a HeegaardBuilder object detects that the stored filling
    bouquet is not (topologically) admissible.
    """
    pass


class TransversePetals(NotTopologicallyAdmissible):
    """
    Raised when a HeegaardBuilder object encounters a pair of edges that form
    resolved filling petals that meet transversely.
    """
    def __init__( self, myEdgeInd, yourEdgeInd ):
        msg = ( "Edges {} and {} form a ".format( myEdgeInd, yourEdgeInd ) +
                "pair of resolved filling petals that meet transversely." )
        super().__init__(msg)


class DisconnectedComplement(NotTopologicallyAdmissible):
    """
    Raised when a HeegaardBuilder object detects that cutting along the
    stored filling bouquet would disconnect the boundary surface of the
    stored triangulation.
    """
    def __init__( self, numComponents ):
        msg = ( "After cutting along the filling bouquet, the boundary " +
                "surface splits into {} components.".format(numComponents) )
        super().__init__(msg)


class NormalCurveAfterResolving(BadBouquet):
    """
    Raised when a HeegaardBuilder object detects a normal curve that is left
    over after resolving all filling petals.
    """
    def __init__(self):
        msg = ( "A normal curve was left over after resolving all " +
                "filling petals." )
        super().__init__(msg)


class IsotopicPetals(BadBouquet):
    """
    Raised when a HeegaardBuilder object detects an edge that is incident to
    two or more filling petals that are all isotopic to each other.
    """
    def __init__( self, edgeIndex, numIsotopicPetals ):
        msg = ( "Edge {} meets {} ".format( edgeIndex, numIsotopicPetals ) +
                "filling petals that are isotopic to each other." )


class NoReducibleEdge(BadBouquet):
    """
    Raised when a HeegaardBuilder object cannot find a reducible edge, which
    can only happen if there is a normal curve.
    """
    def __init__(self):
        msg = ( "Could not find a reducible edge, which means that there " +
                "must be a normal curve." )
        super().__init__(msg)
