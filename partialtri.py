"""
A class to represent a partial triangulation built by layering tetrahedra
according to a boundary curve.
"""
from regina import *
from normalcurve import NormalCurve


class PartialTriangulation(NormalCurve):
    """
    A partial triangulation with a boundary curve, represented as a union
    of a normal curve and a collection of boundary edges.
    """
    def __init__( self, tri, weights, resolvedEdgeIndices=set() ):
        """
        Create a partial triangulation with a boundary curve, described by a
        normal curve and (optionally) some boundary edges.
        """
        super().__init__( tri, weights )

        # Set of edge indices to which we have resolved curves.
        self._resolvedEdges = set(resolvedEdgeIndices)

    def countUnresolved(self):
        """
        Returns the number of unresolved components of this normal curve.
        """
        return len( self._roots )

    def countResolved(self):
        """
        Returns the number of resolved components of this normal curve.
        """
        return len( self._resolvedEdges )

    def countComponents(self):
        """
        Returns the number of components of this normal curve.
        """
        return self.countUnresolved() + self.countResolved()

    def getResolvedEdgeIndices(self):
        """
        Returns a copy of the set of resolved edge indices.
        """
        return set( self._resolvedEdges )

    #TODO Every time we change the triangulation, we need to update resolved
    #   edge indices.
    def layerOn( self, e ):
        """
        Layer a new tetrahedron across the edge e.

        Pre-condition:
        --> e is a boundary edge of self.triangulation().
        """
        # Keep track of how resolved edges get renumbered.
        temp = []
        #TODO
        pass

    def flipEdge( self, e ):
        #TODO
        pass

    def resolveComponent( self, index ):
        """
        Checks whether it is possible to resolve the requested component, and
        if so resolves this component.
        """
        resEdgeInds = self._resolvableEdges[index]
        if not resEdgeInds:
            return False

        # Clear the edge, and then resolve the requested component by simply
        # deleting its intersection points.
        edgeInd = resEdgeInds[0]
        self._clearEdge( self._tri.edge(edgeInd) )
        newWeights = list( self._weights )
        for p, _, _ in self._traverseComponentImpl(index):
            newWeights[ p[0] ] -= 1
        self._resolvedEdges.add(edgeInd)
        self._processWeights(newWeights)
        #TODO Double-check and test.
    #TODO
    pass


# Test code.
if __name__ == "__main__":
    #TODO
