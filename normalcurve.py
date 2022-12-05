"""
A class to represent a normal curve in a boundary surface.
"""
from regina import *


class NormalCurve:
    """
    A normal curve in the boundary of a 3-manifold triangulation.
    """
    def __init__( self, tri, weights ):
        """
        Create a normal curve with the given edge weights in the boundary of
        the given triangulation.
        """
        # Check that we have the correct number of weights.
        if len(weights) != tri.countEdges():
            raise ValueError( "Wrong number of edge weights!" )

        self._tri = tri
        self._processWeights(weights)

    def _processWeights( self, weights ):
        self._weights = tuple(weights)

        # Check that the given weights satisfy the matching constraints.
        if not self._checkMatching():
            raise ValueError(
                    "Edge weights fail the matching constraints!" )

        # Count components using union-find.
        self._points = set()
        self._initialisePoints()
        self._parent = { p: p for p in self._points }
        self._size = { p: 1 for p in self._points }
        # As a by-product of this, we will also:
        # - compute coordinates for this normal curve in terms of its arcs;
        # - assign an index to each component of this normal curve; and
        # - record precisely how this normal curve traverses the boundary of
        #   our triangulation.
        self._arcCoords = [None] * self._tri.countTriangles()
        self._roots = list( self._points )
        self._adjPoints = { p: list() for p in self._points }
        self._adjCutVerts = { p: list() for p in self._points }
        self._countComponentsImpl()

        # Find where all the switches happen.
        self._switch = { p: False for p in self._points }
        self._findSwitches()

    def _checkMatching(self):
        for face in self._tri.triangles():
            if not face.isBoundary():
                continue

            # Check that the edge weights around this boundary triangle
            # satisfy the matching constraints:
            # - Sum of weights is even.
            # - No edge has weight too large to match up with the other two
            #   edges.
            weights = [ self._weights[ face.edge(i).index() ]
                    for i in range(3) ]
            if sum(weights) % 2 != 0:
                return False
            for i in range(3):
                if weights[i] > weights[ (i-1)%3 ] + weights[ (i+1)%3 ]:
                    return False

        return True

    def _initialisePoints(self):
        # Represent each intersection point using a pair consisting of the
        # edge index and the order of the point along the edge.
        for edgeInd in range( self._tri.countEdges() ):
            if not self._tri.edge(edgeInd).isBoundary():
                continue
            for p in range( self._weights[edgeInd] ):
                self._points.add( ( edgeInd, p ) )

    def _find( self, p ):
        # Find operation for union-find.
        while self._parent[p] != p:
            # Use path-halving: replace the current parent of p with its
            # current grandparent, and then go to this grandparent.
            self._parent[p] = self._parent[ self._parent[p] ]
            p = self._parent[p]
        return p

    def _union( self, p, q ):
        # Union operation for union-find.
        p = self._find(p)
        q = self._find(q)
        if p == q:
            return

        # Make sure that the set containing p has size greater than or equal
        # to the set containing q.
        if self._size[p] < self._size[q]:
            p, q = q, p

        # Take the union by making p the new root. Also update the size and
        # the set of roots.
        self._parent[q] = p
        self._size[p] = self._size[p] + self._size[q]
        self._roots.remove(q)

    def _countComponentsImpl(self):
        # Piece together the components of our normal curve by looking at one
        # boundary triangle at a time.
        for face in self._tri.triangles():
            if not face.isBoundary():
                continue
            verts = face.embedding(0).vertices()
            tet = face.embedding(0).tetrahedron()

            # For each i in {0,1,2}, compute the following:
            # - The index of the edge opposite vertex i of this face.
            # - The weight of this edge.
            # - The start vertex of this edge.
            # - The embedding of this edge in tet.
            edgeInds = []
            weights = []
            startVerts = []
            for i in range(3):
                # Get the edge opposite vertex i of this face, and remember
                # both its index and its weight.
                endpoints = { verts[ii] for ii in range(3) }
                endpoints.remove( verts[i] )
                edgeInd = tet.edge( *endpoints ).index()
                edgeInds.append(edgeInd)
                weights.append( self._weights[edgeInd] )

                # To get the start vertex, find the correct corresponding
                # edge embedding.
                edge = self._tri.edge(edgeInd)
                embs = [ edge.embedding(0),
                        edge.embedding( edge.degree() - 1 ) ]
                for emb in embs:
                    if emb.tetrahedron() != tet:
                        continue
                    if endpoints == { emb.vertices()[j] for j in range(2) }:
                        startVerts.append( emb.vertices()[0] )
                        break

            # For each i in {0,1,2}, use the arcs that cut off vertex i of
            # this face to continue piecing together the components.
            arcCounts = []
            for i in range(3):
                ip = (i + 1) % 3
                im = (i - 1) % 3

                # Combinatorially, if we take away all the arcs that go into
                # the edge opposite vertex i, then all the remaining points
                # must be paired together to give arcs cutting off vertex i.
                arcCounts.append(
                        ( weights[ip] + weights[im] - weights[i] ) // 2 )
                for j in range( arcCounts[i] ):
                    # Consider the jth closest point to vertex i along the
                    # edge with index edgeInds[ip], and the jth closest point
                    # to vertex i along the edge with index edgeInds[im].
                    # These two points are connected by an arc.
                    if startVerts[ip] == verts[i]:
                        pp = ( edgeInds[ip], j )
                    else:
                        pp = ( edgeInds[ip], weights[ip] - 1 - j )
                    if startVerts[im] == verts[i]:
                        pm = ( edgeInds[im], j )
                    else:
                        pm = ( edgeInds[im], weights[im] - 1 - j )
                    self._union( pp, pm )

                    # Record the fact that the points pp and pm are adjacent,
                    # including information about which vertex gets cut off
                    # (so that we can count switches later on).
                    self._adjPoints[pp].append(pm)
                    self._adjPoints[pm].append(pp)
                    if startVerts[ip] == verts[i]:
                        self._adjCutVerts[pp].append(0)
                    else:
                        self._adjCutVerts[pp].append(1)
                    if startVerts[im] == verts[i]:
                        self._adjCutVerts[pm].append(0)
                    else:
                        self._adjCutVerts[pm].append(1)

    def _findSwitches(self):
        for e in self._tri.edges():
            edgeInd = e.index()
            for i in range( self._weights[edgeInd] ):
                p = ( edgeInd, i )
                cutVerts = self._adjCutVerts[p]
                self._switch[p] = (
                        self._adjCutVerts[p][0] != self._adjCutVerts[p][1] )

    def triangulation(self):
        """
        Returns a copy of the triangulation whose boundary contains this
        normal curve.
        """
        return Triangulation3( self._tri )

    def countComponents(self):
        """
        Returns the number of components of this normal curve.
        """
        return len( self._roots )

    def traverseComponent( self, index ):
        """
        Iterates through the intersection points of the requested component
        of this normal curve.
        """
        startPt = self._roots[index]
        currentPt, nextPt = startPt, self._adjPoints[startPt][0]
        yield currentPt
        while nextPt != startPt:
            adjNext = self._adjPoints[nextPt]
            indCurrent = adjNext.index(currentPt)
            currentPt, nextPt = nextPt, adjNext[ 1-indCurrent ]
            yield currentPt

    def countSwitchesInComponent( self, index ):
        """
        Returns the number of switches in the requested component of this
        normal curve.
        """
        total = 0
        for p in self.traverseComponent(index):
            if self._switch[p]:
                total += 1
        return total

    def _flipWeight( self, e ):
        # Compute the new weight after flipping the edge e.
        # Start by counting "non-corner arcs" that meet e, as these will also
        # meet the new edge.
        newWeight = 0
        ind = e.index()
        wt = self._weights[ind]
        for i in range(wt):
            if self._switch[ (ind, i) ]:
                # Switch implies non-corner.
                newWeight += 1

        # Now count "corner arcs" that will meet the new edge.
        emb = [ e.embedding(0), e.embedding( e.degree() - 1 ) ]
        for i in range(2):
            tet = emb[i].tetrahedron()
            verts = emb[i].vertices()
            otherArcs = 0
            for j in range(2):
                otherInd = tet.edge( verts[j], verts[2+i] ).index()
                otherArcs += self._weights[otherInd]
            otherArcs = ( otherArcs - wt ) // 2
            newWeight += otherArcs

        return newWeight

    def layerOn( self, e ):
        """
        Layer a new tetrahedron across the edge e.

        Pre-condition:
        --> e is a boundary edge of self.triangulation().
        """
        temp = []
        for edge in self._tri.edges():
            if edge == e:
                continue
            emb = edge.embedding(0)
            temp.append(
                    ( emb.tetrahedron(), emb.vertices(),
                        self._weights[ edge.index() ] ) )
        flipWeight = self._flipWeight(e)
        #TODO TEST
        print( "Flip weight: {}.".format(flipWeight) )

        # Layer on new tetrahedron.
        newTet = self._tri.layerOn(e)
        temp.append( ( newTet, Perm4(2,3,0,1), flipWeight ) )

        # Use temp to compute new weight coordinates.
        newWeights = [0] * self._tri.countEdges()
        for tet, verts, wt in temp:
            edgeInd = tet.edge( verts[0], verts[1] ).index()
            newWeights[edgeInd] = wt
        self._processWeights(newWeights)

    def flipEdge( self, e ):
        """
        Flip the edge e by either removing a layered tetrahedron, or layering
        a new tetrahedron across this edge.

        Pre-condition:
        --> e is a boundary edge of self.triangulation().
        """
        #TODO
        pass

    #TODO What's the best way to provide user access to arc coordinates?
    #TODO
    pass


# Test code.
if __name__ == "__main__":
    initTri = Triangulation3.fromIsoSig( "eHbecadjk" )
    msg = "{} components. Switches: {}."
    testCases = [
            (0,0,0,2,2,2,3,1,1),    # 1-component
            (1,0,1,2,3,2,3,2,1),    # 2-component
            (3,5,4,2,3,5,4,5,3) ]   # 3-component
    for w in testCases:
        tri = Triangulation3(initTri)
        curve = NormalCurve( tri, w )
        comps = curve.countComponents()
        print( msg.format( comps,
            [ curve.countSwitchesInComponent(i) for i in range(comps) ] ) )
        curve.layerOn( tri.edge(0) )
        print( msg.format( comps,
            [ curve.countSwitchesInComponent(i) for i in range(comps) ] ) )
        print()
