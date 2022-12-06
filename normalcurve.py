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
        #TODO Do we really need arc coordinates?
        # As a by-product of this, we will also:
        # - compute coordinates for this normal curve in terms of its arcs;
        #   and
        # - assign an index to each component of this normal curve.
        self._arcCoords = [None] * self._tri.countTriangles()
        self._roots = list( self._points )
        # This is also a good opportunity to record useful information about
        # how each component of this normal curve traverses the boundary of
        # our triangulation.
        self._adjPoints = { p: list() for p in self._points }
        self._adjCutVerts = { p: list() for p in self._points }
        self._adjOppEdgeInds = { p: list() for p in self._points }
        self._countComponentsImpl()

        # Find where all the switches happen.
        self._switch = { p: False for p in self._points }
        self._findSwitches()

        # Cache the following whenever they are first computed:
        # - Resolvable components.
        # - Lengths of components.
        self._resolveEdges = [None] * self.countComponents()
        self._lengths = [None] * self.countComponents()

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
                    self._adjOppEdgeInds[pp].append( edgeInds[i] )
                    self._adjOppEdgeInds[pm].append( edgeInds[i] )
                    if startVerts[ip] == verts[i]:
                        self._adjCutVerts[pp].append(0)
                    else:
                        self._adjCutVerts[pp].append(1)
                    if startVerts[im] == verts[i]:
                        self._adjCutVerts[pm].append(0)
                    else:
                        self._adjCutVerts[pm].append(1)
            #TODO Cache arcCounts (as a tuple) in self._arcCoords, or just
            #   get rid of self._arcCoords entirely?

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

    def weights(self):
        """
        Returns the tuple of edge weights that define this normal curve.
        """
        return self._weights

    def countComponents(self):
        """
        Returns the number of components of this normal curve.
        """
        return len( self._roots )

    def _traverseCurve( self, startPt, d ):
        # Either traverse in "direction 0" or "direction 1", depending on
        # whether d is 0 or 1.
        currentPt, nextPt, nextOppEdgeInd = ( startPt,
                self._adjPoints[startPt][d],
                self._adjOppEdgeInds[startPt][d] )
        yield currentPt, nextPt, nextOppEdgeInd
        while nextPt != startPt:
            adjNext = self._adjPoints[nextPt]
            indCurrent = adjNext.index(currentPt)
            currentPt, nextPt, nextOppEdgeInd = ( nextPt,
                    adjNext[ 1-indCurrent ],
                    self._adjOppEdgeInds[nextPt][ 1-indCurrent ] )
            yield currentPt, nextPt, nextOppEdgeInd

    def _traverseComponentImpl( self, index ):
        for t in self._traverseCurve( self._roots[index], 0 ):
            yield t

    def traverseComponent( self, index ):
        """
        Iterates through the intersection points of the requested component
        of this normal curve.
        """
        for p, _, _ in self._traverseComponentImpl(index):
            yield p

    #TODO Cache?
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

    def componentLength( self, index ):
        """
        Returns the length of the requested component of this normal curve.
        """
        if self._lengths[index] is None:
            self._lengths[index] = 0
            for p, _, _ in self._traverseComponentImpl(index):
                self._lengths[index] += 1
        return self._lengths[index]

    def recogniseResolvable( self, index ):
        """
        Recognise the edge indices to which we can resolve the requested
        component of this normal curve.
        """
        if self._resolveEdges[index] is None:
            switches = 0
            resEdgeInds = []
            for info in self._traverseComponentImpl(index):
                currentPt, nextPt, nextOppEdgeInd = info
                if self._switch[currentPt]:
                    switches += 1
                    if self._switch[nextPt]:
                        resEdgeInds.append(nextOppEdgeInd)
            if switches != 2:
                resEdgeInds = []
            self._resolveEdges[index] = tuple(resEdgeInds)
        return self._resolveEdges[index]

    def _bubblePoints( self, e, i ):
        # Find all intersection points that form a bubble around vertex i of
        # edge e.
        if self._weights[ e.index() ] == 0:
            return set()
        bubble = set()
        if i == 0:
            startPt = ( e.index(), 0 )
        else:
            startPt = ( e.index(), self._weights[ e.index() ] - 1 )
        bubble.add(startPt)

        # Traverse in both directions until we turn away from vertex i.
        for d in range(2):
            if self._adjCutVerts[startPt][d] == 1 - i:
                # Curve immediately turns away from vertex i.
                continue

            # Curve initially turns towards vertex i, and it will continue
            # doing so until we reach a switch.
            for _, p, _ in self._traverseCurve( startPt, d ):
                bubble.add(p)
                if self._switch[p]:
                    break

        # If bubble only contains startPt, then the curve turns away from
        # vertex i in both directions.
        if len(bubble) == 1:
            return set()
        else:
            return bubble

    def _isotopeOffEdge( self, e, i ):
        # Try to isotope the arc that goes around vertex i of edge e, and
        # return True if and only if this was possible.
        # Start by ruling out special cases:
        # - We have no bubble around vertex i because the curve turns away
        #   from vertex i in both directions.
        # - The bubble meets every intersection point of a component, which
        #   means that we have either: (1) instead of isotoping, we should
        #   resolve the curve; or (2) the curve is a vertex link.
        bubble = self._bubblePoints( e, i )
        lb = len(bubble)
        if lb == 0:
            return False
        edgeInd = e.index()
        if i == 0:
            pt = ( edgeInd, 0 )
        else:
            pt = ( edgeInd, self._weights[edgeInd] - 1 )
        if lb == self.componentLength(
                self._roots.index( self._find(pt) ) ):
            return False

        # Work out how weights would change after isotopy.
        newWeights = list( self._weights )
        for edge in self._tri.edges():
            if not edge.isBoundary():
                continue
            edgeInd = edge.index()

            # How do the weights change on this edge?
            oldWeight = self._weights[edgeInd]
            if oldWeight == 1:
                if ( edgeInd, 0 ) not in bubble:
                    newWeights[edgeInd] += 2
            else:
                if ( edgeInd, 0 ) in bubble:
                    newWeights[edgeInd] -= 1
                else:
                    newWeights[edgeInd] += 1
                if ( edgeInd, oldWeight - 1 ) in bubble:
                    newWeights[edgeInd] -= 1
                else:
                    newWeights[edgeInd] += 1
        self._processWeights(newWeights)
        return True

    def _clearEdge( self, e ):
        # Assuming that a component of this normal curve resolves to the edge
        # e, isotope everything off e to make resolving possible.
        tet = e.embedding(0).tetrahedron()
        verts = e.embedding(0).vertices()
        while self._weights[ e.index() ] > 0:
            self._isotopeOffEdge( e, 0 )
            e = tet.edge( verts[0], verts[1] )

    def resolveComponent( self, index, check=True, perform=True ):
        """
        Checks the eligibility of and/or performs the operation of resolving
        the requested component of this normal curve.

        If this routine is asked to both check and perform, then it only
        performs the operation if the check shows it is legal.

        If check is True, then this routine returns True if and only if
        resolving the requested component is legal. If check is False, then
        this routine always returns True.

        Pre-condition:
        --> If this routine is asked to perform the operation without first
            running the check, then it must be known in advance that
            resolving the requested component is legal.
        """
        #TODO
        pass

    #TODO What's the best way to provide user access to arc coordinates?
    #TODO
    pass


# Test code.
if __name__ == "__main__":
    initTri = Triangulation3.fromIsoSig( "eHbecadjk" )
    compsMsg = "{} component(s):"
    edgeMsg = "After layering on edge {}:"
    subMsg = "    Switches: {}. Recognise resolvable: {}."
    #TODO
    msg = "{} components. Switches: {}."
    testCases = [
            ( (0,0,0,2,2,2,3,1,1), 0 ),     # 1-component
            ( (1,0,1,2,3,2,3,2,1), 4 ),     # 2-component
            ( (3,5,4,2,3,5,4,5,3), 0 ),     # 3-component
            ( (1,1,0,1,0,0,0,0,0), 3 ),     # 1-component
            ( (5,5,4,0,5,5,5,4,0), 0 ) ]    # 3-component
    for w, e in testCases:
        tri = Triangulation3(initTri)
        curve = NormalCurve( tri, w )

        # Test components.
        comps = curve.countComponents()
        print( compsMsg.format(comps) )
        for c in range(comps):
            print( subMsg.format( curve.countSwitchesInComponent(c),
                curve.recogniseResolvable(c) ) )

        # Test layering.
        curve.layerOn( tri.edge(e) )
        print( edgeMsg.format(e) )
        for c in range(comps):
            print( subMsg.format( curve.countSwitchesInComponent(c),
                curve.recogniseResolvable(c) ) )

        # Test bubble/isotopy.
        tri = Triangulation3(initTri)
        curve = NormalCurve( tri, w )
        iso = curve._isotopeOffEdge( tri.edge(0), 0 )
        print( "Isotopy: {}.".format(iso) )
        if iso:
            print( curve.weights() )
            comps = curve.countComponents()
            for c in range(comps):
                print( subMsg.format( curve.countSwitchesInComponent(c),
                    curve.recogniseResolvable(c) ) )
        #print( "Bubble: {}.".format(
        #    curve._bubblePoints( tri.edge(0), 0 ) ) )
        print()

    # Test clearing edge 7 of testCases[4].
    print( "Clear edge." )
    tri = Triangulation3(initTri)
    curve = NormalCurve( tri, testCases[4][0] )
    curve._clearEdge( tri.edge(7) )
    print( curve.weights() )
    print()
    #TODO
