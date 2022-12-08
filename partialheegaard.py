"""
A class to represent a partial triangulation of a Heegaard splitting.
"""
from regina import *


class PartialHeegaardSplitting:
    """
    A partial triangulation of a Heegaard splitting, with Heegaard curves
    represented as a union of a normal curve (possibly with many components)
    and a collection of boundary edges.

    Each Heegaard curve given by a component of the normal curve is
    *unresolved*, and each Heegaard curve given by a boundary edge is
    *resolved*.
    """
    def __init__( self, tri, weights, resolvedEdgeIndices=set() ):
        """
        Create a partial triangulation with a boundary curve, described by
        the edge weights of a normal curve together (optionally) with some
        boundary edges.
        """
        # Check that we have the correct number of weights.
        if len(weights) != tri.countEdges():
            raise ValueError( "Wrong number of edge weights!" )

        self._tri = tri
        self._resolvedEdges = set(resolvedEdgeIndices)
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
        self._resolvableEdges = [None] * self.countUnresolved()
        self._lengths = [None] * self.countUnresolved()

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

    def getResolvedEdgeIndices(self):
        """
        Returns a copy of the set of resolved edge indices.
        """
        return set( self._resolvedEdges )

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
            temp.append( (
                emb.tetrahedron(),
                emb.vertices(),
                self._weights[ edge.index() ],
                edge.index() ) )
        flipWeight = self._flipWeight(e)

        # Layer on new tetrahedron.
        newTet = self._tri.layerOn(e)
        temp.append( ( newTet, Perm4(2,3,0,1), flipWeight, None ) )

        # Use temp to compute new weight coordinates.
        newWeights = [0] * self._tri.countEdges()
        newResEdges = set()
        for tet, verts, wt, oldInd in temp:
            newInd = tet.edge( verts[0], verts[1] ).index()
            newWeights[newInd] = wt
            if oldInd in self._resolvedEdges:
                newResEdges.add(newInd)
        self._processWeights(newWeights)
        self._resolvedEdges = newResEdges

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
            for _ in self._traverseComponentImpl(index):
                self._lengths[index] += 1
        return self._lengths[index]

    def recogniseResolvable( self, index ):
        """
        Recognise the edge indices to which we can resolve the requested
        component of this normal curve.
        """
        if self._resolvableEdges[index] is None:
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
            self._resolvableEdges[index] = tuple(resEdgeInds)
        return self._resolvableEdges[index]

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

    def resolveComponent( self, index ):
        """
        Checks whether it is possible to resolve the requested component, and
        if so resolves this component.
        """
        resEdgeInds = self.recogniseResolvable(index)
        if not resEdgeInds:
            return False

        # First compute how edge weights would reduce as a result of
        # resolving the requested component.
        weightChanges = [0] * self._tri.countEdges()
        for p, _, _ in self._traverseComponentImpl(index):
            weightChanges[ p[0] ] -= 1

        # Clear the edge, and then resolve the requested component by simply
        # deleting its contributions to edge weight and adding its index to
        # the set of resolved edge indices.
        edgeInd = resEdgeInds[0]
        self._clearEdge( self._tri.edge(edgeInd) )
        self._resolvedEdges.add(edgeInd)
        newWeights = [ self._weights[i] + weightChanges[i] for
                i in range( self._tri.countEdges() ) ]
        self._processWeights(newWeights)
        return True

    #TODO What's the best way to provide user access to arc coordinates?
    #TODO
    pass


# Test code.
if __name__ == "__main__":
    initTri = Triangulation3.fromIsoSig( "eHbecadjk" )
    compsMsg = "{} component(s):"
    edgeMsg = "After layering on edge {}:"
    subMsg = "    Switches: {}. Recognise resolvable: {}."
    msg = "{} components. Switches: {}."
    testCases = [
            ( (0,0,0,2,2,2,3,1,1), 0 ),     # 1-component
            ( (1,0,1,2,3,2,3,2,1), 4 ),     # 2-component
            ( (3,5,4,2,3,5,4,5,3), 0 ),     # 3-component
            ( (1,1,0,1,0,0,0,0,0), 3 ),     # 1-component
            ( (5,5,4,0,5,5,5,4,0), 0 ) ]    # 3-component

    # Basic operation tests.
    for w, e in testCases:
        tri = Triangulation3(initTri)
        phs = PartialHeegaardSplitting( tri, w )

        # Test components.
        comps = phs.countUnresolved()
        print( compsMsg.format(comps) )
        for c in range(comps):
            print( subMsg.format( phs.countSwitchesInComponent(c),
                phs.recogniseResolvable(c) ) )

        # Test layering.
        phs.layerOn( tri.edge(e) )
        print( edgeMsg.format(e) )
        for c in range(comps):
            print( subMsg.format( phs.countSwitchesInComponent(c),
                phs.recogniseResolvable(c) ) )

        # Test bubble/isotopy.
        tri = Triangulation3(initTri)
        phs = PartialHeegaardSplitting( tri, w )
        iso = phs._isotopeOffEdge( tri.edge(0), 0 )
        print( "Isotopy: {}.".format(iso) )
        if iso:
            print( phs.weights() )
            comps = phs.countUnresolved()
            for c in range(comps):
                print( subMsg.format( phs.countSwitchesInComponent(c),
                    phs.recogniseResolvable(c) ) )
        #print( "Bubble: {}.".format(
        #    phs._bubblePoints( tri.edge(0), 0 ) ) )
        print()

    # Test clearing edge 7 of testCases[4].
    print( "Clear edge." )
    tri = Triangulation3(initTri)
    phs = PartialHeegaardSplitting( tri, testCases[4][0] )
    phs._clearEdge( tri.edge(7) )
    print( phs.weights() )
    print()

    # Test resolving (to edge 7 of testCases[4]).
    print( "Resolve components." )
    tri = Triangulation3(initTri)
    phs = PartialHeegaardSplitting( tri, testCases[4][0] )
    for i in range(3):
        res = phs.resolveComponent(0)
        print( "Resolve component {}: {}.".format( i, res ) )
        if res:
            print( "Resolved edge indices: {}.".format(
                phs.getResolvedEdgeIndices() ) )
            print( "Weights: {}.".format( phs.weights() ) )
            break
    print()
    #TODO
