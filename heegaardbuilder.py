"""
A class for building triangulations from a Heegaard diagram.
"""
from regina import *


class HeegaardBuilder:
    #TODO Proofread documentation.
    #TODO So far, the code doesn't really check for valid inputs.
    #TODO For checking handlebody, don't do a conclusive test. Just use
    #   knowsHandlebody() to weed out obviously bad cases. This obviously
    #   needs to be clearly documented.
    #TODO We need to store a boolean that says whether we have closed up (or
    #   started closing up?) the triangulation, because otherwise the
    #   handlebody checks may fail when cloning.
    #TODO Should we perhaps implement *simultaneous* resolving once we know
    #   that all components as "nice"?
    """
    Implements the various subroutines needed to execute an algorithm for
    constructing a 3-manifold triangulation from a set of Heegaard curves on
    a given handlebody.

    The handlebody must be represented using a one-vertex triangulation
    (preferably using a layered construction, for the reason outlined below).
    The Heegaard curves are represented as a union of:
    --> a normal curve (typically with many components) in the boundary
        surface; and
    --> a collection of boundary edges.
    Each Heegaard curve given by a component of the normal curve is
    *unresolved*, and each Heegaard curve given by a boundary edge is
    *resolved*.

    If, in addition, the initial handlebody uses a layered construction---as
    described in Jaco and Rubinstein's unpublished 2006 paper titled
    "Layered-triangulations of 3-manifolds"---then the method of construction
    implemented by this class guarantees that the cutwidth of the
    triangulation is always bounded above by 4g-2, where g is the genus of
    the initial handlebody. Note that because cutwidth is an upper bound for
    treewidth, this means that the treewidth of the triangulation is also
    always bounded above by 4g-2.
    """
    def __init__( self, tri, weights, resolvedEdgeIndices=set() ):
        """
        Initialises a HeegaardBuilder object that provides the subroutines
        needed to turn tri into a one-vertex triangulation of the closed
        3-manifold corresponding to the Heegaard curves represented by the
        given weights and resolvedEdgeIndices.

        The given triangulation tri should be a one-vertex triangulation of
        a handlebody. As described in the class notes, it is preferable for
        tri to use a layered construction, as this will guarantee an upper
        bound on the cutwidth of the triangulations that this class
        constructs.

        The given weights should describe a normal curve c in the boundary
        surface of tri as follows:
        --> We should have len(weights) == tri.countEdges().
        --> For each i such that e = tri.edge(i) is a boundary edge,
            weights[i] should give the edge weight of c at e.
        The components of this curve c correspond to unresolved Heegaard
        curves.
        
        The resolvedEdgeIndices give indices of boundary edges of tri; for
        each i in resolvedEdgeIndices, the edge tri.edge(i) gives a resolved
        Heegaard curve.

        Since our set of Heegaard curves must be disjoint, note that:
        --> For each i in resolvedEdgeIndices, weights[i] must be 0.
        --> Since the edges corresponding to resolved curves all meet at the
            vertex of tri, we insist that they all meet tangentially (rather
            than transversely).
        """
        # Check that we have the correct number of weights.
        if len(weights) != tri.countEdges():
            raise ValueError( "Wrong number of edge weights!" )

        self._tri = Triangulation3(tri)
        self._resolvedEdges = set(resolvedEdgeIndices)
        self._processWeights(weights)

    def clone(self):
        """
        Creates and returns a clone of this HeegaardBuilder instance.

        Note that the underlying triangulation of the cloned instance will be
        a *copy* of the underlying triangulation of this instance.
        """
        return HeegaardBuilder(
                Triangulation3( self._tri ),
                self._weights,
                set(self._resolvedEdges) )

    def triangulation(self):
        """
        Returns a copy of the underlying triangulation.
        """
        return Triangulation3( self._tri )

    def _processWeights( self, weights ):
        """
        Update internal data to reflect a new set of edge weights.
        """
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
        #   and
        # - assign an index to each component of this normal curve.
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
        # - Components that are isotopic to off-diagonals.
        # - Lengths of components.
        # - Number of switches.
        self._resolvableEdges = [None] * self.countUnresolved()
        self._offDiagEdges = [None] * self.countUnresolved()
        self._lengths = [None] * self.countUnresolved()
        self._switchCounts = [None] * self.countUnresolved()

    def _checkMatching(self):
        """
        Do the underlying edge weights satisfy the matching constraints?
        """
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
        """
        Initialise intersection points for the underlying normal curve, so
        that we can compute components using union-find.
        """
        # Represent each intersection point using a pair consisting of the
        # edge index and the order of the point along the edge.
        for edgeInd in range( self._tri.countEdges() ):
            if not self._tri.edge(edgeInd).isBoundary():
                continue
            for p in range( self._weights[edgeInd] ):
                self._points.add( ( edgeInd, p ) )

    def _find( self, p ):
        """
        Find the current root for the given intersection point p.
        """
        while self._parent[p] != p:
            # Use path-halving: replace the current parent of p with its
            # current grandparent, and then go to this grandparent.
            self._parent[p] = self._parent[ self._parent[p] ]
            p = self._parent[p]
        return p

    def _union( self, p, q ):
        """
        Take the union of the components containing the given intersection
        points p and q.
        """
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
        """
        Union-find implementation for counting the components of the
        underlying normal curve.
        """
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

    def _findSwitches(self):
        """
        Find switches in the underlying normal curve.
        """
        for e in self._tri.edges():
            edgeInd = e.index()
            for i in range( self._weights[edgeInd] ):
                p = ( edgeInd, i )
                cutVerts = self._adjCutVerts[p]
                self._switch[p] = (
                        self._adjCutVerts[p][0] != self._adjCutVerts[p][1] )

    def countUnresolved(self):
        """
        Returns the number of unresolved components of the underlying curve.
        """
        return len( self._roots )

    def countResolved(self):
        """
        Returns the number of resolved components of the underlying curve to
        which we haven't yet attached a disc.
        """
        return len( self._resolvedEdges )

    def countComponents(self):
        """
        Returns the number of components of the underlying curve to which we
        haven't yet attached a disc.
        """
        return self.countUnresolved() + self.countResolved()

    def _traverseCurve( self, startPt, d ):
        """
        Traverse the unresolved Heegaard curve containing the given startPt.

        Traversal begins at startPt, and proceeds in the direction towards
        self._adjPoints[startPt][d].

        At each intersection point P along the curve that we traverse, this
        routine yields a 4-tuple consisting of:
        (0) The point that occurs before P in the direction of traversal.
        (1) The point P.
        (2) The point Q that occurs after P in the direction of traversal.
        (3) The index of the edge parallel to the normal arc that connects
            the points P and Q.
        """
        # Either traverse in "direction 0" or "direction 1", depending on
        # whether d is 0 or 1.
        prevPt, currentPt, nextPt, nextOppEdgeInd = (
                self._adjPoints[startPt][1 - d],
                startPt,
                self._adjPoints[startPt][d],
                self._adjOppEdgeInds[startPt][d] )
        yield prevPt, currentPt, nextPt, nextOppEdgeInd
        while nextPt != startPt:
            indNext = 1 - self._adjPoints[nextPt].index(currentPt)
            prevPt, currentPt, nextPt, nextOppEdgeInd = (
                    currentPt,
                    nextPt,
                    self._adjPoints[nextPt][indNext],
                    self._adjOppEdgeInds[nextPt][indNext] )
            yield prevPt, currentPt, nextPt, nextOppEdgeInd

    def _traverseComponent( self, index ):
        """
        Traverses the requested unresolved component.

        At each intersection point P that we encounter during this traveral,
        this routine yields a 4-tuple consisting of:
        (0) The point that occurs before P in the direction of traversal.
        (1) The point P.
        (2) The point Q that occurs after P in the direction of traversal.
        (3) The index of the edge parallel to the normal arc that connects
            the points P and Q.
        """
        for t in self._traverseCurve( self._roots[index], 0 ):
            yield t

    def _countSwitchesInComponent( self, index ):
        """
        Returns the number of switches in the requested unresolved component
        of the underlying curve.

        The first call to this routine caches the result, so subsequent calls
        with the same index will not need to recompute the answer.
        """
        if self._switchCounts[index] is None:
            self._switchCounts[index] = 0
            for _, p, _, _ in self._traverseComponent(index):
                if self._switch[p]:
                    self._switchCounts[index] += 1
        return self._switchCounts[index]

    def _flipWeightChange( self, e ):
        """
        Compute how the weight changes after flipping the given edge e.
        """
        # Start by counting "corner arcs" that meet e; we will lose the
        # corresponding intersection points after flipping e.
        weightChange = 0
        ind = e.index()
        wt = self._weights[ind]
        for i in range(wt):
            if not self._switch[ (ind, i) ]:
                # No switch means we have a "corner arc".
                weightChange -= 1

        # Now count "corner arcs" that will meet the new edge.
        emb = [ e.embedding(0), e.embedding( e.degree() - 1 ) ]
        for i in range(2):
            tet = emb[i].tetrahedron()
            ver = emb[i].vertices()
            otherArcs = 0
            for j in range(2):
                otherInd = tet.edge( ver[j], ver[2+i] ).index()
                otherArcs += self._weights[otherInd]
            otherArcs = ( otherArcs - wt ) // 2
            weightChange += otherArcs

        return weightChange

    def _flipWeight( self, e ):
        """
        Compute the weight of the new edge that would result from flipping
        the given edge e.
        """
        return self._weights[ e.index() ] + self._flipWeightChange(e)

    def _recogniseLayer( self, e ):
        """
        Checks whether e is the new boundary edge that results from layering
        a tetrahedron t across an edge ee, and if so returns details of this
        layering; otherwise, returns None.

        Specifically, if we do have such a layering, then this routine
        returns a triple consisting of:
        (0) the tetrahedron t,
        (1) another tetrahedron tt (distinct from t), and
        (3) an edge number n,
        such that tt.edge(n) == ee.
        """
        # The edge e must be a boundary edge of degree 1.
        if not e.isBoundary() or e.degree() != 1:
            return None

        # The opposite edge must be internal, and the opposite triangles must
        # be internal and distinct.
        tet = e.embedding(0).tetrahedron()
        ver = e.embedding(0).vertices()
        oppEdge = tet.edge( ver[2], ver[3] )
        if oppEdge.isBoundary():
            return None
        oppTri = [ tet.triangle( ver[i] ) for i in range(2) ]
        if ( oppTri[0].isBoundary() or oppTri[1].isBoundary() or
                oppTri[0] == oppTri[1] ):
            return None

        # Find a (tetrahedron, edge number) pair corresponding to oppEdge.
        for emb in oppEdge.embeddings():
            if emb.tetrahedron() == tet:
                continue
            oppTet, oppEdgeNum = emb.tetrahedron(), emb.edge()
            break
        return ( tet, oppTet, oppEdgeNum )

    def _flipEdgeImpl( self, e, aux=None ):
        """
        Flips the given boundary edge e and returns the new boundary edge
        that results from this flip.

        This routine automatically updates edge weights and resolved edge
        indices. In the special case where e is itself one of the resolved
        edges, this routine simply removes e from the set of resolved edges
        (the intention is that the only time we should ever flip a resolved
        edge is when we want to attach a disc by folding along this edge).

        If aux is None (the default), then this routine has no additional
        functionality.

        However, if aux is not None, then it must be a 2-element list
        consisting of the following items:
        (0) A list F of 2-tuples, each of which consists of a tetrahedron t
            and an integer n such that t.edge(n) is a boundary edge of the
            underlying triangulation.
        (1) A dictionary D that maps tetrahedron indices i to the set of
            locations in the list F that reference tetrahedron i of the
            underlying triangulation.
        In this case, this routine provides the additional functionality of
        adjusting these two items subject to the following conditions:
        --> Each pair (t, n) in F represents the same edge as before, but we
            may have changed the tetrahedron t to ensure that it is still a
            tetrahedron that exists in the underlying triangulation (this is
            necessary because flipping the edge e sometimes involves removing
            a tetrahedron).
        --> The dictionary D is adjusted so that the correspondence between
            tetrahedron indices and elements of F, as described above, still
            holds.

        Pre-condition:
        --> e is a boundary edge of the underlying triangulation.
        --> If aux is not None, then it must be a 2-element list whose
            elements satisfy the description given above.
        """
        removeInfo = self._recogniseLayer(e)
        if removeInfo is None:
            doomed = None
        else:
            doomed, newTet, newEdgeNum = removeInfo

        # If we are going to perform this flip by removing a "doomed"
        # tetrahedron, and if we are given some aux data that needs to be
        # adjusted to account for this, then perform these adjustments now.
        if ( doomed is not None ) and ( aux is not None ):
            flip, refTetInds = aux

            # Make sure we stop referencing the doomed tetrahedron.
            if doomed.index() in refTetInds:
                for loc in refTetInds[ doomed.index() ]:
                    if loc >= len(flip):
                        continue

                    # Find a new (tetrahedron, edge number) pair that avoids
                    # the doomed tetrahedron.
                    oldRefTet, oldRefNum = flip[loc]
                    refEdge = oldRefTet.edge(oldRefNum)
                    for emb in refEdge.embeddings():
                        newRefTet = emb.tetrahedron()
                        if newRefTet != doomed:
                            newRefNum = emb.edge()
                            break

                    # Update both flip and refTetInds.
                    flip[loc] = ( newRefTet, newRefNum )
                    if newRefTet.index() in refTetInds:
                        refTetInds[ newRefTet.index() ].add(loc)
                    else:
                        refTetInds[ newRefTet.index() ] = {loc}

            # Account for the renumbering of tetrahedra.
            newRefTetInds = dict()
            for tetInd in refTetInds:
                if tetInd > doomed.index():
                    newRefTetInds[ tetInd - 1 ] = refTetInds[tetInd]
                elif tetInd < doomed.index():
                    newRefTetInds[tetInd] = refTetInds[tetInd]
            refTetInds.clear()
            refTetInds.update(newRefTetInds)

        # Before we perform the flip, we need to record the information
        # needed to update the edge weights and the resolved edge indices.
        temp = []
        for edge in self._tri.edges():
            if edge == e:
                # Since e will no longer be a boundary edge after flipping,
                # we don't want to remember anything about e.
                continue

            # Get a (tetrahedron, edge number) pair corresponding to edge,
            # taking care that we do not choose the doomed tetrahedron.
            for emb in edge.embeddings():
                tet = emb.tetrahedron()
                if tet != doomed:
                    edgeNum = emb.edge()
                    break
            temp.append( (
                tet,
                edgeNum,
                self._weights[ edge.index() ],
                edge.index() ) )
        if self.countUnresolved() == 0:
            flipWeight = 0
        else:
            flipWeight = self._flipWeight(e)

        # We are now ready to perform the flip. Do this by either removing a
        # doomed tetrahedron, or layering on a new tetrahedron.
        if doomed is None:
            newTet = self._tri.layerOn(e)
            newEdgeNum = 5
        else:
            self._tri.removeTetrahedron(doomed)

        # Before finishing up, we need to update the underlying curves.
        temp.append( ( newTet, newEdgeNum, flipWeight, None ) )
        newWeights = [0] * self._tri.countEdges()
        newResEdges = set()
        for tet, edgeNum, wt, oldInd in temp:
            newInd = tet.edge(edgeNum).index()
            newWeights[newInd] = wt
            if oldInd in self._resolvedEdges:
                newResEdges.add(newInd)
        self._resolvedEdges = newResEdges
        self._processWeights(newWeights)
        return newTet.edge(newEdgeNum)

    def flipEdge( self, e ):
        """
        If the given edge e is a boundary edge that does not form one of the
        resolved components of the underlying curve, then flips e and returns
        the new boundary edge that results from this flip; otherwise, simply
        returns None.

        If we have at least unresolved curve, then this routine also
        automatically updates the edge weights accordingly.

        If e happens to be the new boundary edge that results from layering a
        tetrahedron t, then we flip e by simply removing t. Otherwise, we
        flip e by layering a new tetrahedron across e.
        """
        if ( not e.isBoundary() ) or ( e.index() in self._resolvedEdges ):
            return None
        else:
            return self._flipEdgeImpl(e)

    def _componentLength( self, index ):
        """
        Returns the length of the requested unresolved component of the
        underlying curve.

        The first call to this routine caches the result, so subsequent calls
        with the same index will not need to recompute the answer.
        """
        if self._lengths[index] is None:
            self._lengths[index] = 0
            for _ in self._traverseComponent(index):
                self._lengths[index] += 1
        return self._lengths[index]

    def recogniseResolvable( self, index ):
        """
        Recognise the edge indices to which we can resolve the requested
        unresolved component of the underlying curve.

        The returned object will be a tuple of edge indices; the length of
        this tuple will always be no greater than two.

        The first call to this routine caches the result, so subsequent calls
        with the same index will not need to recompute the answer.
        """
        if self._resolvableEdges[index] is None:
            switches = 0
            resEdgeInds = []
            for info in self._traverseComponent(index):
                _, currentPt, nextPt, nextOppEdgeInd = info
                if self._switch[currentPt]:
                    switches += 1
                    if self._switch[nextPt]:
                        resEdgeInds.append(nextOppEdgeInd)
            if switches == 2:
                self._resolvableEdges[index] = tuple(resEdgeInds)
            else:
                self._resolvableEdges[index] = tuple()
        return self._resolvableEdges[index]

    def _recogniseOffDiagonal( self, index ):
        """
        Recognise the edge indices across which the requested unresolved
        component forms an off-diagonal.

        The first call to this routine caches the result, so subsequent calls
        with the same index will not need to recompute the answer.
        """
        if self._offDiagEdges[index] is None:
            switches = 0
            crossDiagInds = []
            for info in self._traverseComponent(index):
                prevPt, currentPt, nextPt, nextOppEdgeInd = info
                if self._switch[prevPt]:
                    switches += 1
                    if not self._switch[currentPt] and self._switch[nextPt]:
                        # We are in the middle of a bubble of length 2.
                        crossDiagInds.append( currentPt[0] )
            if switches == 2 and crossDiagInds:
                self._offDiagEdges[index] = tuple(crossDiagInds)
            else:
                self._offDiagEdges[index] = tuple()
        return self._offDiagEdges[index]

    def _bubblePoints( self, e, i ):
        """
        Find all intersection points that form a bubble around vertex i of
        edge e.
        """
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
            for _, _, p, _ in self._traverseCurve( startPt, d ):
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
        """
        Try to isotope the arc that goes around vertex i of edge e, and
        return True if and only if this was possible.
        """
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
        if lb == self._componentLength(
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
        """
        Assuming that a component of the underlying normal curve resolves to
        the edge e, isotope everything off e to make resolving possible.
        """
        tet = e.embedding(0).tetrahedron()
        verts = e.embedding(0).vertices()
        while self._weights[ e.index() ] > 0:
            #TODO This check shouldn't be necessary, but leave this here
            #   until I have implemented checks in more appropriate places.
            if not self._isotopeOffEdge( e, 0 ):
                #TODO Sometimes we can't isotope because we actually have to
                #   resolve another curve first.
                #TODO What about off-diagonals?
                # Print enough info so I can look at this by hand.
                for be in self._tri.edges():
                    if not be.isBoundary():
                        continue
                    print( "Ind: {}. Emb0: {}. Emb1: {}. Wt: {}.".format(
                        be.index(),
                        be.embedding(0),
                        be.embedding( be.degree() - 1 ),
                        self._weights[ be.index() ] ) )
                raise RuntimeError( "Isotoping unexpectedly failed." )
            e = tet.edge( verts[0], verts[1] )

    def resolveComponent( self, index, edgeInd=None ):
        """
        Checks whether it is possible to resolve the requested unresolved
        component, and if so resolves this component.

        This routine takes an additional edgeInd argument, which is handled
        as follows:
        --> If edgeInd is None (the default), then this routine returns True
            if and only if the requested component can be resolved to any
            particular boundary edge of the underlying triangulation. In the
            case where such a resolution is possible, this routine makes an
            arbitrary choice for which one to perform (we may have up to two
            choices).
        --> If edgeInd is not None, then it must be the index an edge of the
            underlying triangulation, and this routine returns True if and
            only if the requested component can be resolved to the specific
            edge corresponding to this given edge index. In the case where
            such a resolution is possible, this routine guarantees to resolve
            to the edge specified by edgeInd (even when another choice of
            edge is available).

        Pre-condition:
        --> If edgeInd is not None, then it must be an integer corresponding
            to the index of an edge of the underlying triangulation.
        """
        # Is the requested resolution possible?
        resEdgeInds = self.recogniseResolvable(index)
        if edgeInd is None:
            if not resEdgeInds:
                return False
            else:
                edgeInd = resEdgeInds[0]
        else:
            if edgeInd not in resEdgeInds:
                return False

        # We can definitely resolve. Before we perform the resolution,
        # compute the subsequent reductions to the edge weights.
        weightChanges = [0] * self._tri.countEdges()
        for _, p, _, _ in self._traverseComponent(index):
            weightChanges[ p[0] ] -= 1

        # Clear the edge, and then resolve the requested component by simply:
        # - deleting its contributions to edge weight, and
        # - adding edgeInd to the set of resolved edge indices.
        self._clearEdge( self._tri.edge(edgeInd) )
        self._resolvedEdges.add(edgeInd)
        newWeights = [ self._weights[i] + weightChanges[i] for
                i in range( self._tri.countEdges() ) ]
        self._processWeights(newWeights)
        return True

    #TODO Rename?
    def _isHeegaard(self):
        """
        Does the underlying curve give a valid Heegaard diagram?
        """
        #TODO Check handlebody, number of curves, and separating?
        pass

    #TODO Make public?
    def _allComponentsNice(self):
        """
        Are all components of the underlying curve "nice" in the sense that
        they are either resolved, resolvable or isotopic to an off-diagonal?
        """
        for i in range( self.countUnresolved() ):
            if ( not self.recogniseResolvable(i) and
                    not self._recogniseOffDiagonal(i) ):
                return False
        return True

    def _isMaximalWeight( self, e ):
        """
        Is the edge e maximal-weight?

        If e is not even a boundary edge, then this routine always returns
        False. Otherwise, this routine returns True if and only if the weight
        of e is maximal among boundary edges.
        """
        return ( self._weights[ e.index() ] == max( self._weights ) )

    def _isReducible( self, e ):
        """
        Is the edge e reducible?
        """
        return ( self._flipWeightChange(e) < 0 )

    def _maxReducibleEdges(self):
        """
        Yields all maximal-weight reducible edges.
        """
        for e in self._tri.edges():
            if ( e.isBoundary() and self._isMaximalWeight(e) and
                    self._isReducible(e) ):
                yield e

    def _resolutions(self):
        """
        Yields all possible ways we can resolve a component of the underlying
        normal curve to a boundary edge of the underlying triangulation.

        In more detail, this routine yields all possible pairs (i, e), where:
        --> i is the index of a resolvable component c; and
        --> e is the index of an edge to which we can resolve c.
        """
        for i in range( self.countUnresolved() ):
            for e in self.recogniseResolvable(i):
                yield (i, e)

    def _maxReduceSwitchEdges(self):
        """
        Yields all maximal-weight edges such that flipping the edge reduces
        the total number of switches.
        """
        for e in self._tri.edges():
            if not e.isBoundary() or not self._isMaximalWeight(e):
                continue

            # Check whether flipping this maximal-weight edge e would reduce
            # the number of switches.
            for i in range( self._weights[ e.index() ] ):
                pt = ( e.index(), i )

                # Look for "alternating" pattern.
                if ( self._switch[pt] and
                        self._switch[ self._adjPoints[pt][0] ] and
                        self._switch[ self._adjPoints[pt][1] ] ):
                    # The current edge e is one of the edges that we're
                    # looking for, so we can yield and then move on to the
                    # next edge.
                    yield e
                    break

    def _resolveAllImpl(self):
        #TODO Why is it safe to resolve everything in series? Isn't it
        #   possible to make a curve unresolvable after clearing an edge?
        """
        Assuming that every component of the underlying curve is nice,
        resolves all these components.

        A component is "nice" if it is either resolved, resolvable, or
        isotopic to an off-diagonal.

        Pre-condition:
        --> Every component of the underlying curve must be nice in the sense
            described above.
        """
        # We assume that every component is nice, and we want to resolve all
        # the components that are either resolvable or isotopic to an
        # off-diagonal. Handle the resolvable components first.
        res = True
        while res:
            res = False
            for i in range( self.countUnresolved() ):
                #TODO We make a choice here. Do I want to allow
                #   experimentation with this?
                if self.resolveComponent(i):
                    res = True
                    break

        # Now handle the off-diagonals.
        diag = True
        while diag:
            diag = False
            for i in range( self.countUnresolved() ):
                diagEdgeInds = self._recogniseOffDiagonal(i)
                if diagEdgeInds:
                    diag = True
                    #TODO We make a choice here. Do I want to allow
                    #   experimentation with this?
                    self.flipEdge(
                            self._tri.edge( diagEdgeInds[0] ) )
                    break
        while self.countUnresolved() > 0:
            #TODO Do we need this test?
            if not self.resolveComponent(0):
                raise RuntimeError( "Unexpected unresolvable component." )
            #TODO
            pass

    #TODO Version of this that resolves until we have multiple choices. There
    #   might be a lot of common code, in which case in might make sense to
    #   have a private _resolveAllImpl() routine.
    #TODO Possibly a more sensible way to minimise repeated code is to have
    #   separate implementations for things like finding maximal-weight
    #   reducible edges.
    def resolveAll(self):
        """
        Resolves all components of the underlying normal curve.
        """
        while not self._allComponentsNice():
            # Try to flip a maximal-weight reducible edge.
            foundMaxRed = False
            for e in self._maxReducibleEdges():
                foundMaxRed = True
                self.flipEdge(e)
                break
            if foundMaxRed:
                continue

            # If there is no maximal-weight reducible edge, then try to
            # resolve a component.
            resolved = False
            for i, e in self._resolutions():
                self.resolveComponent(i, e)
                resolved = True
                break
            if resolved:
                continue

            # If we can't resolve either, then we must be able to reduce the
            # number of switches by flipping a maximal-weight edge.
            foundSwitch = False
            for e in self._maxReduceSwitchEdges():
                foundSwitch = True
                self.flipEdge(e)
                break
            if not foundSwitch:
                raise RuntimeError( "Algorithm failed unexpectedly." )

        #TODO Why is it safe to resolve everything in sequence? Isn't it
        #   possible to make a curve unresolvable after clearing an edge?
        #TODO Replace this with self._resolveAllImpl()?
        # Once we break out of the above loop, we know that every component
        # is either resolvable or isotopic to an off-diagonal. Handle the
        # resolvable components first.
        res = True
        while res:
            res = False
            for i in range( self.countUnresolved() ):
                if self.resolveComponent(i):
                    res = True
                    break

        # Now handle the off-diagonals.
        diag = True
        while diag:
            diag = False
            for i in range( self.countUnresolved() ):
                diagEdgeInds = self._recogniseOffDiagonal(i)
                if diagEdgeInds:
                    diag = True
                    self.flipEdge(
                            self._tri.edge( diagEdgeInds[0] ) )
                    break
        while self.countUnresolved() > 0:
            #TODO Do we need this test?
            if not self.resolveComponent(0):
                raise RuntimeError( "Unexpected unresolvable component." )

    def resolveUntilChoice(self):
        #TODO Document return value. Return None if we finish resolving.
        """
        Attempts to resolve all components of the underlying normal curve,
        but stops when we reach a situation where we have multiple choices.
        """
        while not self._allComponentsNice():
            # Do we have one or more maximal-weight reducible edges?
            maxRed = []
            for e in self._maxReducibleEdges():
                maxRed.append(e)
            if maxRed:
                # We found at least one maximal-weight reducible edge.
                if len(maxRed) == 1:
                    self.flipEdge( maxRed[0] )
                    continue
                else:
                    return ( "Reduce", maxRed )

            # If there is no maximal-weight reducible edge, then check
            # whether we have one or more options for resolving a component.
            res = []
            for r in self._resolutions():
                res.append(r)
            if res:
                # We found at least one resolution.
                if len(res) == 1:
                    self.resolveComponent( *res[0] )
                    continue
                else:
                    return ( "Resolve", res )

            # If we can't resolve either, then we must be at least one
            # maximal-weight edge e such that flipping e reduces the number
            # of switches.
            switch = []
            for e in self._maxReduceSwitchEdges():
                switch.append(e)
            if switch:
                # We found at least one edge that we can flip to reduce the
                # number of switches.
                if len(switch) == 1:
                    self.flipEdge(e)
                    continue
                else:
                    return ( "Switch",  switch )
            else:
                raise RuntimeError( "Algorithm unexpectedly failed." )

        #TODO Why is it safe to resolve everything in sequence? Isn't it
        #   possible to make a curve unresolvable after clearing an edge?
        #TODO Use self._resolveAllImpl()?
        # If we break out of the above loop, then we know that every
        # component is either resolvable or isotopic to an off-diagonal.
        #TODO
        pass

    def _fold( self, e ):
        #TODO This has a bunch of pre-conditions.
        """
        If possible, folds about the given boundary edge e; returns True if
        and only if folding is successful.
        """
        #TODO Following can create invalid vertices. This should just be
        #   handled in the pre-conditions, since the only times I ever call
        #   this routine is when I know it is safe to do so.
        emb = [ e.embedding(0), e.embedding( e.degree() - 1 ) ]
        tet = [ x.tetrahedron() for x in emb ]
        ver = [ x.vertices() for x in emb ]
        if tet[0] == tet[1] and ver[0][3] == ver[1][2]:
            return False
        else:
            tet[0].join( ver[0][3], tet[1],
                    ver[1] * Perm4(2,3) * ver[0].inverse() )
            return True

    #TODO I might want to just delete this. It should be easy to reinstate
    #   later if I change my mind.
    def foldAlongCurve( self, e ):
        """
        Folds along the resolved curve corresponding to the given edge e.

        Topologically, this has the effect of attaching a thickened disc
        along this resolved curve.

        Raises ValueError if e does not correspond to a resolved curve.

        Pre-condition:
        --> e is a boundary edge of the underlying triangulation.
        """
        # Make sure that e corresponds to a resolved curve.
        if e.index() not in self._resolvedEdges:
            raise ValueError(
                    "Given edge does not correspond to a resolved curve." )

        # Perform the fold.
        self._fold( self._flipEdgeImpl(e) )

    def _attemptFold(self):
        """
        Assuming boundary component 0 of the underlying triangulation is a
        2-sphere boundary, tries to either partially or completely fill this
        2-sphere with a 3-ball by folding a pair of boundary faces together.

        Returns True if and only if such a fold is successfully performed.
        """
        # The idea is to iterate through each boundary edge e, and consider
        # the square formed by the two boundary triangles incident to e.
        # We could potentially fold this square across either of its two
        # diagonals (for one of these folds, we first have to do an edge
        # flip), so we check whether these folds are "safe".
        bc = self._tri.boundaryComponent(0)
        built = bc.build()
        flip = None # Edge that we can fold after flipping, if necessary.
        for e in built.edges():
            fac = [ emb.triangle() for emb in e.embeddings() ]
            if fac[0] == fac[1]:
                continue
            ver = [ emb.vertices() for emb in e.embeddings() ]

            # First check whether it is safe to fold directly over e.
            if fac[0].vertex( ver[0][2] ) == fac[1].vertex( ver[1][2] ):
                # The fold is safe, so perform it.
                #NOTE Regina promises an edge-index correspondence between bc
                #   and built, but this sometimes fails for some reason.
                #self._fold( bc.edge( e.index() ) )
                self._fold(
                        bc.triangle( fac[0].index() ).edge( ver[0][2] ) )
                return True

            # Now, if necessary, check whether it would be safe to fold after
            # flipping e.
            if ( flip is None ) and ( e.vertex(0) == e.vertex(1) ):
                # Only perform this fold later on, if it is really necessary.
                #NOTE Regina promises an edge-index correspondence between bc
                #   and built, but this sometimes fails for some reason.
                #flip = bc.edge( e.index() )
                flip = bc.triangle( fac[0].index() ).edge( ver[0][2] )

        # We couldn't directly fold, but maybe we can fold after flipping.
        if flip is None:
            return False
        else:
            # It is safe to fold after flipping.
            newEdge = self.flipEdge(flip)
            self._fold(newEdge)
            return True

    def _fillBall(self):
        """
        Assuming the underlying triangulation has 2-sphere boundary, fills
        this 2-sphere with a 3-ball.
        """
        while self._tri.hasBoundaryTriangles():
            if self._attemptFold():
                continue

            # Since folding was not possible, we know (in particular) that we
            # have no vertices of degree less than three. Our goal now is to
            # keep flipping edges until folding is possible; at worst, this
            # is guaranteed to succeed once we have reduced the degree of a
            # boundary vertex to two. Start by finding a boundary vertex of
            # minimum degree.
            bc = self._tri.boundaryComponent(0)
            built = bc.build()
            minVert = built.vertex(0)
            for i in range( 1, built.countVertices() ):
                v = built.vertex(i)
                deg = v.degree()
                if deg < minVert.degree():
                    minVert = v
                    if deg == 3:
                        break

            # Pre-compute a list of tet-edgeNum pairs that we flip.
            flip = []
            refTetInds = dict()
            for i in range( minVert.degree() - 2 ):
                vertEmb = minVert.embedding(i)
                #NOTE Regina promises an edge-index correspondence between bc
                #   and built, but this sometimes fails for some reason.
                #ind = vertEmb.triangle().edge(
                #        vertEmb.vertices()[2] ).index()
                #edgeEmb = bc.edge(ind).embedding(0)
                facInd = vertEmb.triangle().index()
                edgeNum = vertEmb.vertices()[2]
                edgeEmb = bc.triangle(facInd).edge(edgeNum).embedding(0)

                # Update both flip and refTetInds.
                refTet = edgeEmb.tetrahedron()
                flip.append(
                        ( refTet, edgeEmb.edge() ) )
                if refTet.index() in refTetInds:
                    refTetInds[ refTet.index() ].add( len(flip) - 1 )
                else:
                    refTetInds[ refTet.index() ] = { len(flip) - 1 }

            # Keep flipping until we can fold.
            while True:
                tet, edgeNum = flip.pop()
                self._flipEdgeImpl(
                        tet.edge(edgeNum), [ flip, refTetInds ] )
                if self._attemptFold():
                    break

    def attachOtherHandlebody(self):
        """
        If our Heegard curves are all resolved, then performs the topological
        operation of attaching another handlebody according to these curves.

        Returns True if and only if all our Heegaard curves resolved; this
        routine modifies the underlying triangulation when and only when it
        returns True.
        """
        if self.countUnresolved() > 0:
            return False

        # Initially, we just know that we need to attach a disc to every
        # resolved component by flipping and then folding. However, if this
        # requires us to flip two edges of a single triangle f, then we first
        # need to split f into two triangles by flipping the third edge of f.
        # We use a stack to keep track of all the edges (stored as pairs
        # consisting of a tetrahedron and an edge number) that we need to
        # flip, as well as the order in which we need to flip these edges.
        flipSet = set()
        flipStack = []
        refTetInds = dict()
        for e in self._resolvedEdges:
            edge = self._tri.edge(e)
            flipSet.add( edge.index() )
            emb = edge.embedding(0)
            refTet = emb.tetrahedron()
            flipStack.append( ( refTet, emb.edge() ) )
            if refTet.index() in refTetInds:
                refTetInds[ refTet.index() ].add( len(flipStack) - 1 )
            else:
                refTetInds[ refTet.index() ] = { len(flipStack) - 1 }
        splitFace = True
        splitIndices = set()
        while splitFace:
            splitFace = False

            # For each boundary face f that we have not yet split, check
            # whether it is not necessary to split f.
            for f in self._tri.triangles():
                if ( not f.isBoundary() ) or ( f.index() in splitIndices ):
                    continue

                # If we need to flip two of the edges of f, then we first
                # need to split f by flipping its third edge.
                unflippedEdgeNums = {0,1,2}
                for e in range(3):
                    if f.edge(e).index() in flipSet:
                        unflippedEdgeNums.remove(e)
                if len(unflippedEdgeNums) == 1:
                    edge = f.edge( unflippedEdgeNums.pop() )
                    flipSet.add( edge.index() )
                    emb = edge.embedding(0)

                    # Update both flipStack and refTetInds.
                    refTet = emb.tetrahedron()
                    flipStack.append( ( refTet, emb.edge() ) )
                    if refTet.index() in refTetInds:
                        refTetInds[ refTet.index() ].add(
                                len(flipStack) - 1 )
                    else:
                        refTetInds[ refTet.index() ] = {
                                len(flipStack) - 1 }
                    splitFace = True
                    splitIndices.add( f.index() )

        # Flip every edge in flipStack. Since we know that the last few flips
        # correspond to the resolved curves, we can also fold these now.
        #NOTE The number of components of our underlying curve decreases as
        #   we attach discs, so we need to remember the *initial* number.
        curveCount = self.countComponents()
        while flipStack:
            tet, edgeNum = flipStack.pop()
            newEdge = self._flipEdgeImpl(
                    tet.edge(edgeNum), [ flipStack, refTetInds ] )
            if len(flipStack) < curveCount:
                self._fold(newEdge)

        # To complete the construction, fold until we have no boundary left.
        self._fillBall()
        return True


# Test code.
if __name__ == "__main__":
    #TODO Tidy up tests.
    initTri = Triangulation3.fromIsoSig( "eHbecadjk" )
#    compsMsg = "{} component(s):"
#    edgeMsg = "After layering on edge {}:"
#    subMsg = "    Switches: {}. Resolvable: {}. Off-diagonal: {}."
    #TODO Test invalid inputs too.
    testCases = [
#            (0,0,0,2,2,2,3,1,1),    # 1-component
            (1,0,1,2,3,2,3,2,1),    # 2-component
#            (3,5,4,2,3,5,4,5,3),    # 3-component
#            (1,1,0,1,0,0,0,0,0),    # 1-component
#            (5,5,4,0,5,5,5,4,0),    # 3-component
            (1,2,3,4,5,2,3,4,1) ]   # 2-component
    print()
    for w in testCases:
        oldTri = Triangulation3(initTri)
        hb = HeegaardBuilder( oldTri, w )
        hb.resolveAll()
        hb.attachOtherHandlebody()
        newTri = hb.triangulation()

        #TODO
        # The triangulation newTri should be closed.
        print( "Valid: {}. Closed: {}. Orbl: {}.".format(
            newTri.isValid(), newTri.isClosed(), newTri.isOrientable() ) )
        sim = Triangulation3(newTri)
        sim.intelligentSimplify()
        sim.intelligentSimplify()
        print( "Final size: {}. Original: {}. Simplified: {}.".format(
            newTri.size(), newTri.isoSig(), sim.isoSig() ) )
        print()
#    # Basic operation tests.
#    for w, e in testCases:
#        tri = Triangulation3(initTri)
#        hb = HeegaardBuilder( tri, w )
#
#        # Test components.
#        comps = hb.countUnresolved()
#        print( compsMsg.format(comps) )
#        for c in range(comps):
#            print( subMsg.format( hb._countSwitchesInComponent(c),
#                hb.recogniseResolvable(c), hb._recogniseOffDiagonal(c) ) )
#
#        # Test layering.
#        hb.flipEdge( tri.edge(e) )
#        print( edgeMsg.format(e) )
#        for c in range(comps):
#            print( subMsg.format( hb._countSwitchesInComponent(c),
#                hb.recogniseResolvable(c), hb._recogniseOffDiagonal(c) ) )
#
#        # Test bubble/isotopy.
#        tri = Triangulation3(initTri)
#        hb = HeegaardBuilder( tri, w )
#        iso = hb._isotopeOffEdge( tri.edge(0), 0 )
#        print( "Isotopy: {}.".format(iso) )
#        if iso:
#            print( hb._weights )
#            comps = hb.countUnresolved()
#            for c in range(comps):
#                print( subMsg.format( hb._countSwitchesInComponent(c),
#                    hb.recogniseResolvable(c),
#                    hb._recogniseOffDiagonal(c) ) )
#        #print( "Bubble: {}.".format(
#        #    hb._bubblePoints( tri.edge(0), 0 ) ) )
#        print()
#
#    # Test clearing edge 7 of testCases[4].
#    print( "Clear edge." )
#    tri = Triangulation3(initTri)
#    hb = HeegaardBuilder( tri, testCases[4][0] )
#    hb._clearEdge( tri.edge(7) )
#    print( hb._weights )
#    print()
#
#    # Test resolving (to edge 7 of testCases[4]).
#    print( "Resolve one component." )
#    tri = Triangulation3(initTri)
#    hb = HeegaardBuilder( tri, testCases[4][0] )
#    for i in range(3):
#        res = hb.resolveComponent(0)
#        print( "Resolve component {}: {}.".format( i, res ) )
#        if res:
#            print( "Resolved edge indices: {}.".format(
#                hb._resolvedEdges ) )
#            print( "Weights: {}.".format( hb._weights ) )
#            break
#    print()
#
#    # Test full construction (on testCases[1] and testCases[5]).
#    for c in {1,5}:
#        print( "Resolve all components, case {}.".format(c) )
#        tri = Triangulation3(initTri)
#        hb = HeegaardBuilder( tri, testCases[c][0] )
#        hb.resolveAll()
#        #TODO Rename HeegaardBuilder.triangulation()?
#        print( "Final size: {}. Resolved edges: {}.".format(
#            hb.triangulation().size(), hb._resolvedEdges ) )
#        print()
#        print( "Full construction, case {}.".format(c) )
#        mfd = hb.attachOtherHandlebody()
#        print( "Valid: {}. Closed: {}. Orbl: {}.".format(
#            mfd.isValid(), mfd.isClosed(), mfd.isOrientable() ) )
#        sim = Triangulation3(mfd)
#        sim.intelligentSimplify()
#        sim.intelligentSimplify()
#        print( "Final size: {}. Original: {}. Simplified: {}.".format(
#            mfd.size(), mfd.isoSig(), sim.isoSig() ) )
#        print()
