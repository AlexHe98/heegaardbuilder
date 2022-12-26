"""
A class for building triangulations from a Heegaard diagram.
"""
from regina import *


#TODO Should probably make almost all of the methods underscored, since most
#   of them should never need to be accessed by users.
class HeegaardBuilder:
    #TODO Rewrite documentation.
    #TODO So far, the code doesn't really check for valid inputs.
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
        # - Components that are isotopic to off-diagonals.
        # - Lengths of components.
        # - Number of switches.
        self._resolvableEdges = [None] * self.countUnresolved()
        self._offDiagEdges = [None] * self.countUnresolved()
        self._lengths = [None] * self.countUnresolved()
        self._switchCounts = [None] * self.countUnresolved()

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
        Returns the tuple of edge weights that define the unresolved
        components of the underlying curve.
        """
        return self._weights

    def resolvedEdgeIndices(self):
        """
        Returns a copy of the set of resolved edge indices.
        """
        return set( self._resolvedEdges )

    def countUnresolved(self):
        """
        Returns the number of unresolved components of the underlying curve.
        """
        return len( self._roots )

    def countResolved(self):
        """
        Returns the number of resolved components of the underlying curve.
        """
        return len( self._resolvedEdges )

    def countComponents(self):
        """
        Returns the number of components of the underlying curve.
        """
        return self.countUnresolved() + self.countResolved()

    def _traverseCurve( self, startPt, d ):
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

    def _traverseComponentImpl( self, index ):
        for t in self._traverseCurve( self._roots[index], 0 ):
            yield t

    def traverseComponent( self, index ):
        """
        Iterates through the intersection points of the requested unresolved
        component of the underlying curve.
        """
        for _, p, _, _ in self._traverseComponentImpl(index):
            yield p

    def countSwitchesInComponent( self, index ):
        """
        Returns the number of switches in the requested unresolved component
        of the underlying curve.

        The first call to this routine caches the result, so subsequent calls
        with the same index will not need to recompute the answer.
        """
        if self._switchCounts[index] is None:
            self._switchCounts[index] = 0
            for p in self.traverseComponent(index):
                if self._switch[p]:
                    self._switchCounts[index] += 1
        return self._switchCounts[index]

    def _flipWeightChange( self, e ):
        # Compute how the weight changes after flipping the edge e.
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
        return self._weights[ e.index() ] + self._flipWeightChange(e)

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
        #TODO Optimise!
        self.layerOn(e)

    def componentLength( self, index ):
        """
        Returns the length of the requested unresolved component of the
        underlying curve.

        The first call to this routine caches the result, so subsequent calls
        with the same index will not need to recompute the answer.
        """
        if self._lengths[index] is None:
            self._lengths[index] = 0
            for _ in self._traverseComponentImpl(index):
                self._lengths[index] += 1
        return self._lengths[index]

    #TODO What's the best way to provide user access to arc coordinates?

    def recogniseResolvable( self, index ):
        """
        Recognise the edge indices to which we can resolve the requested
        unresolved component of the underlying curve.

        The first call to this routine caches the result, so subsequent calls
        with the same index will not need to recompute the answer.
        """
        if self._resolvableEdges[index] is None:
            switches = 0
            resEdgeInds = []
            for info in self._traverseComponentImpl(index):
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

    def recogniseOffDiagonal( self, index ):
        """
        Recognise the edge indices across which the requested unresolved
        component forms an off-diagonal.

        The first call to this routine caches the result, so subsequent calls
        with the same index will not need to recompute the answer.
        """
        if self._offDiagEdges[index] is None:
            switches = 0
            crossDiagInds = []
            for info in self._traverseComponentImpl(index):
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
            #TODO This check shouldn't be necessary, but leave this here
            #   until I have implemented checks in more appropriate places.
            if not self._isotopeOffEdge( e, 0 ):
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

    def resolveComponent( self, index ):
        """
        Checks whether it is possible to resolve the requested unresolved
        component, and if so resolves this component.
        """
        resEdgeInds = self.recogniseResolvable(index)
        if not resEdgeInds:
            return False

        # First compute how edge weights would reduce as a result of
        # resolving the requested component.
        weightChanges = [0] * self._tri.countEdges()
        for _, p, _, _ in self._traverseComponentImpl(index):
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

    def isHeegaard(self):
        """
        Does the underlying curve give a valid Heegaard diagram?
        """
        #TODO Check handlebody, number of curves, and separating?
        pass

    def allComponentsNice(self):
        """
        Are all components of the underlying curve "nice" in the sense that
        they are either resolved, resolvable or isotopic to an off-diagonal?
        """
        for i in range( self.countUnresolved() ):
            if ( not self.recogniseResolvable(i) and
                    not self.recogniseOffDiagonal(i) ):
                return False
        return True

    def isMaximalWeight( self, e ):
        """
        Is the edge e maximal-weight?

        If e is not even a boundary edge, then this routine always returns
        False. Otherwise, this routine returns True if and only if the weight
        of e is maximal among boundary edges.
        """
        return ( self._weights[ e.index() ] == max( self._weights ) )

    def isReducible( self, e ):
        """
        Is the edge e reducible?
        """
        return ( self._flipWeightChange(e) < 0 )

    def resolveAll(self):
        """
        Resolves all components of the underlying curve.
        """
        while not self.allComponentsNice():
            # Try to flip a maximal-weight reducible edge.
            foundMinRed = False
            for e in self._tri.edges():
                if ( not e.isBoundary() or not self.isMaximalWeight(e) or
                        not self.isReducible(e) ):
                    continue
                foundMinRed = True
                self.flipEdge(e)
                break
            if foundMinRed:
                continue

            # If there is no maximal-weight reducible edge, then try to
            # resolve a component.
            resolved = False
            for i in range( self.countUnresolved() ):
                if self.resolveComponent(i):
                    resolved = True
                    break
            if resolved:
                continue

            # If we can't resolve either, then we must be able to reduce the
            # number of switches by flipping a maximal-weight edge.
            foundSwitch = False
            for e in self._tri.edges():
                if not e.isBoundary() or not self.isMaximalWeight(e):
                    continue

                # Check whether flipping this maximal-weight edge e would
                # reduce the total number of switches.
                for i in range( self._weights[ e.index() ] ):
                    pt = ( e.index(), i )
                    if ( self._switch[pt] and
                            self._switch[ self._adjPoints[pt][0] ] and
                            self._switch[ self._adjPoints[pt][1] ] ):
                        foundSwitch = True
                        self.flipEdge(e)
                        break
                if foundSwitch:
                    break

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
        #TODO Optimise?
        diag = True
        while diag:
            diag = False
            for i in range( self.countUnresolved() ):
                diagEdgeInds = self.recogniseOffDiagonal(i)
                if diagEdgeInds:
                    diag = True
                    self.flipEdge( self._tri.edge( diagEdgeInds[0] ) )
                    break
        while self.countUnresolved() > 0:
            #TODO Do we need this test?
            if not self.resolveComponent(0):
                raise RuntimeError( "Unexpected unresolvable component." )

    #TODO This is effectively a static method, so possibly this doesn't
    #   belong here.
    def fold( self, e ):
        #TODO This has a bunch of pre-conditions.
        """
        If possible, folds about the given boundary edge e; returns True if
        and only if folding is successful.
        """
        #TODO Following can create invalid vertices.
        emb = [ e.embedding(0), e.embedding( e.degree() - 1 ) ]
        tet = [ x.tetrahedron() for x in emb ]
        ver = [ x.vertices() for x in emb ]
        if tet[0] == tet[1] and ver[0][3] == ver[1][2]:
            return False
        else:
            tet[0].join( ver[0][3], tet[1],
                    ver[1] * Perm4(2,3) * ver[0].inverse() )
            return True

    #TODO This is awkward with the current implementation. I should probably
    #   be careful to distinguish different "states" for self._tri depending
    #   on whether we have started "closing up".
    def _attemptFold(self):
        # Assuming boundary component 0 of the underlying triangulation is a
        # 2-sphere boundary, tries to (at least partially) fill this 2-sphere
        # with a 3-ball by folding a pair of boundary faces together.
        #
        # Returns True if and only if such a fold is successfully performed.
        #
        # The idea is to iterate through each boundary edge e, and consider
        # the square formed by the two boundary triangles incident to e.
        # We could potentially fold this square across either of its two
        # diagonals (for one of these folds, we have to layer one tetrahedron
        # first), so we check whether these folds are "safe".
        bc = self._tri.boundaryComponent(0)
        built = bc.build()
        layer = None # Edge that we can fold after layering, if necessary.
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
                #self.fold( bc.edge( e.index() ) )
                self.fold(
                        bc.triangle( fac[0].index() ).edge( ver[0][2] ) )
                return True

            # Now, if necessary, check whether it would be safe to fold after
            # layering across e.
            if ( layer is None ) and ( e.vertex(0) == e.vertex(1) ):
                # Only perform this fold later on, if it is really necessary.
                #NOTE Regina promises an edge-index correspondence between bc
                #   and built, but this sometimes fails for some reason.
                #layer = bc.edge( e.index() )
                layer = bc.triangle( fac[0].index() ).edge( ver[0][2] )

        # We couldn't directly fold, but maybe we can fold after layering.
        if layer is None:
            return False
        else:
            newTet = self._tri.layerOn(layer)
            self.fold( newTet.edge(5) )
            return True

    def fillBall(self):
        """
        Assuming the underlying triangulation has 2-sphere boundary, fills
        this 2-sphere with a 3-ball.
        """
        while self._tri.hasBoundaryTriangles():
            if self._attemptFold():
                continue

            # Since folding was not possible, we know (in particular) that we
            # have no vertices of degree less than three. Our goal now is to
            # repeatedly layer until folding is possible; at worst, this is
            # guaranteed to succeed once we have reduced the degree of a
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

            # Pre-compute a list of tet-edgeNum pairs across which we layer.
            layer = []
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
                layer.append(
                        ( edgeEmb.tetrahedron(), edgeEmb.edge() ) )

            # Keep layering until we can fold.
            while True:
                tet, edgeNum = layer.pop()
                self._tri.layerOn( tet.edge(edgeNum) )
                if self._attemptFold():
                    break

    #TODO This is awkward with the current implementation. I should probably
    #   be careful to distinguish different "states" for self._tri depending
    #   one whether we have started "closing up".
    def constructManifold(self):
        """
        Construct a triangulation of the closed 3-manifold represented by
        the underlying Heegaard diagram.
        """
        # Initially, we just know that we need to flip every edge that
        # corresponds to a resolved component. However, if we need to flip
        # two edges of a single triangle f, then we first need to split f
        # into two triangles by flipping the third edge of f. We use a stack
        # to keep track of all the edges (stored as pairs consisting of a
        # tetrahedron and an edge number) that we need to flip, as well as
        # the order in which we need to flip these edges.
        flipSet = set()
        flipStack = []
        for e in self.resolvedEdgeIndices():
            edge = self._tri.edge(e)
            flipSet.add( edge.index() )
            emb = edge.embedding(0)
            flipStack.append( ( emb.tetrahedron(), emb.edge() ) )
        splitFace = True
        splitIndices = set()
        while splitFace:
            splitFace = False

            # For each boundary face f that we have not yet split, check
            # whether it is now necessary to split f.
            for f in self._tri.triangles():
                if ( not f.isBoundary() ) or ( f.index() in splitIndices ):
                    continue

                # If we need to flip two of the edges of f, then we need to
                # split f by first flipping its third edge.
                unflippedEdgeNums = {0,1,2}
                for e in range(3):
                    edge = f.edge(e)
                    if edge.index() in flipSet:
                        unflippedEdgeNums.remove(e)
                if len(unflippedEdgeNums) == 1:
                    edge = f.edge( unflippedEdgeNums.pop() )
                    flipSet.add( edge.index() )
                    emb = edge.embedding(0)
                    flipStack.append( ( emb.tetrahedron(), emb.edge() ) )
                    splitFace = True
                    splitIndices.add( f.index() )

        # Flip every edge in flipStack. Since we know that the last few flips
        # correspond to the resolved curves, we can also fold these now.
        while flipStack:
            tet, edgeNum = flipStack.pop()
            #TODO Optimise!
            newTet = self._tri.layerOn( tet.edge(edgeNum) )
            if len(flipStack) < self.countComponents():
                self.fold( newTet.edge(5) )

        # To complete the construction, fold until we have no boundary left.
        self.fillBall()
        #TODO Return the triangulation, or just a copy?
        return self._tri


# Test code.
if __name__ == "__main__":
    #TODO Tidy up tests.
    initTri = Triangulation3.fromIsoSig( "eHbecadjk" )
    compsMsg = "{} component(s):"
    edgeMsg = "After layering on edge {}:"
    subMsg = "    Switches: {}. Resolvable: {}. Off-diagonal: {}."
    testCases = [
            ( (0,0,0,2,2,2,3,1,1), 0 ),     # 1-component
            ( (1,0,1,2,3,2,3,2,1), 4 ),     # 2-component
            ( (3,5,4,2,3,5,4,5,3), 2 ),     # 3-component
            ( (1,1,0,1,0,0,0,0,0), 3 ),     # 1-component
            ( (5,5,4,0,5,5,5,4,0), 0 ),     # 3-component
            ( (1,2,3,4,5,2,3,4,1), 0 ) ]    # 2-component

    # Basic operation tests.
    for w, e in testCases:
        tri = Triangulation3(initTri)
        hb = HeegaardBuilder( tri, w )

        # Test components.
        comps = hb.countUnresolved()
        print( compsMsg.format(comps) )
        for c in range(comps):
            print( subMsg.format( hb.countSwitchesInComponent(c),
                hb.recogniseResolvable(c), hb.recogniseOffDiagonal(c) ) )

        # Test layering.
        hb.layerOn( tri.edge(e) )
        print( edgeMsg.format(e) )
        for c in range(comps):
            print( subMsg.format( hb.countSwitchesInComponent(c),
                hb.recogniseResolvable(c), hb.recogniseOffDiagonal(c) ) )

        # Test bubble/isotopy.
        tri = Triangulation3(initTri)
        hb = HeegaardBuilder( tri, w )
        iso = hb._isotopeOffEdge( tri.edge(0), 0 )
        print( "Isotopy: {}.".format(iso) )
        if iso:
            print( hb.weights() )
            comps = hb.countUnresolved()
            for c in range(comps):
                print( subMsg.format( hb.countSwitchesInComponent(c),
                    hb.recogniseResolvable(c),
                    hb.recogniseOffDiagonal(c) ) )
        #print( "Bubble: {}.".format(
        #    hb._bubblePoints( tri.edge(0), 0 ) ) )
        print()

    # Test clearing edge 7 of testCases[4].
    print( "Clear edge." )
    tri = Triangulation3(initTri)
    hb = HeegaardBuilder( tri, testCases[4][0] )
    hb._clearEdge( tri.edge(7) )
    print( hb.weights() )
    print()

    # Test resolving (to edge 7 of testCases[4]).
    print( "Resolve one component." )
    tri = Triangulation3(initTri)
    hb = HeegaardBuilder( tri, testCases[4][0] )
    for i in range(3):
        res = hb.resolveComponent(0)
        print( "Resolve component {}: {}.".format( i, res ) )
        if res:
            print( "Resolved edge indices: {}.".format(
                hb.resolvedEdgeIndices() ) )
            print( "Weights: {}.".format( hb.weights() ) )
            break
    print()

    # Test full construction (on testCases[1] and testCases[5]).
    for c in {1,5}:
        print( "Resolve all components, case {}.".format(c) )
        tri = Triangulation3(initTri)
        hb = HeegaardBuilder( tri, testCases[c][0] )
        hb.resolveAll()
        print( "Final size: {}. Resolved edges: {}.".format(
            hb.triangulation().size(), hb.resolvedEdgeIndices() ) )
        print()
        print( "Full construction, case {}.".format(c) )
        mfd = hb.constructManifold()
        print( "Valid: {}. Closed: {}. Orbl: {}.".format(
            mfd.isValid(), mfd.isClosed(), mfd.isOrientable() ) )
        sim = Triangulation3(mfd)
        sim.intelligentSimplify()
        sim.intelligentSimplify()
        print( "Final size: {}. Original: {}. Simplified: {}.".format(
            mfd.size(), mfd.isoSig(), sim.isoSig() ) )
        print()
