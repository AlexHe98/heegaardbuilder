"""
A class for building triangulations from a Heegaard bouquet.
"""
from regina import *


#############################################################################
#                              Helper routines                              #
#                              ===============                              #
#############################################################################


def _faceNumberings(e):
    """
    Returns the triangle indices and vertex numberings for the boundary
    triangles incident to the given boundary edge e.

    In detail, this routine returns a pair (f, v) such that:
    --> f[0] is the index of the boundary triangle given by
            e.embedding(0).tetrahedron().triangle(3);
    --> f[1] is the index of the boundary triangle given by
            e.embedding( e.degree() - 1 ).tetrahedron().triangle(2);
    --> for i,j in {0,1}, v[i][j] is the vertex number of triangle f[i]
        corresponding to vertex j of the edge e; and
    --> for i in {0,1}, v[i][2] is the vertex number of triangle f[i]
        that is opposite the edge e.

    Pre-condition;
    --> e is a boundary edge of a 3-manifold triangulation.
    """
    emb = [ e.embedding(0), e.embedding( e.degree() - 1 ) ]
    f = []
    v = []
    p = [ Perm4(), Perm4(2,3) ]
    for i in range(2):
        face = emb[i].tetrahedron().triangle( emb[i].vertices()[3-i] )
        f.append( face.index() )
        perm = ( face.embedding(0).vertices().inverse() *
                emb[i].vertices() * p[i] )
        v.append( [ perm[0], perm[1], perm[2] ] )
    return (f, v)


def _checkTri(tri):
    """
    Checks that the triangulation tri is valid, orientable, one-vertex
    and has exactly one boundary component.
    """
    if not tri.isValid():
        raise ValueError( "Triangulation is not valid." )
    if not tri.isOrientable():
        raise ValueError( "Triangulation is not orientable." )
    if tri.countVertices() != 1:
        raise ValueError( "Triangulation is not one-vertex." )
    bdryComps = tri.countBoundaryComponents()
    if bdryComps != 1:
        raise ValueError( "Triangulation must have exactly one boundary " +
                "component, not {}.".format(bdryComps) )


def _checkBouquet( tri, weights, resolved ):
    """
    Checks that weights and resolved are formatted correctly.

    Specifically, weights and resolved must satisfy the following conditions:
    --> weights must be a list or tuple consisting of n non-negative
        integers, where n = tri.countEdges();
    --> for each i, weights[i] must be 0 if tri.edge(i) is internal;
    --> resolved must be a (possibly empty) set consisting of indices of
        boundary edges of tri; and
    --> for each i in resolved, weights[i] must be 0.
    These restrictions on the format are the same as conditions (2), (3),
    (5) and (6) in the documentation for the HeegaardBuilder.setBouquet()
    routine.

    Pre-condition:
    --> tri is valid, orientable, one-vertex, and has exactly one boundary
        component.
    """
    # Check that weights and resolved have the correct format.
    givenNum = len(weights)
    requiredNum = tri.countEdges()
    if givenNum != requiredNum:
        raise ValueError( "Given {} edge weights, but ".format(givenNum) +
                "the triangulation has {} edges.".format(requiredNum) )
    for i, wt in enumerate(weights):
        if wt < 0:
            raise ValueError( "Edge weights must be non-negative, but " +
                    "edge {} has weight {}.".format( i, wt ) )
        elif wt > 0:
            if not tri.edge(i).isBoundary():
                raise ValueError( "Edge {} is internal, so ".format(i) +
                        "its weight must be 0, not {}.".format(wt) )
    for i in resolved:
        if ( i < 0 or i >= tri.countEdges() or
                not tri.edge(i).isBoundary() ):
            raise ValueError( "Edge {} is internal, so it ".format(i) +
                    "cannot be assigned as a resolved edge." )
        if weights[i] > 0:
            raise ValueError( "Edge {} is resolved, so its ".format(i) +
                    "weight must be 0, not {}.".format( weights[i] ) )


def _computeArcCoords( tri, weights ):
    """
    Computes the normal and root arc coordinates defined by the given
    edge weights, and returns these coordinates together with some
    auxiliary information.

    The arc coordinates only make sense if the given edge weights satisfy the
    matching constraints for a Heegaard bouquets. Specifically, the matching
    constraints require that for each triangular face F in the boundary of
    tri, either:
    (a) one edge of F contributes more than half of the total weight of all
        three edges of F; or
    (b) the total weight of the three edges of F is even.
    If condition (a) holds, then F contains one or more root arcs. Otherwise,
    if condition (a) fails, then F must contain only normal arcs, in which
    case condition (b) must hold. This routine raises ValueError if the
    matching constraints fail.

    The auxiliary information consists of a tuple with the following entries:
    (0) The total number of root arcs.
    (1) A Boolean that is True if and only if at least one arc coordinate
        is nonzero.

    Pre-condition:
    --> tri is valid, orientable, one-vertex, and has exactly one
        boundary component.
    --> weights is formatted correctly, as specified in the documentation for
        the _checkBouquet() routine.
    """
    arcCoords = []
    totalRootArcs = 0
    hasNonzero = False
    for face in tri.triangles():
        if not face.isBoundary():
            # No arcs in an internal face.
            arcCoords.append(None)
            continue

        wt = [ weights[ face.edge(i).index() ] for i in range(3) ]
        totalWt = sum(wt)
        normal = [0,0,0]    # Normal arc coordinates.
        root = [0,0,0]      # Root arc coordinates.

        # If one edge of this face contributes more than half of the
        # total weight, then this face contains root arcs.
        maxWt = 0
        maxInd = 0
        for i in range(3):
            if wt[i] > maxWt:
                maxWt = wt[i]
                maxInd = i
        rootArcs = 2*maxWt - totalWt
        if rootArcs > 0:
            root[maxInd] = rootArcs
            normal[ maxInd - 1 ] = wt[ maxInd - 2 ]
            normal[ maxInd - 2 ] = wt[ maxInd - 1 ]

            # Update arcCoords, totalRootArcs and hasNonzero.
            arcCoords.append( normal + root )
            totalRootArcs += rootArcs
            hasNonzero = True
            continue

        # Otherwise, this face must contain only normal arcs, in which
        # case the total weight must be even.
        if totalWt % 2 != 0:
            raise ValueError( "Edge weights fail to satisfy the matching " +
                    "constraints in triangle {}.".format( face.index() ) )
        for i in range(3):
            normal[i] = ( wt[i-1] + wt[i-2] - wt[i] ) // 2
            if normal[i] > 0:
                hasNonzero = True
        arcCoords.append( normal + root )

    # We exited the loop, which means that the matching constraints are
    # satisfied and we successfully compute arc coordinates.
    return ( arcCoords, ( totalRootArcs, hasNonzero ) )


def _checkTangential( tri, resolved ):
    """
    Letting R denote the collection of boundary edges of tri given by the
    list of resolved edge indices, this routine checks that all the edges in
    R meet tangentially at the vertex of tri.

    Pre-condition:
    --> tri is valid, orientable, one-vertex, and has exactly one boundary
        component.
    --> resolved is formatted correctly, as specified in the documentation
        for the _checkBouquet() routine.
    --> The Heegaard bouquet given by the edges in R is combinatorially
        admissible.
    """
    stack = []
    bc = tri.boundaryComponent(0)
    built = bc.build()
    for emb in built.vertex(0).embeddings():
        triFace = bc.triangle( emb.triangle().index() )
        edgeInd = triFace.edge( emb.vertices()[2] ).index()
        if edgeInd in resolved:
            if edgeInd in stack:
                topInd = stack.pop()
                if topInd != edgeInd:
                    raise ValueError( "Edges {} and ".format(topInd) +
                            "{} form Heegaard petals ".format(edgeInd) +
                            "that meet transversely." )
            else:
                stack.append(edgeInd)


def _checkComplement( tri, resolved ):
    """
    Letting R denote the collection of boundary edges of tri given by the
    list of resolved edge indices, this routine checks that the boundary
    surface of tri remains connected after cutting along the edges in R.

    Pre-condition:
    --> tri is valid, orientable, one-vertex, and has exactly one boundary
        component.
    --> resolved is formatted correctly, as specified in the documentation
        for the _checkBouquet() routine.
    --> The Heegaard bouquet given by the edges in R is combinatorially
        admissible.
    """
    bc = tri.boundaryComponent(0)

    # Use union-find to count components.
    parent = dict()
    size = dict()
    components = bc.countTriangles()
    for f in bc.triangles():
        ind = f.index()
        parent[ind] = ind
        size[ind] = 1

    # For each boundary edge e that does not form a resolved Heegaard
    # petal, take the union of the two components on either side of e.
    for e in bc.edges():
        if e.index() in resolved:
            continue

        # Find representatives for the components on either side of e.
        f, _ = _faceNumberings(e)
        for i in range(2):
            while parent[ f[i] ] != f[i]:
                # Use path-halving: replace the current parent of f[i]
                # with its current grandparent, and then go to this
                # grandparent.
                parent[ f[i] ] = parent[ parent[ f[i] ] ]
                f[i] = parent[ f[i] ]

        # Since we are *not* cutting along the edge e, take the union of
        # the components on either side.
        if f[0] == f[1]:
            continue
        if size[ f[0] ] < size[ f[1] ]:
            f[0], f[1] = f[1], f[0]
        parent[ f[1] ] = f[0]
        size[ f[0] ] = size[ f[0] ] + size[ f[1] ]
        components -= 1

        # If we are already down to one component, then we are done.
        if components == 1:
            return

    # If we exit the above loop, then we never managed to get the number
    # of components down to one.
    raise ValueError( "After cutting along the Heegaard bouquet, the " +
            "boundary surface should remain connected. However, it " +
            "actually has {} components.".format(components) )


def _checkAdmissible( tri, resolved ):
    """
    Letting R denote the collection of boundary edges of tri given by the
    list of resolved edge indices, this routine checks that R describes a
    (topologically) admissibly Heegaard bouquet in the boundary of tri.

    Pre-condition:
    --> tri is valid, orientable, one-vertex, and has exactly one boundary
        component.
    --> resolved is formatted correctly, as specified in the documentation
        for the _checkBouquet() routine.
    --> The Heegaard bouquet given by the edges in R is combinatorially
        admissible.
    """
    _checkTangential( tri, resolved )
    _checkComplement( tri, resolved )


def _recogniseLayer(e):
    """
    Checks whether e is the new boundary edge that results from layering a
    tetrahedron t across an edge ee, and if so returns details of this
    layering; otherwise, returns None.

    Specifically, if we do have such a layering, then this routine returns a
    a 2-tuple containing the following items:
    (0) the tetrahedron t; and
    (1) an edge-embedding of the edge ee that references a tetrahedron tt
        such that tt != t.
    """
    # The edge e must be a boundary edge of t of degree 1.
    if not e.isBoundary() or e.degree() != 1:
        return None

    # The opposite edge must be internal, and the opposite triangles must be
    # internal and distinct.
    tet = e.embedding(0).tetrahedron()
    ver = e.embedding(0).vertices()
    oppEdge = tet.edge( ver[2], ver[3] )
    if oppEdge.isBoundary():
        return None
    oppTri = [ tet.triangle( ver[i] ) for i in range(2) ]
    if ( oppTri[0].isBoundary() or oppTri[1].isBoundary() or
            oppTri[0] == oppTri[1] ):
        return None

    # Find an edge-embedding of oppEdge that references a tetrahedron other
    # than tet.
    for oppEmb in oppEdge.embeddings():
        if oppEmb.tetrahedron() != tet:
            break
    return ( tet, oppEmb )


def _addEdgeReference( edge, edgeRefs ):
    """
    Update edgeRefs to include the given edge.

    edgeRefs must be a 2-element list consisting of the following items:
    (0) A list F of edge-embeddings of distinct boundary edges of some
        3-manifold triangulation T. (The intention is that F is a stack of
        edges that we wish to flip one by one, which requires that we keep
        track of these edges as we perform flips.)
    (1) A dictionary D that maps tetrahedron indices i to the set of
        locations in the list F that reference tetrahedron i of the
        triangulation T. We allow D to reference locations that no longer
        exist in F because the corresponding element was popped from the end
        of F at some point.

    This routine adds an edge-embedding of the given edge to the list F, and
    updates the dictionary D accordingly.
    """
    flip, refTetInds = edgeRefs
    emb = edge.embedding(0)
    tetInd = emb.tetrahedron().index()
    if tetInd in refTetInds:
        refTetInds[tetInd].add( len(flip) )
    else:
        refTetInds[tetInd] = { len(flip) }
    flip.append(emb)


def _adjustEdgeReferences( edgeRefs, doomed ):
    """
    Adjust edgeRefs so that it no longer references the doomed tetrahedron.

    edgeRefs must be a 2-element list consisting of the following items:
    (0) A list F of edge-embeddings of distinct boundary edges of some
        3-manifold triangulation T. (The intention is that F is a stack of
        edges that we wish to flip one by one, which requires that we keep
        track of these edges as we perform flips.)
    (1) A dictionary D that maps tetrahedron indices i to the set of
        locations in the list F that reference tetrahedron i of the
        triangulation T. We allow D to reference locations that no longer
        exist in F because the corresponding element was popped from the end
        of F at some point.

    This routine adjusts the edge-embeddings in the list F so that none of
    these edge-embeddings references the doomed tetrahedron, and updates the
    dictionary D accordingly.
    """
    flip, refTetInds = edgeRefs

    # Make sure we stop referencing the doomed tetrahedron.
    if doomed.index() in refTetInds:
        for loc in refTetInds[ doomed.index() ]:
            if loc >= len(flip):
                # This can happen when flip is intended to be a stack, in
                # which case refTetInds may reference an element that
                # has been removed from flip.
                continue

            # Find new edge-embedding that avoids the doomed tetrahedron.
            oldEmb = flip[loc]
            oldRefTet = oldEmb.tetrahedron()
            oldRefNum = oldEmb.edge()
            refEdge = oldRefTet.edge(oldRefNum)
            for newEmb in refEdge.embeddings():
                newRefTet = newEmb.tetrahedron()
                if newRefTet != doomed:
                    break

            # Update both flip and refTetInds.
            flip[loc] = newEmb
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


#############################################################################
#                     Implementation of HeegaardBuilder                     #
#                     =================================                     #
#############################################################################


class HeegaardBuilder:
    """
    Implements the various subroutines needed to execute an algorithm for
    taking an orientable one-vertex triangulation T with a single boundary
    component, and filling with a handlebody according to a system of curves
    given by a Heegaard bouquet in the boundary of T.

    A Heegaard bouquet is a bouquet of circles embedded in the boundary of T
    such that the vertex of the bouquet lies on the vertex of T, and each
    petal p (called a Heegaard petal) is a circle embedded in the boundary of
    T in one of two ways:
    (1) The petal p is resolved, meaning that it coincides with an edge of T.
    (2) The petal p is unresolved, meaning that it is embedded so that:
        --> away from the vertex, p meets edges of T transversely; and
        --> p meets the interior of each triangular face F in the boundary of
            T in finitely many normal arcs, together with up to two root arcs
            that go from a vertex of F to the interior of the edge opposite
            this vertex.

    We call a Heegaard bouquet "(topologically) admissible" if all of its
    petals meet tangentially at the vertex (in other words, the petals can be
    made disjoint after isotoping them away from the vertex), and these
    petals form a system of Heegaard curves on the boundary of T. For a
    weaker condition that is much simpler to test, call a Heegaard bouquet
    "combinatorially admissible" if it has exactly 2*(g - r) root arcs, where
    g is the genus of the boundary surface of T and r is the number of
    resolved petals.

    If T is a triangulation of a genus-g handlebody, then T with an
    admissible Heegaard bouquet describes a genus-g Heegaard splitting of an
    orientable closed 3-manifold M, and filling with a handlebody gives a
    triangulation T' of M in which the Heegaard surface appears as a
    subcomplex of the 2-skeleton of T'.

    Moreover, if T is actually a *layered* triangulation of a genus-g
    handlebody (as described in Jaco and Rubinstein's unpublished 2006 paper
    titled "Layered-triangulations of 3-manifolds"), then the method of
    construction implemented by this class guarantees that the triangulation
    T' will have cutwidth at most 4g-2. Since cutwidth is an upper bound for
    treewidth, this means that T' will have treewidth at most 4g-2.
    """
    def __init__(self):
        """
        Initialises an empty HeegaardBuilder.
        """
        self.reset()

    def reset(self):
        """
        Resets self to an empty HeegaardBuilder.
        """
        self._tri = None
        self._weights = None
        self._resolved = None
        self._arcCoords = None
        self._genus = None

    def isEmpty(self):
        """
        Is this HeegaardBuilder currently empty?
        """
        return ( self._tri is None )

    def triangulation(self):
        """
        Returns the triangulation currently stored by this HeegaardBuilder,
        or None if this HeegaardBuilder is empty.

        Warning:
        --> You are allowed to modify the returned triangulation using a
            routine other than one of the routines provided by this
            HeegaardBuilder, but doing so will invalidate the Heegaard
            bouquet currently stored by this HeegaardBuilder; as a result,
            subsequent attempts to work with this Heegaard bouquet will lead
            to undefined behaviour.
        """
        return self._tri

    def _setBouquetImpl( self, tri, weights, resolved, arcCoords,
            arcInfo=None ):
        """
        Set a Heegaard bouquet in the given triangulation tri with the given
        edge weights, resolved edge indices, and arc coordinates.

        The optional arcInfo variable, if provided, should be a 2-tuple
        containing auxiliary information as returned by the
            HeegaardBuilder._computeArcCoords()
        routine. In the case where this auxiliary information is provided,
        this routine will raise ValueError if there are no root arcs, but
        there is still at least one normal arc (since this would imply that
        we have a normal curve).

        This routine will also raise ValueError if every Heegaard petal is
        resolved, and the Heegaard bouquet is not (topologically) admissible.
        """
        # Check for normal curves if it is easy to do so. In particular, this
        # check only occurs if every Heegaard petal is resolved.
        if arcInfo is not None:
            rootArcs, hasNonzero = arcInfo
            if rootArcs == 0 and hasNonzero:
                raise ValueError( "Found a normal curve after resolving " +
                "all Heegaard petals." )

        # If necessary, check that the Heegaard bouquet is admissible.
        if tri.isValid() and tri.countBoundaryComponents() == 1:
            genus = ( 2 - tri.boundaryComponent(0).eulerChar() ) // 2

            # Check admissibility if every Heegaard petal is resolved.
            if len(resolved) == genus:
                _checkAdmissible( tri, resolved )
        else:
            genus = None

        # Only set Heegaard bouquet if all checks pass.
        self._tri = tri
        self._weights = weights
        self._resolved = resolved
        self._arcCoords = arcCoords
        self._genus = genus

    def setBouquet( self, tri, weights, resolved=set() ):
        """
        Set a combinatorially admissible Heegaard bouquet in the given
        triangulation tri by directly specifying edge weights and resolved
        edge indices.

        The input must be "reasonable" in the following sense:
        (1) tri must be valid, orientable, one-vertex, and have exactly one
            boundary component;
        (2) weights must be a list or tuple consisting of n non-negative
            integers, where n = tri.countEdges();
        (3) for each i, weights[i] must be 0 if tri.edge(i) is internal;
        (4) if we assign the value weights[i] to tri.edge(i) for every i such
            that tri.edge(i) is boundary, then these values satisfy the
            matching constraints (described below) for a Heegaard bouquet in
            the boundary of tri;
        (5) resolved must be a set (empty by default) consisting of indices
            of boundary edges tri; and
        (6) for each i in resolved, weights[i] must be 0.
        If the input is not reasonable, then this routine raises ValueError,
        and does not change the currently stored Heegaard bouquet.

        The matching constraints require that for each triangular face F in
        the boundary of tri, either:
        (a) one edge of F contributes more than half of the total weight of
            all three edges of F; or
        (b) the total weight of the three edges of F is even.
        If condition (a) holds, then F contains one or more root arcs.
        Otherwise, if condition (a) fails, then F must contain only normal
        arcs, in which case condition (b) must hold.

        Assuming the input is reasonable, we obtain either:
        --> a Heegaard bouquet; or
        --> a union of a Heegaard bouquet with a (possibly disconnected)
            normal curve.
        In general, this routine does not check for the presence of normal
        curves. The only exception is when the Heegaard bouquet consists
        entirely of resolved petals, in which case normal curves exist if and
        only if we have one or more nonzero edge weights; this routine raises
        ValueError if it detects normal curves in this particular case.

        This routine also raises ValueError if the Heegaard bouquet is not
        combinatorially admissible. However, for (topological) admissibility,
        this routine only raises ValueError if this condition fails in the
        special case where every Heegaard petal is resolved.

        Warning:
        --> This class stores direct references to tri, weights and resolved.
            You are allowed to modify these objects using a routine other
            than one of the routines provided by this HeegaardBuilder, but
            doing so will invalidate the Heegaard bouquet; as a result,
            subsequent attempts to work with this Heegaard bouquet will lead
            to undefined behaviour.
        """
        # If the input is not "reasonable", then the following routines will
        # raise ValueError.
        _checkTri(tri)
        _checkBouquet( tri, weights, resolved )
        arcCoords, arcInfo = _computeArcCoords( tri, weights )

        # Also raise ValueError if the Heegaard bouquet is not
        # combinatorially admissible.
        rootArcs = arcInfo[0]
        chi = tri.boundaryComponent(0).eulerChar() # Equal to (2 - 2*g)
        if rootArcs != 2 - chi - 2*len(resolved):
            # Equivalently, rootArcs is not equal to 2*(g - r), where g is
            # the genus of the boundary surface of tri, and r is the number
            # of resolved petals.
            raise ValueError( "Heegaard bouquet is not combinatorially " +
                    "admissible." )

        # Since we have checked all the necessary conditions, we can now set
        # the Heegaard bouquet. Note that the following routine will, if
        # possible, use the optional arcInfo variable to detect the presence
        # of normal curves.
        self._setBouquetImpl( tri, weights, resolved, arcCoords, arcInfo )

    def setClone( self, other ):
        """
        Set a Heegaard bouquet by cloning the triangulation and bouquet
        currently stored by the other HeegaardBuilder.

        If the other HeegaardBuilder is empty, then this is equivalent to
        resetting this HeegaardBuilder to be empty.
        """
        if other.isEmpty():
            self.reset()
        else:
            # Pass clones of everything to the _setBouquetImpl() routine.
            arcCoords = []
            for a in other._arcCoords:
                if a is None:
                    arcCoords.append(None)
                else:
                    arcCoords.append( list(a) )
            self._setBouquetImpl(
                    Triangulation3( other._tri ), list( other._weights ),
                    set( other._resolved ), arcCoords )

    def countResolved(self):
        """
        Returns the number of resolved petals in the current Heegaard
        bouquet, or None if this HeegaardBuilder is empty.
        """
        if self.isEmpty():
            return None
        else:
            return len( self._resolved )

    def flipWeight( self, e ):
        """
        Computes the weight of the new edge that would result from flipping
        the given edge e.

        Pre-condition;
        --> e is a boundary edge of self.triangulation().
        """
        f, v = _faceNumberings(e)

        # Weight of the new edge receives contributions of three types:
        # - normal arcs opposite e;
        # - root arcs disjoint from the interior of e; and
        # - pairs of normal arcs forming segments that join opposite edges of
        #   the square given by the two triangles incident to e.
        # Contributions of the first two types are relatively straightforward
        # to deal with.
        newWt = 0
        for i in range(2):
            # Contribution from normal arcs opposite e in triangle f[i].
            newWt += self._arcCoords[ f[i] ][ v[i][2] ]
            for ii in range(2):
                # Contribution from root arcs through vertex v[i][ii] of
                # triangle f[i].
                newWt += self._arcCoords[ f[i] ][ 3 + v[i][ii] ]

        # Contributions of the third type.
        maxArcs = 0
        maxi = 0
        maxii = 0
        for i in range(2):
            for ii in range(2):
                arcs = self._arcCoords[ f[i] ][ v[i][ii] ]
                if arcs > maxArcs:
                    maxArcs = arcs
                    maxi = i
                    maxii = ii
        other = 1 - maxi
        otherNormal = self._arcCoords[ f[other] ][ v[other][maxii] ]
        otherRoot = self._arcCoords[ f[other] ][ 3 + v[other][2] ]
        contribution = maxArcs - otherNormal - otherRoot
        if contribution > 0:
            newWt += contribution
        return newWt

    def flipWeightChange( self, e ):
        """
        Returns the change in total weight that would result from flipping
        the edge e.

        Pre-condition;
        --> e is a boundary edge of self.triangulation().
        """
        return ( self.flipWeight(e) - self._weights[ e.index() ] )

    def isReducible( self, e ):
        """
        Is the given edge e reducible?

        Pre-condition:
        --> e is a boundary edge of self.triangulation().
        """
        return ( self.flipWeight(e) < self._weights[ e.index() ] )

    def isResolvable( self, e ):
        """
        Is the given edge e "resolvable" in the sense that flipping e would
        create a new edge that forms a resolved petal?

        This routine could detect two Heegaard petals that are isotopic; if
        it does detect such petals, then it raises ValueError.

        Pre-condition:
        --> e is a boundary edge of self.triangulation().
        """
        f, v = _faceNumberings(e)

        # Triangles f[0] and f[1] must both contain root arcs that meet e.
        root = []
        for i in range(2):
            root.append( self._arcCoords[ f[i] ][ 3 + v[i][2] ] )
            if root[i] == 0:
                return False

        # Check whether we have a pair of root arcs that "match up" to give a
        # petal that would become resolved after flipping the edge e.
        normal = []
        for i in range(2):
            n = []
            for ii in range(2):
                n.append( self._arcCoords[ f[i] ][ v[i][ii] ] )
            normal.append(n)
        maxArcs = 0
        maxi = 0
        maxii = 0
        for i in range(2):
            for ii in range(2):
                if normal[i][ii] > maxArcs:
                    maxArcs = normal[i][ii]
                    maxi = i
                    maxii = ii
        other = 1 - maxi
        leftoverRoots = ( normal[other][maxii] + root[other]
                - normal[maxi][maxii] )
        if leftoverRoots < 1:
            return False

        # We can only have exactly one petal p that would become resolved
        # after flipping the edge e; this is because any other such petal
        # would be isotopic to p, which would violate admissibility.
        newResolvedPetals = min( leftoverRoots, root[maxi] )
        if newResolvedPetals == 1:
            return True
        else:
            raise ValueError( "Edge {} meets ".format( e.index() ) +
                    "{} Heegaard petals ".format(newResolvedPetals) +
                    "that are isotopic to each other." )

    def _flipEdgeImpl( self, e, edgeRefs=None ):
        """
        Flips the given boundary edge e and returns the new boundary edge
        that results from this flip.

        If e is not a resolved Heegaard petal, then this routine
        automatically updates the stored Heegaard bouquet so that it remains
        topologically equivalent to the Heegaard bouquet that we had before
        the flip. However, if e *is* a resolved petal, then this routine
        removes e from the set of resolved petals; the intention is that the
        only time we should ever flip a resolved edge is when we want to
        attach a disc by folding along this edge.

        If flipping the edge e causes every Heegaard petal to become
        resolved, then this routine conclusively checks that the following
        conditions are satisfied:
        (a) There are no normal curves.
        (b) The Heegaard bouquet is (topologically) admissible.
        This routine raises ValueError if either of these conditions fails.

        If edgeRefs is None (the default), then this routine has no
        additional functionality.

        However, if edgeRefs is not None, then it must be a 2-element list
        consisting of the following items:
        (0) A list F of edge-embeddings of distinct boundary edges of
            self.triangulation(). (The intention is that F is a stack of
            edges that we wish to flip one by one, which requires that we
            keep track of these edges as we perform flips.)
        (1) A dictionary D that maps tetrahedron indices i to the set of
            locations in the list F that reference tetrahedron i of
            self.triangulation(). We allow D to reference locations that no
            longer exist in F because the corresponding element of F was
            popped from the end of F at some point.
        In this case, this routine provides the additional functionality of
        adjusting these two items subject to the following conditions:
        --> Each edge-embedding in F represents the same edge as before, but
            we may have changed the actual edge-embedding to ensure that we
            are still referencing a tetrahedron that actually exists in
            self.triangulation() (this is necessary because flipping the edge
            e sometimes involves removing a tetrahedron).
        --> The dictionary D is adjusted so that the correspondence between
            tetrahedron indices and locations in F still holds.

        Pre-condition:
        --> e is a boundary edge of self.triangulation().
        --> If edgeRefs is not None, then it must be a 2-element list whose
            elements satisfy the description given above.
        """
        removeInfo = _recogniseLayer(e)
        if removeInfo is None:
            doomed = None
        else:
            doomed, newEmb = removeInfo
            newTet = newEmb.tetrahedron()
            newEdgeNum = newEmb.edge()

        # If we are going to perform this flip by removing a "doomed"
        # tetrahedron, and if we are given some edgeRefs that need to be
        # adjusted to account for this, then perform these adjustments now.
        if ( doomed is not None ) and ( edgeRefs is not None ):
            _adjustEdgeReferences( edgeRefs, doomed )

        # Before we perform the flip, we need to record the information
        # needed to update the edge weights and the resolved edge indices.
        updateInfo = []
        for edge in self._tri.boundaryComponent(0).edges():
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
            updateInfo.append( (
                tet,
                edgeNum,
                self._weights[ edge.index() ],
                edge.index() ) )
        flipWeight = self.flipWeight(e)
        isResolvable = self.isResolvable(e)

        # We are now ready to perform the flip. Do this by either removing a
        # doomed tetrahedron, or layering on a new tetrahedron.
        if doomed is None:
            newTet = self._tri.layerOn(e)
            newEdgeNum = 5
        else:
            self._tri.removeTetrahedron(doomed)
        newEdge = newTet.edge(newEdgeNum)

        # Compute new edge weights and new resolved edge indices.
        newWeights = [0] * self._tri.countEdges()
        updateInfo.append( ( newTet, newEdgeNum, flipWeight, None ) )
        newResEdges = set()
        if isResolvable:
            newResEdges.add( newEdge.index() )
        for tet, edgeNum, wt, oldInd in updateInfo:
            newInd = tet.edge(edgeNum).index()
            newWeights[newInd] = wt
            if oldInd in self._resolved:
                newResEdges.add(newInd)

        # Update the Heegaard bouquet. Note that self._setBouquetImpl() will
        # check normal curves and admissibility if necessary.
        arcCoords, arcInfo = _computeArcCoords( self._tri, newWeights )
        self._setBouquetImpl(
                self._tri, newWeights, newResEdges, arcCoords, arcInfo )
        return newEdge

    def flipEdge( self, e ):
        """
        If the given edge e is a boundary edge that does not form a resolved
        petal, then flips e and returns the new boundary edge that results
        from this flip; otherwise, simply returns None.

        This routine automatically updates the stored Heegaard bouquet so
        that after the flip, it remains topologically equivalent to the
        currently stored Heegaard bouquet.

        If e happens to be the new boundary edge that results from layering a
        tetrahedron t, then we flip e by simply removing t. Otherwise, we
        flip e by layering a new tetrahedron across e.

        If flipping the edge e causes every Heegaard petal to become
        resolved, then this routine conclusively checks that the following
        conditions are satisfied:
        (a) There are no normal curves.
        (b) The Heegaard bouquet is (topologically) admissible.
        This routine raises ValueError if either of these conditions fails.

        Pre-condition:
        --> e is an edge of self.triangulation().
        """
        if ( not e.isBoundary() ) or ( e.index() in self._resolved ):
            return None
        else:
            return self._flipEdgeImpl(e)

    def reducibleEdges(self):
        """
        Yields all reducible edges.

        This routine raises ValueError if this HeegaardBuilder is currently
        empty.
        """
        if self.isEmpty():
            raise ValueError( "Empty HeegaardBuilder." )
        for e in self._tri.boundaryComponent(0).edges():
            if self.isReducible(e):
                yield e

    def resolveAllPetals(self):
        """
        Flips edges until all of the Heegaard petals are resolved.

        This routine chooses reducible edges arbitrarily.

        This routine raises ValueError if:
        --> this HeegaardBuilder is currently empty; or
        --> while flipping edges, we reach a situation where there are no
            reducible edges.
        """
        if self.isEmpty():
            raise ValueError( "Empty HeegaardBuilder." )
        while self.countResolved() < self._genus:
            reduced = False
            for e in self.reducibleEdges():
                # If flipping the edge e would leave behind only resolved
                # petals, then flipEdge() will check that:
                # - there are no normal arcs; and
                # - the Heegaard bouquet is (topologically) admissible.
                self.flipEdge(e)
                reduced = True
                break
            if not reduced:
                # Since we can always reduce when we have a Heegaard bouquet
                # (regardless of whether it is admissible), the only way we
                # could run out of reducible edges is if we have a normal
                # curve.
                raise ValueError( "Could not find a reducible edge, " +
                        "which means there must be a normal curve." )

    def resolveGreedily(self):
        """
        Greedily flips edges until all of the Heegaard petals are resolved.

        The word "greedily" refers to the method by which this routine
        chooses the edge to flip: it always chooses an edge that reduces the
        total edge weight as much as possible.

        This routine raises ValueError if:
        --> this HeegaardBuilder is currently empty; or
        --> while flipping edges, we reach a situation where there are no
            reducible edges.
        """
        if self.isEmpty():
            raise ValueError( "Empty HeegaardBuilder." )
        while self.countResolved() < self._genus:
            e = None
            reduction = 0
            for candidate in self.reducibleEdges():
                r = self.flipWeightChange(candidate)
                if r < reduction:
                    reduction = r
                    e = candidate
            if e is None:
                # Since we can always reduce when we have a Heegaard bouquet
                # (regardless of whether it is admissible), the only way we
                # could run out of reducible edges is if we have a normal
                # curve.
                raise ValueError( "Could not find a reducible edge, " +
                        "which means there must be a normal curve." )
            # If flipping the edge e would leave behind only resolved
            # petals, then flipEdge() will check that:
            # - there are no normal arcs; and
            # - the Heegaard bouquet is (topologically) admissible.
            self.flipEdge(e)

    def resolveUntilChoice( self, greedy=True ):
        """
        Attempts to flip edges until all of the Heegaard petals are resolved,
        but stops when we reach a situation where we have multiple choices
        for which edge to flip.
        
        If greedy is True (the default), then this routine only considers
        reducible edges that reduce the total edge weight as much as
        possible. You can force this routine to consider *all* reducible
        edges by setting greedy to False, but be warned that this could
        significantly increase the number of available choices.

        This routine returns None if it successfully resolves all Heegaard
        petals without ever having to arbitrarily choose among multiple
        reducible edges. However, if such a choice is required, then this
        routine stops flipping edges, and instead returns a list containing
        all of the reducible edges that could have been flipped.

        This routine raises ValueError if:
        --> this HeegaardBuilder is currently empty; or
        --> while flipping edges, we reach a situation where there are no
            reducible edges.
        """
        if self.isEmpty():
            raise ValueError( "Empty HeegaardBuilder." )
        while self.countResolved() < self._genus:
            reduction = 0 # Only used if greedy is True.
            redEdges = []
            for e in self.reducibleEdges():
                if greedy:
                    change = self.flipWeightChange(e)
                    if change > reduction:
                        continue
                    if change < reduction:
                        reduction = change
                        redEdges = []
                redEdges.append(e)
            if redEdges:
                # We found at least one reducible edge.
                if len(redEdges) == 1:
                    # If flipping the edge e would leave behind only resolved
                    # petals, then flipEdge() will check that:
                    # - there are no normal arcs; and
                    # - the Heegaard bouquet is (topologically) admissible.
                    self.flipEdge( redEdges[0] )
                else:
                    return redEdges
            else:
                # Since we can always reduce when we have a Heegaard bouquet
                # (regardless of whether it is admissible), the only way we
                # could run out of reducible edges is if we have a normal
                # curve.
                raise ValueError( "Could not find a reducible edge, " +
                        "which means there must be a normal curve." )

        # If we exit the above loop, then we successfully resolved all petals
        # without needing to make an arbitrary choice.
        return None

    def resolveInAllWays( self, greedy=True ):
        """
        Yields instances of HeegaardBuilder corresponding to all possible
        ways to resolve the stored Heegaard bouquet.
        
        If greedy is True (the default), then this routine only considers
        reducible edges that reduce the total edge weight as much as
        possible. You can force this routine to consider *all* reducible
        edges by setting greedy to False, but be warned that this could
        significantly increase the number of available choices.

        This routine never modifies this instance of HeegaardBuilder.

        This routine raises ValueError if:
        --> this HeegaardBuilder is currently empty; or
        --> while flipping edges, we reach a situation where there are no
            reducible edges.

        Warning:
        --> In total, this routine could yield an enormous number of
            instances of HeegaardBuilder, even if greedy is set to True.
        """
        if self.isEmpty():
            raise ValueError( "Empty HeegaardBuilder." )
        hb = HeegaardBuilder()
        hb.setClone(self)
        redEdges = hb.resolveUntilChoice(greedy)
        if redEdges is None:
            yield hb
        else:
            for e in redEdges:
                hbClone = HeegaardBuilder()
                hbClone.setClone(hb)
                # If flipping the edge e would leave behind only resolved
                # petals, then flipEdge() will check that:
                # - there are no normal arcs; and
                # - the Heegaard bouquet is (topologically) admissible.
                hbClone.flipEdge(
                        hbClone.triangulation().edge( e.index() ) )
                for r in hbClone.resolveInAllWays(greedy):
                    yield r

    def _fold( self, e ):
        """
        Folds about the given boundary edge e of self.triangulation().
        
        Such a fold is only possible if the boundary triangles on either side
        of e are distinct. This routine raises RuntimeError if the given edge
        e does not satisfy this requirement.

        This routine assumes that all edge weights are 0; if this condition
        fails, then this routine will erase all non-zero edge weights. This
        routine also assumes that none of the four boundary edges opposite e
        form resolved petals (however, e itself may form a resolved petal).

        Provided that the two boundary vertices opposite the edge e are not
        pinched together, performing the requested fold will not create any
        ideal vertices in self.triangulation().

        Pre-condition:
        --> e is a boundary edge of self.triangulation().
        --> All edge weights are 0.
        --> None of the four boundary edges opposite e form resolved petals.
        """
        tet = []
        ver = []
        for i in [ 0, e.degree() - 1 ]:
            emb = e.embedding(i)
            tet.append( emb.tetrahedron() )
            ver.append( emb.vertices() )
        if tet[0] == tet[1] and ver[0][3] == ver[1][2]:
            # Cannot fold if the boundary triangles on either side of the
            # edge e are actually the same triangle.
            raise RuntimeError( "Attempted to fold a triangle " +
                    "onto itself." )
        else:
            tet[0].join( ver[0][3], tet[1],
                    ver[1] * Perm4(2,3) * ver[0].inverse() )

            # Update the Heegaard bouquet.
            wts = [0] * self.triangulation().countEdges()
            self._resolved.discard( e.index() )
            arcCoords, arcInfo = _computeArcCoords( self._tri, wts )
            self._setBouquetImpl(
                    self._tri, wts, self._resolved, arcCoords, arcInfo )

    def _attemptFold(self):
        """
        Assuming that boundary component 0 of self.triangulation() is a
        2-sphere, attempts to either partially or completely fill this
        2-sphere with a 3-ball by folding a pair of boundary faces together.

        Returns True if and only if such a fold is successfully performed.

        Pre-condition:
        --> This HeegaardBuilder is currently non-empty.
        --> The surface self.triangulation().boundaryComponent(0) is a
            2-sphere.
        """
        # The idea is to iterate through each boundary edge e, and consider
        # the square formed by the two boundary triangles incident to e. We
        # could potentially fold this square across either of its two
        # diagonals (for one of these folds, we would first need to perform
        # an edge flip), so we check whether these folds are "safe".
        bc = self._tri.boundaryComponent(0)
        built = bc.build()
        flip = None # Edge that we can fold after flipping, if necessary.
        for e in built.edges():
            face = []
            ver = []
            for emb in e.embeddings():
                face.append( emb.triangle() )
                ver.append( emb.vertices() )
            if face[0] == face[1]:
                continue

            # First check whether it is safe to fold directly over e.
            if face[0].vertex( ver[0][2] ) == face[1].vertex( ver[1][2] ):
                # The fold is safe, so perform it.
                self._fold( bc.triangle( face[0].index() ).edge( ver[0][2] ) )
                return True

            # Now, if necessary, check whether it would be safe to fold after
            # flipping e.
            if ( flip is None ) and ( e.vertex(0) == e.vertex(1) ):
                # Only perform this fold later on, if it is really necessary.
                flip = bc.triangle( face[0].index() ).edge( ver[0][2] )

        # We couldn't directly fold, but maybe we can fold after flipping.
        if flip is None:
            return None
        else:
            # It is safe to fold after flipping.
            self._fold( self._flipEdgeImpl(flip) )
            return True

    def _fillBall(self):
        """
        Assuming that self.triangulation() has 2-sphere boundary, fills this
        2-sphere with a 3-ball.

        Pre-condition:
        --> This HeegaardBuilder is currently non-empty.
        --> The boundary surface of self.triangulation() is a 2-sphere.
        """
        while self._tri.hasBoundaryTriangles():
            if self._attemptFold():
                continue

            # Since folding was not possible, we know (in particular) that
            # the boundary surface S has no vertices of degree less than 3.
            # Our goal now is to keep flipping edges until it *is* possible
            # to perform a fold; provided we take care to always reduce the
            # minimum vertex degree in S, this procedure must eventually
            # terminate. Thus, we start by finding a vertex of minimum degree
            # in S.
            bc = self._tri.boundaryComponent(0)
            built = bc.build()
            minVert = built.vertex(0)
            for i in range( 1, built.countVertices() ):
                v = built.vertex(i)
                deg = v.degree()
                if deg < minVert.degree():
                    minVert = v

                    # We know that built contains no vertices of degree less
                    # than 3, so if we have reduced minVert.degree() down to
                    # 3, then we have already found the minimum.
                    if deg == 3:
                        break

            # Create a stack of edge-embeddings such that flipping the
            # corresponding edges will reduce the degree of minVert.
            flip = []
            refTetInds = dict()
            edgeRefs = [ flip, refTetInds ]
            for i in range( minVert.degree() - 2 ):
                vertEmb = minVert.embedding(i)
                faceInd = vertEmb.triangle().index()
                edgeNum = vertEmb.vertices()[2]
                _addEdgeReference(
                        bc.triangle(faceInd).edge(edgeNum), edgeRefs )

            # Keep flipping until we can fold.
            while True:
                emb = flip.pop()
                self._flipEdgeImpl(
                        emb.tetrahedron().edge( emb.edge() ), edgeRefs )
                if self._attemptFold():
                    break

    def fillHandlebody(self):
        """
        If every Heegaard petal is resolved, then fills self.triangulation()
        with a handlebody H such that each Heegaard petal bounds a disc in H.

        If every Heegaard petal is indeed resolved, then this routine also
        returns True (after it has performed the operation of filling with a
        handlebody). Otherwise, this routine returns False and leaves
        self.triangulation() unchanged.

        This routine raises ValueError if this HeegaardBuilder is currently
        empty.
        """
        if self.isEmpty():
            raise ValueError( "Empty HeegaardBuilder." )
        if self.countResolved() != self._genus:
            return False

        # We know that every Heegaard petal is resolved. Our goal is to
        # attach a disc to each petal by flipping and then folding. However,
        # if this requires us to flip two edges of a single triangle f, then
        # we first need to split f into two triangles by flipping the third
        # edge of f. We use list of edge-embeddings to keep track of all the
        # edges that we need to flip, and we treat this list as a stack to
        # keep track of the order in which we need to flip the edges.
        flipStack = []
        refTetInds = dict()
        edgeRefs = [ flipStack, refTetInds ]
        flipSet = set()
        for e in self._resolved:
            flipSet.add(e)
            _addEdgeReference( self._tri.edge(e), edgeRefs )
        splitFace = True
        splitIndices = set()
        while splitFace:
            splitFace = False

            # For each boundary face f that we have not yet split, check
            # whether it has become necessary to split f.
            for f in self._tri.boundaryComponent(0).triangles():
                if f.index() in splitIndices:
                    continue

                # If we need to flip two of the edges of f, then we need to
                # split f by flipping its third edge.
                unflippedEdgeNums = {0,1,2}
                for e in range(3):
                    if f.edge(e).index() in flipSet:
                        unflippedEdgeNums.remove(e)
                if len(unflippedEdgeNums) == 1:
                    edge = f.edge( unflippedEdgeNums.pop() )
                    flipSet.add( edge.index() )
                    _addEdgeReference( edge, edgeRefs )
                    splitFace = True
                    splitIndices.add( f.index() )

        # Flip every edge in flipStack. Since we know that the last few flips
        # correspond to the resolved petals, we can also fold these now.
        foldCount = self.countResolved()
        while flipStack:
            emb = flipStack.pop()
            newEdge = self._flipEdgeImpl(
                    emb.tetrahedron().edge( emb.edge() ), edgeRefs )
            if len(flipStack) < foldCount:
                self._fold(newEdge)

        # We should now be left with sphere boundary. Complete the
        # construction by filling with a ball.
        self._fillBall()
        return True
