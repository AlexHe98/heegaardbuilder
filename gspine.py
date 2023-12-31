"""
Find g-spines in a triangulation.
"""
from regina import *


def _nextCombination(c):
    """
    Given a list c of distinct non-negative integers arranged in ascending
    order, returns the "next" combination.

    In detail, the "next" combination is defined as follows:
    --> If c == [ 0, ..., len(c) - 1 ], then the next combination is None.
    --> Otherwise, let i be the smallest index such that c[i] != i. In this
        case, the next combination is given by
            [ c[i] - (i+1), ..., c[i] - 1, c[i+1], ..., c[n-1] ].
        In this case, this routine modifies c in-place, and returns c after
        these modifications have been completed.
    """
    # Find i such that c[i] > i.
    found = False
    for i, entry in enumerate(c):
        if entry > i:
            found = True
            break
    if not found:
        return None

    # Generate next combination.
    newEntry = entry - i - 1
    for ii in range(i+1):
        c[ii] = newEntry
        newEntry += 1
    return c


def combinations( n, k ):
    """
    Yields all combinations of k elements from the set { 0, ..., n-1 }.
    """
    c = [ n-k+i for i in range(k) ]
    while True:
        yield c
        if _nextCombination(c) is None:
            return


#TODO Add functionality to extend beyond minimal complete spine layerings
#   (which we need to extend beyond minimal layered handlebodies).
class SpineLayering:
    #TODO More detailed class notes. In particular, should define precisely
    #   what constitutes an "extension" of a spine layering.
    """
    Represents a layering of tetrahedra onto a "spine" in the 2-skeleton of a
    3-manifold triangulation.
    """
    def __init__( self, spineFaces,
            _spine=None, _layeredTetInds=None, _bdryFaceInds=None,
            _tetIncidences=None, _unlayeredEdges=None, _ungluedEdges=None ):
        """
        Initialises a trivial (zero-tetrahedron) layering onto the spine
        formed by the given triangular faces.

        Pre-condition:
        --> spineFaces must be a collection of distinct triangular faces of a
            single 3-manifold triangulation.

        Warning:
        --> Modifying the triangulation that contains the given triangular
            faces will invalidate this spine layering, as well as all
            extensions of this spine layering.

        Internal-use parameters:
            The parameters _spine, _layeredTetInds, _bdryFaceInds,
            _tetIncidences, _unlayeredEdges and _ungluedEdges are for
            internal use only. End-users should never change these parameters
            from their default values.
        """
        self._spineFaces = list(spineFaces)
        self._spineSize = len(spineFaces)
        self._tri = spineFaces[0].triangulation()
        self._detail = None
        if _spine is None:
            # Set this spine layering to be trivial.
            self._spine = Triangulation2()
            self._spine.newTriangles( self._spineSize )
            self._layeredTetInds = []
            self._bdryFaceInds = dict()
            self._tetIncidences = dict()
            self._unlayeredEdges = dict()
            self._ungluedEdges = dict()
            for i in range( self._spineSize ):
                ###
                # self._bdryFaceInds: For each boundary face f, maps the
                #   index of f to the number of tetrahedra incident to f.
                face = self._spineFaces[i]
                faceInd = face.index()
                self._bdryFaceInds[faceInd] = 2
                # self._tetIncidences: For each tetrahedron t incident to the
                #   boundary faces of this spine layering, maps the index of
                #   t to a set of pairs of the form ( face index, embedding
                #   index ); each such pair describes one of the ways in
                #   which t is incident to this spine layering.
                for embInd in range(2):
                    if embInd == 1 and face.isBoundary():
                        continue
                    incidence = ( faceInd, embInd )
                    tetInd = face.embedding(embInd).tetrahedron().index()
                    if tetInd in self._tetIncidences:
                        self._tetIncidences[tetInd].add(incidence)
                    else:
                        self._tetIncidences[tetInd] = {incidence}
                # self._unlayeredEdges and self._ungluedEdges: These mappings
                #   keep track of which edges are available for layering.
                for direction in range(2):
                    for edgeNum in range(3):
                        unlayered = ( faceInd, direction, edgeNum )
                        self._unlayeredEdges[unlayered] = (
                                i, direction, Perm4() )
                        unglued = ( i, direction, edgeNum )
                        self._ungluedEdges[unglued] = unlayered
        else:
            # Internal use: Set a non-trivial spine when constructing
            #   extensions.
            self._spine = _spine
            self._layeredTetInds = _layeredTetInds
            self._bdryFaceInds = _bdryFaceInds
            self._tetIncidences = _tetIncidences
            self._unlayeredEdges = _unlayeredEdges
            self._ungluedEdges = _ungluedEdges

    def triangulation(self):
        """
        Returns the triangulation that contains this spine layering.

        Warning:
        --> Modifying the returned triangulation will invalidate this spine
            layering, as well as all extensions of this spine layering.
        """
        return self._tri

    def size(self):
        """
        Returns the number of tetrahedra in this spine layering.
        """
        return len( self._layeredTetInds )

    def layeredTetrahedron( self, n ):
        """
        Returns the nth tetrahedron in this spine layering.

        Pre-condition:
        --> n is a non-negative integer strictly less than self.size().
        """
        return self._tri.tetrahedron( self._layeredTetInds[n] )

    def layeredTriangulation(self):
        """
        Constructs and returns the 3-dimensional triangulation described by
        this spine layering, together with a mapping from tetrahedron indices
        in the returned triangulation to tetrahedron indices in the
        triangulation that contains this spine layering.
        """
        layeredTri = Triangulation3( self._tri )

        # Deleting tetrahedra is not sufficient to correctly construct
        # layeredTri, because self._tri might have gluings between faces that
        # we consider to be in the boundary of the layering. Thus, we first
        # unglue all faces that are supposed to be boundary in layeredTri.
        for b in self._bdryFaceInds:
            emb = self._tri.triangle(b).embedding(0)
            tetIndex = emb.tetrahedron().index()
            faceNum = emb.face()
            layeredTri.tetrahedron(tetIndex).unjoin(faceNum)

        # Now that we have fixed the boundary, all we need to do is delete
        # the tetrahedra that are not part of the layering.
        mapping = [ i for i in range( self._tri.size() ) ]
        for i in range( self._tri.size() - 1, -1, -1 ):
            if i not in self._layeredTetInds:
                layeredTri.removeTetrahedronAt(i)
                mapping.remove(i)
        return ( layeredTri,
                { i: mapping[i] for i in range( layeredTri.size() ) } )

    def spineSize(self):
        """
        Returns the number of triangles in the underlying spine.
        """
        return self._spineSize

    def spineFace( self, index ):
        """
        Returns the triangular face at the given index in the underlying
        spine.
        """
        return self._spineFaces[index]

    def spineTriangulation(self):
        """
        Returns the 2-manifold triangulation corresponding to the underlying
        spine.

        In detail, the returned triangulation has the following
        correspondence with the underlying spine:
        --> Let i be a non-negative integer strictly less than
            self.spineSize(). Then triangle i of the returned triangulation
            is a copy of the triangle given by self.spineFace(i), and the
            vertices of these two triangles are numbered in the same way.
        --> Two triangles in the returned triangulation are glued together
            along an edge e if and only if this gluing is "witnessed" by a
            tetrahedron in this spine layering that is layered across the
            edge corresponding to e.

        Warning:
        --> Modifying the returned triangulation will invalidate this spine
            layering.
        """
        return self._spine

    def isComplete(self):
        """
        Is this spine layering complete?

        In detail, we consider this spine layering to be complete if the
        2-manifold triangulation S corresponding to the underlying spine
        satisfies the following properties:
        (1) S is connected; and
        (2) S has at most one unglued edge.
        This routine returns True if and only if these two properties hold.
        """
        return ( self._spine.isConnected and
                len( self._ungluedEdges ) < 4 )

    def _extensionsByPair( self, tetInd, faceInd, embInd, pair ):
        # Extract data about the pair of faces on which we would like to
        # extend this spine layering.
        tet = self._tri.tetrahedron(tetInd)
        ver = []
        for i in range(2):
            face = self._tri.triangle( faceInd[ pair[i] ] )
            emb = face.embedding( embInd[ pair[i] ] )
            ver.append( emb.vertices() )
        edgeNum = []
        unlayered = []
        for i in range(2):
            edgeNum.append( ver[i].inverse()[ ver[1-i][3] ] )
            unlayered.append( (
                faceInd[ pair[i] ], embInd[ pair[i] ], edgeNum[i] ) )

        # We can extend this spine layering if the given pair of incidences
        # "witnesses" a gluing between two unglued edges in the spine.
        witnessedFaceInd = []
        witnessedPerm = []
        witnessedEdgeNum = []
        for i in range(2):
            if unlayered[i] in self._unlayeredEdges:
                wFace, _, wPerm = self._unlayeredEdges[ unlayered[i] ]
                witnessedFaceInd.append(wFace)
                witnessedPerm.append(wPerm)
                witnessedEdgeNum.append( wPerm[ edgeNum[i] ] )
            else:
                break
        if len(witnessedFaceInd) != 2:
            return
        elif ( witnessedFaceInd[0] == witnessedFaceInd[1] and
                witnessedEdgeNum[0] == witnessedEdgeNum[1] ):
            return

        # Make copies of the data from this spine layering, so that we can
        # compute the extended data without modifying this spine layering.
        spine = Triangulation2( self._spine )
        layeredTetInds = list( self._layeredTetInds )
        bdryFaceInds = dict( self._bdryFaceInds )
        tetIncidences = dict()
        for t in self._tetIncidences:
            if t == tetInd:
                # This saves us the work of removing unwanted
                # incidences later on.
                continue
            tetIncidences[t] = set( self._tetIncidences[t] )
        unlayeredEdges = dict( self._unlayeredEdges )
        ungluedEdges = dict( self._ungluedEdges )

        # Compute extended spine: Perform the gluing that is witnessed by the
        #   new layered tetrahedron.
        myFace = spine.triangle( witnessedFaceInd[0] )
        myEdge = witnessedEdgeNum[0]
        yourFace = spine.triangle( witnessedFaceInd[1] )
        gluing = Perm3.contract(
                witnessedPerm[1] *
                ver[1].inverse() *
                Perm4( ver[0][3], ver[1][3] ) *
                ver[0] *
                witnessedPerm[0].inverse() )
        myFace.join( myEdge, yourFace, gluing )

        # Compute extended layeredTetInds: Append the index of the new
        #   layered tetrahedron.
        layeredTetInds.append(tetInd)

        # Compute extended bdryFaceInds: Remove the two faces on which we
        #   layered, and add the two new boundary faces given by the new
        #   layered tetrahedron.
        # We also take this opportunity to collect some data that is
        # necessary for extending:
        #   --> tetIncidences
        #   --> unlayeredEdges
        #   --> ungluedEdges
        oppTetInd = []          # For tetIncidences.
        oppIncidence = []       # For tetIncidences.
        oppDirection = dict()   # For unlayeredEdges and ungluedEdges.
        layerFaceNums = { ver[0][3], ver[1][3] }
        otherFaceNums = {0,1,2,3} - layerFaceNums
        for faceNum in otherFaceNums:
            # After extending this spine layering, tet.triangle(faceNum)
            # becomes a new boundary triangle.
            bf = tet.triangle(faceNum)
            bfInd = bf.index()
            if bfInd in bdryFaceInds:
                bdryFaceInds[bfInd] += 1
            else:
                bdryFaceInds[bfInd] = 1

            # Populate oppTetInd, oppIncidence and oppDirection.
            if bf.isBoundary():
                # Boundary face contributes no incidence.
                oppTetInd.append(None)
                oppIncidence.append(None)
                oppDirection[faceNum] = 1
            else:
                # For an internal face, find the embedding that is opposite
                # the new layered tetrahedron.
                for oppEmbInd in range(2):
                    oppEmb = bf.embedding(oppEmbInd)
                    oppTet = oppEmb.tetrahedron()
                    oppFaceNum = oppEmb.face()
                    if oppTet == tet and oppFaceNum == faceNum:
                        # This is not the opposite embedding.
                        continue
                    elif oppTet.index() in layeredTetInds:
                        # The opposite tetrahedron is already part of this
                        # spine layering, so it doesn't contribute a new
                        # tetrahedron incidence.
                        oppTetInd.append(None)
                        oppIncidence.append(None)
                        oppDirection[faceNum] = oppEmbInd
                        break
                    else:
                        # The opposite tetrahedron contributes a new
                        # tetrahedron incidence.
                        oppTetInd.append( oppTet.index() )
                        oppIncidence.append( ( bfInd, oppEmbInd ) )
                        oppDirection[faceNum] = oppEmbInd
                        break
        for i in range(2):
            # Remove (one copy of) faceInd[ pair[i] ] from bdryFaceInds.
            f = faceInd[ pair[i] ]
            bdryFaceInds[f] -= 1
            if bdryFaceInds[f] == 0:
                bdryFaceInds.pop(f)

        # Compute extended tetIncidences: Add incidences given by the two new
        #   boundary faces.
        # Note that we do not need to remove incidences to the new layered
        # tetrahedron, because we already handled this when we initialised
        # tetIncidences, oppTetInd and oppIncidence.
        for i in range(2):
            if oppTetInd[i] is None:
                continue
            elif oppTetInd[i] in tetIncidences:
                tetIncidences[ oppTetInd[i] ].add( oppIncidence[i] )
            else:
                tetIncidences[ oppTetInd[i] ] = { oppIncidence[i] }

        # Compute extended unlayeredEdges and ungluedEdges: Completely remove
        #   unlayered edge data associated to the edge on which we just
        #   layered, and replace any other unlayered edge data associated to
        #   the layered faces with unlayered edge data associated to the new
        #   boundary faces.
        for i in range(2):
            for direction in range(2):
                # Remove unlayered edge data associated to the edge on which
                # we just layered.
                glued = ( witnessedFaceInd[i],
                        direction,
                        witnessedEdgeNum[i] )
                layered = ungluedEdges.pop(glued)
                unlayeredEdges.pop(layered)
        for faceNum in otherFaceNums:
            newFaceInd = tet.triangle(faceNum).index()
            newDirection = oppDirection[faceNum]
            fMap = tet.triangleMapping(faceNum)
            for i in range(2):
                # If there is old unlayered edge data, then we need to
                # replace it with new data.
                oldFaceInd = faceInd[ pair[i] ]
                oldDirection = embInd[ pair[i] ]
                oldEdgeNum = ver[i].inverse()[faceNum]
                oldUnlayered = ( oldFaceInd, oldDirection, oldEdgeNum )
                if oldUnlayered not in unlayeredEdges:
                    # There is no old data to update.
                    continue

                # There is old data that we need to update.
                # First update unlayeredEdges.
                newEdgeNum = fMap.inverse()[ ver[i][3] ]
                newUnlayered = ( newFaceInd, newDirection, newEdgeNum )
                spineInd, spineDir, oldPerm = unlayeredEdges.pop(
                        oldUnlayered )
                newPerm = ( oldPerm *
                        ver[i].inverse() *
                        Perm4( ver[i][3], faceNum ) *
                        fMap )
                unlayeredEdges[newUnlayered] = (
                        spineInd, spineDir, newPerm )

                # Now update ungluedEdges.
                unglued = ( spineInd, spineDir, oldPerm[oldEdgeNum] )
                ungluedEdges[unglued] = newUnlayered

        # Yield extended spine layering.
        yield SpineLayering( self._spineFaces,
                spine, layeredTetInds, bdryFaceInds, tetIncidences,
                unlayeredEdges, ungluedEdges )

    def _extensionsByTetInd( self, tetInd ):
        # We need two or more incidences to be able to extend.
        incidences = self._tetIncidences[tetInd]
        count = len(incidences)
        if count < 2:
            return

        # Extract indices for the faces and embeddings.
        faceInd = []
        embInd = []
        for fi, ei in incidences:
            faceInd.append(fi)
            embInd.append(ei)

        # Each pair of incidences potentially gives a way to extend.
        for pair in combinations( count, 2 ):
            for extension in self._extensionsByPair(
                    tetInd, faceInd, embInd, pair ):
                yield extension

    def extensionsByOne(self):
        """
        Yields all spine layerings that can be constructed by extending this
        spine layering by one tetrahedron.
        """
        #TODO Allow extensions beyond minimal complete layerings.
        if self.isComplete():
            return

        # For each tetrahedron that is incident to this spine layering, check
        # whether we can extend this spine layering through this tetrahedron.
        for tetInd in self._tetIncidences:
            for extension in self._extensionsByTetInd(tetInd):
                yield extension

    def minimalCompleteExtensions(self):
        """
        Yields all minimal complete spine layerings that can be constructed
        by extending this spine layering.
        """
        if self.isComplete():
            yield self
        else:
            for extensionByOne in self.extensionsByOne():
                for extension in extensionByOne.minimalCompleteExtensions():
                    yield extension

    def detail(self):
        """
        Returns a detailed human-readable description of this spine layering.
        """
        #TODO Make the description more detailed.
        if self._detail is None:
            spineFaceInds = [ f.index() for f in self._spineFaces ]
            self._detail = "Layered tetrahedra: {}\n\n".format(
                    self._layeredTetInds )
            self._detail += "Spine:\n"
            self._detail += " Triangle  |  glued to:"
            self._detail += "     (01)     (02)     (12)\n"
            self._detail += " ----------+-----------"
            self._detail += "---------------------------\n"
            for f in spineFaceInds:
                self._detail += "      {:>3}  |           ".format(f)
                myFace = self._spine.triangle( spineFaceInds.index(f) )
                for e, i, j in ( (2,0,1), (1,0,2), (0,1,2) ):
                    yourFace = myFace.adjacentTriangle(e)
                    if yourFace is None:
                        self._detail += "  unglued"
                    else:
                        gluing = myFace.adjacentGluing(e)
                        self._detail += " {:>3} ({}{})".format(
                                spineFaceInds[ yourFace.index() ],
                                gluing[i], gluing[j] )
                self._detail += "\n"
            self._detail += "\n"
        return self._detail


def _trivialLayerings( tri, k ):
    for spineFaceInds in combinations( tri.countTriangles(), k ):
        spineFaces = [ tri.triangle(i) for i in spineFaceInds ]
        yield SpineLayering(spineFaces)


def spines( tri, k ):
    """
    For each spine consisting of k triangles in the given triangulation,
    yields one possible layering on this spine.
    """
    for layering in _trivialLayerings( tri, k ):
        for extension in layering.minimalCompleteExtensions():
            yield extension
            break


def spineLayerings( tri, k ):
    """
    Yields all minimal complete spine layerings with k spine faces in the
    given triangulation.
    """
    for layering in _trivialLayerings( tri, k ):
        for extension in layering.minimalCompleteExtensions():
            yield extension
