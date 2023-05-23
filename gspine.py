"""
Find g-spines in a triangulation.
"""
from regina import *


def nextCombination(c):
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
        if nextCombination(c) is None:
            return


def surfaceFromFaces( tri, faceIndices ):
    """
    If the given list of face indices gives a surface inside tri, then
    returns a triangulation of this surface; otherwise, returns None.
    """
    # Try to extract edge gluing information from the given face indices.
    edgeGluings = dict()
    for i, f in enumerate(faceIndices):
        face = tri.triangle(f)
        for e in range(3):
            edgeInd = face.edge(e).index()
            if edgeInd in edgeGluings:
                edgeGluings[edgeInd].append(
                        ( i, e, Perm3.contract( face.edgeMapping(e) ) ) )
            else:
                edgeGluings[edgeInd] = [
                        ( i, e, Perm3.contract( face.edgeMapping(e) ) ) ]

    # Try to build a triangulation from the edge gluings.
    surf = Triangulation2()
    size = len(faceIndices)
    triangle = surf.newTriangles(size)
    for g in edgeGluings.values():
        if len(g) > 2:
            return None
        if len(g) < 2:
            continue

        # Glue two edges together in surf.
        myIndex, myEdge, myMapping = g[0]
        yourIndex, yourEdge, yourMapping = g[1]
        triangle[myIndex].join( myEdge, triangle[yourIndex],
                yourMapping * myMapping.inverse() )
    return surf


def surfacesIn2Skeleton( tri, k ):
    """
    Yields all connected surfaces in tri that are built from k triangular
    faces of tri.
    """
    for c in combinations( tri.countTriangles(), k ):
        surf = surfaceFromFaces( tri, c )
        if surf is not None:
            yield ( c, surf )


# Test code.
if __name__ == "__main__":
    testSigs = [
            "eHuGabdes",
            "iLLLPQcbddegfhghabfsccswr",
            "lLLLLPPQcadegigiihjkkjaxreousjnck",
            "mLvLLMQQPaefhikighkjlljxfrtaangkjdj",
            "oLLvzwQMAQccdhikhghjlklmnnhshsaocnhvvnwlj" ]
    for sig in testSigs:
        print()
        print( "==================================================" )
        print(sig)
        print( "--------------------------------------------------" )
        print()
        tri = Triangulation3.fromIsoSig(sig)

        # Test searching for g-spines.
        triangles = tri.countTriangles()
        surfaces = 0
        oneVertex = 0
        for c, surf in surfacesIn2Skeleton( tri, 3 ):
            surfaces += 1
            if surf.countVertices() == 1:
                oneVertex += 1
                print(c)
                print( surf.detail() )
        print()
        print( "Triangles: {}. Surfaces: {}. One-vertex: {}.".format(
            triangles, surfaces, oneVertex ) )
        print()
