"""
Small test suite for the HeegaardBuilder class.
"""
from heegaardbuilder import HeegaardBuilder
from regina import *


if __name__ == "__main__":
    print()
    hb = HeegaardBuilder()
    sigGenus2 = "eHbecadjk"
    sigGenus3 = "hHbLbqiabegeti"

    # Cases that should raise ValueError upon setting bouquet.
    print( "============================================================" )
    print( " Genus 2: ValueError upon setting bouquet" )
    print( "------------------------------------------------------------" )
    print()
    setError = [
            [ (0,0,0,0,0,0,0,0), set(), "Wrong number of weights" ],
            [ (0,0,0,0,0,0,0,-1,0), set(), "Negative weight" ],
            [ (0,1,0,0,0,0,0,0,0), {1}, "Resolved edge with weight > 0" ],
            [ (1,3,2,0,2,3,2,0,2), set(), "Matching constraints fail" ],
            [ (0,0,0,2,2,1,2,0,0), {0,2}, "Wrong number of roots" ],
            [ (0,0,0,2,2,2,3,1,1), {0,2}, "Normal curve" ],
            [ (0,0,0,0,0,0,0,0,0), {0,2}, "Transverse Heegaard petals" ] ]
    tri = Triangulation3.fromIsoSig(sigGenus2)
    for w, r, name in setError:
        print(name)
        try:
            hb.setBouquet( tri, w, r )
        except ValueError as e:
            print( "    {}".format(e) )
        else:
            raise RuntimeError( "Incorrectly set bouquet." )
        print()
    newTet = tri.layerOn( tri.edge(1) )
    intEdgeInd = newTet.edge(0).index()
    print( "Internal edge with weight > 0" )
    wList = [0] * 10
    wList[intEdgeInd] = 1
    try:
        hb.setBouquet( tri, tuple(wList), set() )
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Incorrectly set bouquet." )
    print()
    print( "Internal resolved edge" )
    w = (0,0,0,0,0,0,0,0,0,0)
    try:
        hb.setBouquet( tri, w, {intEdgeInd} )
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Incorrectly set bouquet." )
    print()
    print( "PASSED." )
    print()

    # Test empty HeegaardBuilder.
    print( "============================================================" )
    print( " Genus 2: Empty HeegaardBuilder" )
    print( "------------------------------------------------------------" )
    print()
    hb.reset()
    if not hb.isEmpty():
        raise RuntimeError(
                "HeegaardBuilder incorrectly identified as non-empty." )
    if hb.triangulation() is not None:
        raise RuntimeError(
                "Empty HeegaardBuilder contains a triangulation." )
    if hb.countResolved() is not None:
        raise RuntimeError(
                "Empty HeegaardBuilder contains resolved petals." )
    print( "HeegaardBuilder.resolveAllPetals()" )
    try:
        hb.resolveAllPetals()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to identify empty HeegaardBuilder." )
    print()
    print( "HeegaardBuilder.resolveGreedily()" )
    try:
        hb.resolveGreedily()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to identify empty HeegaardBuilder." )
    print()
    print( "HeegaardBuilder.resolveUntilChoice()" )
    try:
        hb.resolveUntilChoice()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to identify empty HeegaardBuilder." )
    print()
    print( "HeegaardBuilder.resolveInAllWays()" )
    try:
        for _ in hb.resolveInAllWays():
            pass
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to identify empty HeegaardBuilder." )
    print()
    print( "HeegaardBuilder.fillHandlebody()" )
    try:
        hb.fillHandlebody()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to identify empty HeegaardBuilder." )
    print()
    print( "PASSED." )
    print()

    # Testing edge weights and edge-flips.
    print( "============================================================" )
    print( " Genus 2: Edge weights and edge-flips" )
    print( "------------------------------------------------------------" )
    print()
    tri = Triangulation3.fromIsoSig(sigGenus2)
    w = (0,0,0,2,3,1,2,1,0)
    hb.setBouquet( tri, w )
    print( "Flip weight" )
    f = hb.flipWeight( tri.edge(4) )
    if f != 0:
        raise RuntimeError(
                "Flipping edge 4 should give weight 0, not {}.".format(f) )
    f = hb.flipWeight( tri.edge(5) )
    if f != 3:
        raise RuntimeError(
                "Flipping edge 5 should give weight 3, not {}.".format(f) )
    print( "    Passed." )
    print()
    print( "Reducibility" )
    if not hb.isReducible( tri.edge(4) ):
        raise RuntimeError( "Edge 4 incorrectly identified as irreducible." )
    if hb.isReducible( tri.edge(5) ):
        raise RuntimeError( "Edge 5 incorrectly identified as reducible." )
    redEdgeInds = set()
    for e in hb.reducibleEdges():
        redEdgeInds.add( e.index() )
    if redEdgeInds != {3,4}:
        raise RuntimeError(
                "Incorrectly yielded following reducible edges: {}.".format(
                    redEdgeInds ) )
    print( "    Passed." )
    print()
    print( "Flip edge, and resolvability" )
    emb = tri.edge(7).embedding(0)
    hb.flipEdge( tri.edge(4) )
    if not hb.isResolvable( emb.tetrahedron().edge( emb.edge() ) ):
        raise RuntimeError( "Failed to certify resolvable edge." )
    tri = Triangulation3.fromIsoSig(sigGenus2)
    w = (0,2,0,0,0,0,0,0,0)
    hb.setBouquet( tri, w )
    try:
        hb.isResolvable( tri.edge(1) )
    except ValueError as e:
        print( "    Correctly raised ValueError: {}".format(e) )
    else:
        raise RuntimeError( "Failed to detect isotopic petals." )
    print( "    Passed." )
    print()
    print( "PASSED." )
    print()

    # Cases that should raise ValueError while resolving.
    print( "============================================================" )
    print( " Genus 2: ValueError while resolving" )
    print( "------------------------------------------------------------" )
    print()
    print( "No reducible edges" )
    tri = Triangulation3.fromIsoSig(sigGenus2)
    w = (3,3,2,3,0,2,0,0,1)
    hb.setBouquet( tri, w, {4} )
    flipEmbs = [ tri.edge(i).embedding(0) for i in [0,1,5] ]
    for emb in flipEmbs:
        hb.flipEdge( emb.tetrahedron().edge( emb.edge() ) )
    try:
        hb.resolveAllPetals()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError(
                "Failed to detect non-existence of reducible edges." )
    print()
    print( "Normal curve after resolving" )
    tri = Triangulation3.fromIsoSig(sigGenus2)
    w = (3,2,3,0,3,0,1,1,1)
    hb.setBouquet( tri, w )
    try:
        hb.resolveAllPetals()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to detect normal curve." )
    print()
    print( "Transverse petals after resolving" )
    tri = Triangulation3.fromIsoSig(sigGenus2)
    w = (1,1,2,0,1,1,0,0,0)
    hb.setBouquet( tri, w )
    try:
        hb.resolveAllPetals()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to detect transverse petals." )
    print()
    print( "Disconnected complement after resolving" )
    tri = Triangulation3.fromIsoSig(sigGenus2)
    w = (1,0,1,0,1,0,0,0,0)
    hb.setBouquet( tri, w, {1} )
    try:
        hb.resolveAllPetals()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to detect disconnected complement." )
    print()
    print( "PASSED." )
    print()

    # Cases that should succeed.
    print( "============================================================" )
    print( " Genus 2: Cases that should succeed" )
    print( "------------------------------------------------------------" )
    print()
    success = [
            (0,0,0,2,3,1,2,1,0),
            (0,2,2,4,5,1,2,3,0) ]
    for w in success:
        print( "Weights: {}".format(w) )
        print( "------------------------------------" )
        print()

        # Test resolveAllPetals().
        print( "    HeegaardBuilder.resolveAllPetals()" )
        tri = Triangulation3.fromIsoSig(sigGenus2)
        hb.setBouquet( tri, w )
        if hb.fillHandlebody():
            raise RuntimeError( "Incorrectly filled with handlebody." )
        hb.resolveAllPetals()
        if not hb.fillHandlebody():
            raise RuntimeError( "Failed to fill with handlebody." )
        if hb.fillHandlebody():
            raise RuntimeError( "Incorrectly filled closed triangulation." )
        if not tri.isValid():
            raise RuntimeError( "Invalid after filling." )
        if not tri.isClosed():
            raise RuntimeError( "Not closed after filling." )
        if not tri.isOrientable():
            raise RuntimeError( "Non-orientable after filling." )
        print( "        {}: {}; {}.".format(
            "Original", tri.size(), tri.isoSig() ) )
        origSig = tri.isoSig()
        tri.intelligentSimplify()
        tri.intelligentSimplify()
        print( "        {}: {}; {}.".format(
            "...Final", tri.size(), tri.isoSig() ) )
        print()

        # Test resolveGreedily().
        print( "    HeegaardBuilder.resolveGreedily()" )
        tri = Triangulation3.fromIsoSig(sigGenus2)
        hb.setBouquet( tri, w )
        if hb.fillHandlebody():
            raise RuntimeError( "Incorrectly filled with handlebody." )
        hb.resolveGreedily()
        if not hb.fillHandlebody():
            raise RuntimeError( "Failed to fill with handlebody." )
        if hb.fillHandlebody():
            raise RuntimeError( "Incorrectly filled closed triangulation." )
        if not tri.isValid():
            raise RuntimeError( "Invalid after filling." )
        if not tri.isClosed():
            raise RuntimeError( "Not closed after filling." )
        if not tri.isOrientable():
            raise RuntimeError( "Non-orientable after filling." )
        print( "        {}: {}; {}.".format(
            "Original", tri.size(), tri.isoSig() ) )
        origSig = tri.isoSig()
        tri.intelligentSimplify()
        tri.intelligentSimplify()
        print( "        {}: {}; {}.".format(
            "...Final", tri.size(), tri.isoSig() ) )
        print()

        # Test resolveInAllWays().
        print( "    HeegaardBuilder.resolveInAllWays()" )
        tri = Triangulation3.fromIsoSig(sigGenus2)
        hb.setBouquet( tri, w )
        count = 0
        minSize = None
        minSig = None
        maxSize = 0
        maxSig = None
        sigs = set()
        for r in hb.resolveInAllWays():
            count += 1
            r.fillHandlebody()
            rTri = r.triangulation()
            sigs.add( rTri.isoSig() )
            if minSize is None or rTri.size() < minSize:
                minSize = rTri.size()
                minSig = rTri.isoSig()
            if rTri.size() > maxSize:
                maxSize = rTri.size()
                maxSig = rTri.isoSig()
        print( "        Total: {}. Unique: {}.".format( count, len(sigs) ) )
        print( "        MinSize: {}. MinSig: {}.".format(
            minSize, minSig ) )
        print( "        MaxSize: {}. MaxSig: {}.".format(
            maxSize, maxSig ) )
        print()
    print( "PASSED." )
    print()

    # Genus 3 cases that should raise ValueError.
    print( "============================================================" )
    print( " Genus 3: Cases that should raise ValueError" )
    print( "------------------------------------------------------------" )
    print()
    print( "Isotopic petals while resolving" )
    tri = Triangulation3.fromIsoSig(sigGenus3)
    w = (3,1,5,4,4,0,0,2,3,1,1,1,0,1,0)
    hb.setBouquet( tri, w )
    try:
        hb.resolveAllPetals()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to detect isotopic petals." )
    print()
    print( "Tranvserse petals after resolving" )
    tri = Triangulation3.fromIsoSig(sigGenus3)
    w = (2,1,2,3,3,1,1,2,3,0,4,0,0,0,0)
    hb.setBouquet( tri, w )
    try:
        hb.resolveAllPetals()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to detect transverse petals." )
    print()
    print( "Disconnected complement after resolving" )
    tri = Triangulation3.fromIsoSig(sigGenus3)
    w = (2,0,1,2,2,0,0,1,1,0,1,0,0,0,0)
    hb.setBouquet( tri, w )
    try:
        hb.resolveAllPetals()
    except ValueError as e:
        print( "    {}".format(e) )
    else:
        raise RuntimeError( "Failed to detect disconnected complement." )
    print()
    print( "PASSED." )
    print()

    # Genus 3 cases that should succedd.
    print( "============================================================" )
    print( " Genus 3: Cases that should succeed" )
    print( "------------------------------------------------------------" )
    print()
    w = (5,1,2,1,3,1,1,5,1,0,2,0,1,1,0)
    print( "Weights: {}".format(w) )
    print( "------------------------------------------------------" )
    print()
    print( "    HeegaardBuilder.resolveAllPetals()" )
    tri = Triangulation3.fromIsoSig(sigGenus3)
    hb.setBouquet( tri, w )
    if hb.fillHandlebody():
        raise RuntimeError( "Incorrectly filled with handlebody." )
    hb.resolveAllPetals()
    if not hb.fillHandlebody():
        raise RuntimeError( "Failed to fill with handlebody." )
    if hb.fillHandlebody():
        raise RuntimeError( "Incorrectly filled closed triangulation." )
    if not tri.isValid():
        raise RuntimeError( "Invalid after filling." )
    if not tri.isClosed():
        raise RuntimeError( "Not closed after filling." )
    if not tri.isOrientable():
        raise RuntimeError( "Non-orientable after filling." )
    print( "        {}: {}; {}.".format(
        "Original", tri.size(), tri.isoSig() ) )
    origSig = tri.isoSig()
    tri.intelligentSimplify()
    tri.intelligentSimplify()
    print( "        {}: {}; {}.".format(
        "...Final", tri.size(), tri.isoSig() ) )
    print()
    print( "    HeegaardBuilder.resolveGreedily()" )
    tri = Triangulation3.fromIsoSig(sigGenus3)
    hb.setBouquet( tri, w )
    if hb.fillHandlebody():
        raise RuntimeError( "Incorrectly filled with handlebody." )
    hb.resolveGreedily()
    if not hb.fillHandlebody():
        raise RuntimeError( "Failed to fill with handlebody." )
    if hb.fillHandlebody():
        raise RuntimeError( "Incorrectly filled closed triangulation." )
    if not tri.isValid():
        raise RuntimeError( "Invalid after filling." )
    if not tri.isClosed():
        raise RuntimeError( "Not closed after filling." )
    if not tri.isOrientable():
        raise RuntimeError( "Non-orientable after filling." )
    print( "        {}: {}; {}.".format(
        "Original", tri.size(), tri.isoSig() ) )
    origSig = tri.isoSig()
    tri.intelligentSimplify()
    tri.intelligentSimplify()
    print( "        {}: {}; {}.".format(
        "...Final", tri.size(), tri.isoSig() ) )
    print()
    print( "    HeegaardBuilder.resolveInAllWays(), truncated after 50" )
    tri = Triangulation3.fromIsoSig(sigGenus3)
    hb.setBouquet( tri, w )
    count = 0
    minSize = None
    minSig = None
    maxSize = 0
    maxSig = None
    sigs = set()
    for r in hb.resolveInAllWays():
        if count > 49:
            break
        count += 1
        r.fillHandlebody()
        rTri = r.triangulation()
        sigs.add( rTri.isoSig() )
        if minSize is None or rTri.size() < minSize:
            minSize = rTri.size()
            minSig = rTri.isoSig()
        if rTri.size() > maxSize:
            maxSize = rTri.size()
            maxSig = rTri.isoSig()
    print( "        Total: {}. Unique: {}.".format( count, len(sigs) ) )
    print( "        MinSize: {}. MinSig: {}.".format(
        minSize, minSig ) )
    print( "        MaxSize: {}. MaxSig: {}.".format(
        maxSize, maxSig ) )
    print()
    print( "PASSED." )
    print()
