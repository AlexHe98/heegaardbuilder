from heegaardbuilder import *
from gspine import *
from regina import *


def _testLayering( tri, k ):
    for spineFaceInds in combinations( tri.countTriangles(), k ):
        spineFaces = [ tri.triangle(i) for i in spineFaceInds ]
        layering = SpineLayering(spineFaces)
        count = 0
        for extension in layering.minimalCompleteExtensions():
            if count == 0:
                print( spineFaceInds, extension._layeredTetInds )
            count += 1
        if count > 0:
            print(count)
            print()



if __name__ == "__main__":
    print()
    sigGenus2 = "eHbecadjk"

    # Handlebody.
    print( "===============" )
    tri = Triangulation3.fromIsoSig(sigGenus2)
    print( tri.isoSig() )
    _testLayering( tri, 3 )

    # After resolving.
    print( "===============" )
    w = (0,2,2,4,5,1,2,3,0)
    hb = HeegaardBuilder()
    hb.setBouquet( tri, w )
    hb.resolveGreedily()
    print( tri.isoSig() )
    _testLayering( tri, 3 )

    # After filling.
    print( "===============" )
    hb.fillHandlebody()
    print( tri.isoSig() )
    _testLayering( tri, 3 )

    # Directly from iso sig nHuKfvPQPMabdgikhkkjlmhjfmdscnjex.
    print( "===============" )
    tri = Triangulation3.fromIsoSig( "nHuKfvPQPMabdgikhkkjlmhjfmdscnjex" )
    print( tri.isoSig() )
    _testLayering( tri, 3 )

    # Directly from iso sig lLLLLPPQcadegigiihjkkjaxreousjnck.
    print( "===============" )
    tri = Triangulation3.fromIsoSig( "lLLLLPPQcadegigiihjkkjaxreousjnck" )
    print( tri.isoSig() )
    _testLayering( tri, 3 )

    # Another genus 2.
    print( "===============" )
    tri = Triangulation3.fromIsoSig( "eHuGabdes" )
    print( tri.isoSig() )
    _testLayering( tri, 3 )

    # Genus 3.
    print( "===============" )
    sigGenus3 = "hHbLbqiabegeti"
    tri = Triangulation3.fromIsoSig(sigGenus3)
    print( tri.isoSig() )
    _testLayering( tri, 5 )
