"""
Small test suite for the SpineLayering class.
"""
from sys import stdout
from timeit import default_timer as timer
from gspine import *
from regina import *


if __name__ == "__main__":
    testData = [
            ( "eHuGabdes", 2, "Handlebody" ),
            ( "eHbecadjk", 2, "Handlebody" ),
            ( "nHuKfvPQPMabdgikhkkjlmhjfmdscnjex", 2,
                "Handlebody: eHbecadjk plus extra layers" ),
            ( "iLLLPQcbddegfhghabfsccswr", 2,
                "Attached handlebody to eHbecadjk" ),
            ( "lLLLLPPQcadegigiihjkkjaxreousjnck", 2,
                "Attached handlebody to eHbecadjk" ),
            ( "mLvLLMQQPaefhikighkjlljxfrtaangkjdj", 2,
                "Attached handlebody to eHbecadjk" ),
            ( "oLLvzwQMAQccdhikhghjlklmnnhshsaocnhvvnwlj", 2,
                "Attached handlebody to eHbecadjk" ),
            ( "hHbLbqiabegeti", 3, "Handlebody" ) ]
    msg = "Time: {:.6f} seconds. Spines: {}. Layerings: {}."
    for sig, genus, description in testData:
        tri = Triangulation3.fromIsoSig(sig)
        triangles = tri.countTriangles()
        print()
        print( "==========================================================" )
        print( "Iso sig: {}".format(sig) )
        print( "#triangles: {}".format(triangles) )
        print( "Genus {}".format(genus) )
        print(description)
        print( "----------------------------------------------------------" )
        print()
        stdout.flush()

        # Test searching for g-spines.
        spineCount = 0
        layeringCount = 0
        start = timer()
        for spineFaceInds in combinations( triangles, 2*genus-1 ):
            spineFaces = [ tri.triangle(i) for i in spineFaceInds ]
            layering = SpineLayering(spineFaces)
            subCount = 0
            for extension in layering.minimalCompleteExtensions():
                if subCount == 0:
                    detail = extension.detail()
                subCount += 1
                layeringCount += 1
            if subCount > 0:
                spineCount += 1
                print( ".................................................." )
                print( "  Spine faces: {}".format(spineFaceInds) )
                print( "   #layerings: {}".format(subCount) )
                print( " -------------" )
                print( " One layering:" )
                print()
                print(detail)
                stdout.flush()
        end = timer()
        print()
        print( ".................................................." )
        print()
        print( msg.format( end - start, spineCount, layeringCount ) )
        print()
        stdout.flush()
