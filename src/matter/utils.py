import importlib
def almostEqual(*args):
    tol = 10^-7
    for val1 in args:
        for val2 in args:
            if abs(val1-val2)>tol: return False # break if it is ever false
    #otherwise return true
    return True

def isfloat(s):
    """True if argument can be converted to float"""
    try:
        x = float(s)
        return True
    except ValueError:
        pass
    return False


def getModules( package ):
    "get and import a list of python modules under given python package"
    packagename = package.__name__
    path = package.__path__[0]
    import os
    rt = []
    for entry in os.listdir( path ):
        f = os.path.join( path, entry )
        if os.path.isdir( f ): continue
        if entry.find( "__init__.py" ) != -1: continue
        modulename, ext = os.path.splitext( entry )
        if ext not in [".py", '.pyc'] : continue
        # exec("from %s import %s as m" % (packagename, modulename))
        m = importlib.import_module("%s.%s"%(packagename,modulename))
        rt.append( m )
        continue
    return rt



def test_getModules():
    from . import atomic_properties
    print((getModules( atomic_properties )))
    return 


def main():
    test_getModules()
    return 

if __name__ == "__main__": main()
