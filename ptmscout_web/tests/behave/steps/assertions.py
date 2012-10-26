def assertDoesNotContain(needle, haystack):
    assert (haystack.find(needle) == -1)
    
def assertContains(needle, haystack):
    assert (haystack.find(needle) > -1)
    
def assertEqual(expected, actual):
    if actual != expected:
        raise AssertionError("%s != %s" % (str(expected), str(actual)))
