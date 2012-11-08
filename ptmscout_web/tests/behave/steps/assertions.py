import re
def assertDoesNotContain(needle, haystack):
    assert (haystack.find(needle) == -1)
    
def assertContains(needle, haystack):
    if not (haystack.find(needle) > -1):
        print "'%s' not found in '%s'" % (needle, haystack)
        assert False
    
def assertEqual(expected, actual):
    if actual != expected:
        raise AssertionError("%s != %s" % (str(expected), str(actual)))


def assertRegexMatch(regex, haystack):
    m = re.search(regex, haystack)
    if m == None:
        assert False, "'%s' does not match '%s'" % (regex, haystack)