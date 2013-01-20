import re

def assertDoesNotContain(needle, haystack):
    if (haystack.find(needle) > -1):
        assert False, "'%s' was found in '%s'" % (needle, haystack)
    
def assertContains(needle, haystack):
    if not (haystack.find(needle) > -1):
        assert False, "'%s' not found in '%s'" % (needle, haystack)
   

def assertEqual(expected, actual):
    if actual != expected:
        raise AssertionError("\n%s\n    !=\n%s" % (str(expected), str(actual)))

def assertRegexMatch(regex, haystack):
    m = re.search(regex, haystack)
    if m == None:
        assert False, "'%s' does not match '%s'" % (regex, haystack)

def assertIn(needle, haystack):
    if needle not in haystack:
        assert False, "'%s' not found in '%s'" % (str(needle), str(haystack))
