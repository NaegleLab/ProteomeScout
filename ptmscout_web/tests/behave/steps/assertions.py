import re
import time
import difflib

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

def assertAlmostEqual(expected, value, tolerance=0.001):
    test = (expected - tolerance <= value and value <= expected + tolerance)
    if not test:
        assert False, "%f != %f +- %f" % (value, expected, tolerance)

def synchronous_assert_called(mock, limit=1):
    slept_time = 0.0
    while not mock.called and slept_time < limit:
        time.sleep(0.1)
        slept_time += 0.1    
    
    assert mock.called, "Synchronous call to %s timed out" % (str(mock))

def diff(s1, s2):
    s = difflib.SequenceMatcher(a=s1, b=s2)
    output = "\n"
    for tag, i1, i2, j1, j2 in s.get_opcodes():
        output += "%-10s: %s  <-->  %s\n" % (tag, s1[i1:i2], s2[j1:j2])
    return output

def assertTextFilesEqual(fn1, fn2):
    i = 0
    with open(fn1,'r') as f1:
        with open(fn2,'r') as f2:
            for f1l in f1:
                i += 1
                f2l = f2.readline()
                assert f1l == f2l, "Files differ at line: %d\n%s\n" % (i, diff(f1l, f2l)) 