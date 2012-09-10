def assertDoesNotContain(needle, haystack):
    assert (haystack.find(needle) == -1)
    
def assertContains(needle, haystack):
    assert (haystack.find(needle) > -1)