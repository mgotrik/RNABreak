from basicfunc import getPairsFromList

def test_getPairsFromList():
    assert(getPairsFromList([1,2,3,4]) == [(1,2), (3,4)])
    
    assert(getPairsFromList([1,2,3,4,5,6]) == [(1,2), (3,4), (5,6)])

