def IsLess(a, b):
    return a < b

def IsMore(a, b):
    return a > b

def IsEqual(a, b, eps):
    return abs(a - b) < eps

def IsLessEqual(a, b, eps):
    return a < b or IsEqual(a, b, eps)

def IsMoreEqual(a, b, eps):
    return a > b or IsEqual(a, b, eps)