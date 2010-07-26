import math

def logsumexp(lnv):
    """Sum exp(item) for item in lnv (log-normal vector) without overflow."""
    n = lnv[0]
    maxAbs = n
    minN = n
    maxN = n
    c = n
    for item in lnv[1:]:
        n = item
        if n > maxN:
            maxN = n
        if abs(n) > maxAbs:
            maxAbs = abs(n)
        if n < minN:
            minN = n
    if maxAbs > maxN:
        c = minN
    else:
        c = maxN
    return c + math.log(sum([math.exp(i - c) for i in lnv]))

