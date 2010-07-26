import math

M_LN10 = math.log(10)
M_LOG10E = math.log10(math.e)

def phred2ln(qual):
    return M_LN10 * qual * -.1

def ln2phred(prob):
    return -10 * M_LOG10E * prob

def phred2float(qual):
    return math.pow(10, qual * -.1)

def float2phred(prob):
    return min(-10 * math.log10(prob), 99)

