from scipy.special import gamma, gammaln

def factorialln(n):
    if n == 1:
        return 0
    elif n == 0:
        return 0
    elif n < 0:
        raise Exception("factorial is not defined for n < 0")
    else:
        return gammaln(n + 1)
