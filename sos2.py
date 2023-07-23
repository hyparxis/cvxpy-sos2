import cvxpy as cp

def SOS2(n):
    # Gray code table for n bits
    def graycode(n):
        code = np.ndarray((1 << n, n))
        for i in range(1 << n):
            g = format(i ^ (i >> 1), f'0{n}b')
            for j in range(n):
                code[i, j] = int(g[j])
        return code
    
    # Log2 rounded up to nearest int
    def clog2(x):
        return np.ceil(np.log2(x)).astype(int)

    # Knot point k is in segment s
    def segment_contains_knotpoint(s, k):
        return s == k or k == s + 1

    λ = cp.Variable(n, pos=True)
    z = cp.Variable(clog2(n), boolean=True)
    g = graycode(clog2(n))

    constraints = []
    # Loop over binary variables and their possible states
    for i in range(clog2(n)):
        for b in (0, 1):
            k_not_in_gib = {k for k in range(n)}
            # Loop over segments
            for s in range(n - 1):
                # Loop over knot points
                for k in range(n):
                    # Segment contains knot point and g_i == b is in segment's 
                    # gray code parameterization
                    if segment_contains_knotpoint(s, k) and g[s, i] == b:
                        k_not_in_gib -= {k}
            if b == 1:
                constraints += [sum([λ[k] for k in k_not_in_gib]) <= 1 - z[i]]
            else:
                constraints += [sum([λ[k] for k in k_not_in_gib]) <= z[i]]
    constraints += [cp.sum(λ) == 1]
    return λ, constraints
