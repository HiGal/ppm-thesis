import numpy as np


class Alignment:

    def __init__(self, N, M, q, Z, sequence, header, spec_name, spec_id, uniprot_id):
        self.N = N
        self.M = M
        self.q = q
        self.Z = Z
        self.sequence = sequence
        self.header = header
        self.spec_name = spec_name
        self.spec_id = spec_id
        self.uniprot_id = uniprot_id


class HarminizaedAlignments:

    def __init__(self, X1: Alignment, X2: Alignment):
        self.X1 = X1
        self.X2 = X2


class FreqC:

    def __init__(self, Pij, Pi, specs, M, matching, q):
        self.Pij = Pij
        self.Pi = Pi
        self.specs = specs
        self.M = M
        self.matching = matching
        self.q = q

    def FreqC(self, X1: Alignment, X2: Alignment):
        N = X1.N + X2.N
        s = X1.q - 1
        M = X1.M
        return np.zeros((s * N, s * N)), np.zeros(N * s), [], np.zeros((1, 2)), np.zeros(M), X1.q


class FastC:

    def __init__(self, Cij, specs, M):
        self.Cij = Cij
        self.specs = specs
        self.M = M

    def FastC(self, N):
        return np.zeros((N, N)), [], np.zeros((1, 2))
