from .logging import logger
import math


class RnsBasis:
    def __init__(self, moduli):
        RnsBasis.ensure_coprime_moduli(moduli)
        logger.debug(f"New RnsBasis: {moduli}")
        self.moduli = moduli
        self.max = math.prod(moduli)

    @staticmethod
    def ensure_coprime_moduli(moduli):
        """
        Check if the moduli are pairwise coprime

        Args:
        moduli (iterable): Array of moduli for the RNS representation

        Raises:
        ValueError: If the moduli are not pairwise coprime
        """
        for i, m1 in enumerate(moduli):
            for m2 in moduli[i + 1:]:
                if math.gcd(m1, m2) != 1:
                    raise ValueError("Moduli must be pairwise coprime")

    def __repr__(self):
        return f"RnsBasis(moduli={self.moduli})"

    def __len__(self):
        return len(self.moduli)

    def __eq__(self, other):
        return self.moduli == other.moduli

    def to_rns(self, value):
        logger.debug(f"{value} -> {self}")
        residues = [value % m for m in self.moduli]
        return Rns(self, residues)


class Rns:
    def __init__(self, basis, residues):
        self._basis = basis
        self.residues = residues
        logger.debug(f"Created: {self}")

    @property
    def basis(self):
        return self._basis

    @staticmethod
    def from_residues(basis, residues):
        if len(basis) != len(residues):
            raise ValueError("Number of residues must match number of moduli")
        logger.debug(f"Creating RNS decomposition from residues: {residues}")
        return Rns(basis, residues)

    def __repr__(self):
        comps = ', '.join(f'{r} mod {m}' for (r, m) in zip(
            self.residues, self.basis.moduli))
        return f"Rns({comps})"

    def coeffwise_op(self, other, op):
        if not (isinstance(other, Rns) and self.basis == other.basis):
            raise ValueError(f"Incompatible basis: {self.basis}"
                             f" != {other.basis}")
        logger.debug(f"Performing operation {op} on {self} and {other}")
        residues = [op(a, b) % m for a, b, m in zip(
                    self.residues, other.residues, self.basis.moduli)]
        return Rns(self.basis, residues)

    def __getitem__(self, key):
        return self.residues[key]

    def __add__(self, other):
        return self.coeffwise_op(other, lambda a, b: a + b)

    def __sub__(self, other):
        return self.coeffwise_op(other, lambda a, b: a - b)

    def __mul__(self, other):
        return self.coeffwise_op(other, lambda a, b: a * b)

    def to_int(self):
        """
        Convert the RNS representation to an integer

        Returns:
        int: Integer representation of the RNS representation
        """

        M = math.prod(self.basis.moduli)
        total = 0
        for mi, ai in zip(self.basis.moduli, self.residues):
            Mi = M // mi
            inv = pow(int(Mi), -1, int(mi))
            total += ai * Mi * inv
        result = total % M
        logger.debug(f"{self} -> {result}")
        return result
