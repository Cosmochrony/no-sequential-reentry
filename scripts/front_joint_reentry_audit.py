#!/usr/bin/env python3
"""Joint exact audit of three proposed history-fibre re-entry levers.

Archived reproducibility script for the no-sequential-re-entry note.

The two required gates are:
  R. a natural linear map from the erased D4 history fibre to End(V_c);
  O. a natural orientation distinguishing (Z,w) from (-Z,w^-1).

The script certifies the elementary algebra used in Section AH.  Corpus typing
claims (O33 and EBJ) are audited from their deposited statements, not inferred
from this finite model.
"""

from fractions import Fraction


def mat_vec(matrix, vector):
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]


def heis_mul(g, h, q):
    a, b, z = g
    ap, bp, zp = h
    return ((a + ap) % q, (b + bp) % q, (z + zp + a * bp) % q)


def sigma_x(g, q):
    a, b, z = g
    return ((-a) % q, b, (-z) % q)


def audit_conditional_expectation():
    half = Fraction(1, 2)
    projection = [[half, half], [half, half]]
    alternating = [Fraction(1), Fraction(-1)]
    invariant = [Fraction(1), Fraction(1)]
    assert mat_vec(projection, alternating) == [0, 0]
    assert mat_vec(projection, invariant) == invariant

    # Abstract endpoint transport on the X/X^-1 fibre: the two endpoints
    # contribute with opposite odd symbols.  It is non-zero on ker(P).
    endpoint_odd = [[Fraction(1), Fraction(-1)]]
    assert mat_vec(endpoint_odd, alternating) == [Fraction(2)]

    # Representative downstream maps Q.  The theorem is general: Q(Pv)=0
    # for every linear Q whenever Pv=0.
    for downstream in (
        [[Fraction(1), Fraction(0)]],
        [[Fraction(2), Fraction(-3)], [Fraction(5), Fraction(7)]],
    ):
        assert mat_vec(downstream, mat_vec(projection, alternating)) == [0] * len(downstream)

    return {
        "alternating_fibre_vector": alternating,
        "conditional_expectation_kills_it": True,
        "endpoint_odd_map_is_nonzero_on_it": True,
        "factorisation_through_D4_impossible": True,
        "sequential_second_projection_cannot_restore_it": True,
    }


def audit_sigma_x(q=5):
    elements = [(a, b, z) for a in range(q) for b in range(q) for z in range(q)]
    for g in elements:
        for h in elements:
            left = sigma_x(heis_mul(g, h, q), q)
            right = heis_mul(sigma_x(g, q), sigma_x(h, q), q)
            assert left == right
    return {
        "q": q,
        "pairs_checked": len(elements) ** 2,
        "sigma_X_is_Heisenberg_automorphism": True,
        "abstract_cocycle_does_not_orient_the_axis": True,
    }


def main():
    ce = audit_conditional_expectation()
    sx = audit_sigma_x()
    print("== Joint re-entry audit: exact algebra ==")
    print(f"  Conditional expectation kills alternating fibre line: {ce['conditional_expectation_kills_it']}")
    print(f"  Endpoint odd map is non-zero on the same line: {ce['endpoint_odd_map_is_nonzero_on_it']}")
    print(f"  Factorisation through D4 is impossible: {ce['factorisation_through_D4_impossible']}")
    print(f"  Any sequential second projection still kills it: {ce['sequential_second_projection_cannot_restore_it']}")
    print(f"  sigma_X automorphism pairs checked at q={sx['q']}: {sx['pairs_checked']}")
    print("  sigma_X preserves the abstract Heisenberg product: True")
    print("== Corpus-typed verdict ==")
    print("  O33 supplies End(V_c) as an internal operator carrier, but no map from")
    print("  the D4 quotient; its Weyl carrier is phase-space based and parity graded.")
    print("  EBJ supplies a downstream Lorentzian chiral tangent only after P(s), J_Pi,")
    print("  and the odd modulus are given; it supplies no history-to-jet transport.")
    print("  Conditional expectation and every sequential projection tower erase, rather")
    print("  than re-enter, the alternating fibre character.")
    print("Verdict: all three deposited levers fail the joint R+O requirement.  A parallel")
    print("substrate channel with a central/Weyl orientation would be new architecture,")
    print("not a derivation from the current D4 channel.")


if __name__ == "__main__":
    main()
