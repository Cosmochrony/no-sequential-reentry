#!/usr/bin/env python3
"""Axis-selector front: do the deposited enriched channel data select
i(A_X + A_Y) against i(A_X - A_Y)?

Archived reproducibility script for the no-sequential-re-entry note.

Order fixed by Jerome (15 July 2026): (1) stabiliser of the enriched data
(oriented b-shadow, central cocycle, character c, endpoint representation);
(2) equivariant-selector test; (3) special check whether the central register
at lengths 2-3 fixes the relative sign; kill-switch: if the two axes remain
exchanged by a deposited symmetry, close the route; only if an axis is
selected, audit the D4 -> End(V_c) linear re-entry legitimacy.

Author-review calibration: preserving a chosen central coordinate pointwise
is stronger than preserving the Heisenberg cocycle as an abstract structure.
The map sigma_X:(a,b,z)->(-a,b,-z) is a group automorphism, implemented
semilinearly by RK.  Thus the cocycle alone does not remove the exchanger;
only an additional orientation/pinning of the central and Weyl generators
does so.

Objects: A_X = W(X) - W(-X), A_Y = W(Y) - W(-Y), T_pm = i(A_X +- A_Y).
Candidate deposited symmetries acting on transported operators:
  frame rotation  w  : (a,b) -> phase-space rotation  (destroys the shadow
                       labelling: X-letters become b-emitters);
  fibre swap  sigma_X: X <-> X^{-1}, i.e. (a,b,z) -> (-a, b, -z);
  fibre swap  sigma_Y: Y <-> Y^{-1}, i.e. (a,b,z) -> (a, -b, -z);
  parity      R      : (a,b,z) -> (-a,-b, z);
  conjugation K      : antiunitary, W_c(a,b) -> W_c(a,-b);
  RK                 : antiunitary, W_c(a,b) -> W_c(-a,b);
  central flip       : c -> -c (changes the deposited block).
Conventions as in front_hist_fibre_doubling.py (O33 / TrajectoryBranching).
"""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

TOL = 1e-11


def weyl_op(q, c, a, b):
    inv2 = pow(2, -1, q)
    mat = np.zeros((q, q), dtype=complex)
    for x in range(q):
        phase = (b * ((x - a * inv2) % q)) % q
        mat[x, (x - a) % q] = np.exp(2j * np.pi * c * phase / q)
    return mat


def rho_op(q, c, a, b, z):
    mat = np.zeros((q, q), dtype=complex)
    for x in range(q):
        mat[x, (x - a) % q] = np.exp(2j * np.pi * c * ((z + b * (x - a)) % q) / q)
    return mat


def parity_op(q):
    mat = np.zeros((q, q))
    for x in range(q):
        mat[x, (-x) % q] = 1.0
    return mat


def fourier_op(q, c):
    return np.array([[np.exp(2j * np.pi * c * (x * y % q) / q) for y in range(q)]
                     for x in range(q)]) / math.sqrt(q)


def weyl_coeff(q, mat, a, b):
    """Coefficient of W(a,b) in mat (Weyl orthogonality: tr(W^dag W') = q)."""
    return np.trace(weyl_op(q, 1, a, b).conj().T @ mat) / q


def axis_action(name, image_TP, TP, TM):
    for label, ref in (("+T+", TP), ("-T+", -TP), ("+T-", TM), ("-T-", -TM)):
        if np.abs(image_TP - ref).max() < TOL:
            return label
    return "outside plane"


def part1_symmetry_table(q=29, c=1):
    print("== 1. Enriched-data stabiliser and axis action ==")
    AX = weyl_op(q, c, 1, 0) - weyl_op(q, c, -1, 0)
    AY = weyl_op(q, c, 0, 1) - weyl_op(q, c, 0, -1)
    TP, TM = 1j * (AX + AY), 1j * (AX - AY)
    R = parity_op(q)
    F = fourier_op(q, c)
    Finv = np.linalg.inv(F)

    def ad_unitary(U, A):
        return U @ A @ np.linalg.inv(U)

    def ad_antiunitary_K(A):        # K A K = conj(A) in the position basis
        return A.conj()

    def ad_RK(A):
        return R @ A.conj() @ R

    rows = []
    rows.append(("w (frame, unitary F)", "no (X becomes a b-emitter)", "n/a", "yes",
                 axis_action("w", ad_unitary(F, TP), TP, TM)))
    sx = 1j * ((-1) * AX + AY)      # sigma_X: A_X -> -A_X, A_Y -> A_Y
    rows.append(("sigma_X (fibre swap)", "yes (oriented)", "NO (z -> -z)", "yes",
                 axis_action("sx", sx, TP, TM)))
    sy = 1j * (AX - AY)
    rows.append(("sigma_Y (fibre swap)", "NO (orientation)", "NO (z -> -z)", "yes",
                 axis_action("sy", sy, TP, TM)))
    rows.append(("R (parity)", "NO (orientation)", "yes", "yes",
                 axis_action("R", ad_unitary(R, TP), TP, TM)))
    rows.append(("K (antiunitary)", "NO (orientation)", "NO", "yes",
                 axis_action("K", 1j * (ad_antiunitary_K(AX) - -ad_antiunitary_K(AY)) * 0
                             + (-1j) * (ad_antiunitary_K(AX) + ad_antiunitary_K(AY)),
                             TP, TM)))
    rows.append(("RK (antiunitary)", "yes (oriented)", "NO (z -> -z)", "yes",
                 axis_action("RK", (-1j) * (ad_RK(AX) + ad_RK(AY)), TP, TM)))
    rows.append(("c -> -c (pair exchange)", "yes", "yes", "fixed block: NO",
                 "exchanges axes / conjugate factors"))

    print(f"  {'symmetry':26s} {'oriented shadow':17s} {'fixes z coord':11s} {'char c':7s} action on T+")
    for name, shd, coc, cc, act in rows:
        print(f"  {name:26s} {shd:17s} {coc:11s} {cc:7s} {act}")

    swap_check = np.abs(ad_RK(rho_op(q, c, 3, 2, 5)) - rho_op(q, c, -3 % q, 2, -5 % q)).max()
    print(f"  Identity check: Ad_RK rho_c(a,b,z) = rho_c(-a,b,-z): error {swap_check:.2e}")
    print("  Hence the fibre swap sigma_X is implemented on transports EXACTLY by the")
    print("  deposited antiunitary RK.  It reverses the chosen z coordinate, but it")
    print("  PRESERVES the Heisenberg structure: sigma_X(gg')=sigma_X(g)sigma_X(g')")
    print("  because -(z+z'+ab')=(-z)+(-z')+(-a)b'.")
    print("  Therefore RK remains an exchanger of the abstract cocycle-enriched datum.")
    print("  It is excluded only after additionally requiring the named central")
    print("  coordinate (and the named Weyl generator) to be fixed pointwise.")
    print("  The kill-switch fires for strict D4 AND for abstract cocycle enrichment;")
    print("  it becomes silent only in the stronger coordinate-pinned model.")
    return AX, AY, TP, TM, F, Finv, R


def part2_selector(q, c, AX, AY, TP, TM, F, Finv, R):
    print("== 2. Coordinate-pinned selector test ==")
    tsel = 1j * (AX + F @ AX @ Finv)
    print(f"  T_sel = i(A_X + F A_X F^-1) vs T_plus: error {np.abs(tsel - TP).max():.2e}")
    tsel_flip = 1j * (-AX + F @ (-AX) @ Finv)
    print(f"  Global fibre relabelling (A_X -> -A_X): T_sel -> -T_sel "
          f"(error {np.abs(tsel_flip + tsel).max():.2e}): the AXIS is Z2-robust.")
    tsel_inv = 1j * (AX + Finv @ AX @ F)
    print(f"  Choosing F^-1 instead: gives T_minus (error {np.abs(tsel_inv - TM).max():.2e}),")
    print("  In the chosen O12/Korinman coordinates, the concrete kernel")
    print("  F[x,y] = q^(-1/2) psi_c(xy) selects the first axis (Section P match).")
    Fc = fourier_op(q, (q - c) % q)
    tsel_c = 1j * (AX + Fc @ AX @ np.linalg.inv(Fc))
    print(f"  c -> -c: the same construction gives T_minus (error "
          f"{np.abs(tsel_c - TM).max():.2e}): the selector is c-ODD, and c is pinned")
    print("  per deposited block (c = c1+c2+c3 != 0).")
    print("  Calibration: the Trajectory Branching D4 observable includes the conjugate")
    print("  factor q-c.  At that paired level c -> -c exchanges the two factors and")
    print("  further prevents the fixed-c coordinate choice from being a pairwise axis.")
    rkF = R @ F.conj() @ R
    for label, ref in (("F", F), ("F^-1", Finv), ("F^dag", F.conj().T)):
        err = np.abs(rkF - ref).max()
        print(f"  RK F RK vs {label}: {err:.2e}")
    print("  RK sends the named operator F to F^-1.  This excludes RK from the")
    print("  POINTWISE stabiliser of the coordinate presentation, but is also exactly")
    print("  the covariance w -> w^-1.  It therefore does not establish a canonical")
    print("  axis under isomorphisms of the abstract Heisenberg-Weil datum.")


def words(n):
    INV = {0: 1, 1: 0, 2: 3, 3: 2}
    if n == 0:
        yield ()
        return
    for first in range(4):
        stack = [(first,)]
        while stack:
            w = stack.pop()
            if len(w) == n:
                yield w
            else:
                for nxt in range(4):
                    if nxt != INV[w[-1]]:
                        stack.append(w + (nxt,))


B_EMIT = {0: 0, 1: 0, 2: 1, 3: -1}
A_EMIT = {0: 1, 1: -1, 2: 0, 3: 0}


def endpoint(w):
    a = b = z = 0
    for g in w:
        z += a * B_EMIT[g]
        a += A_EMIT[g]
        b += B_EMIT[g]
    return a, b, z


def alt_sign(w):
    sign = 1
    prev = None
    for g in w:
        s = B_EMIT[g]
        if s == 0 and prev != 0 and g == 1:
            sign = -sign
        prev = s
    return sign


def part3_lengths23(q=29, c=1):
    print("== 3. Central register at lengths 2-3 ==")
    R = parity_op(q)
    for n in (2, 3):
        fibres = defaultdict(list)
        for w in words(n):
            fibres[tuple(B_EMIT[g] for g in w)].append(w)
        planeX = planeY = both = 0
        planeX_nz = planeY_nz = 0
        for shd, fib in fibres.items():
            if all(s != 0 for s in shd):
                continue
            for use_z in (True, False):
                if not use_z and n == 2:
                    continue
                op = np.zeros((q, q), dtype=complex)
                for w in fib:
                    a, b, z = endpoint(w)
                    zz = z % q if use_z else 0
                    op += alt_sign(w) * rho_op(q, c, a % q, b % q, zz)
                odd = 0.5 * (op - R @ op @ R)
                cx = abs(weyl_coeff(q, odd, 1, 0)) + abs(weyl_coeff(q, odd, q - 1, 0))
                cy = abs(weyl_coeff(q, odd, 0, 1)) + abs(weyl_coeff(q, odd, 0, q - 1))
                if use_z:
                    if cx > TOL and cy > TOL:
                        both += 1
                    elif cx > TOL:
                        planeX += 1
                        planeX_nz += 1
                    elif cy > TOL:
                        planeY += 1
                        planeY_nz += 1
        print(f"  n={n}: fibres with BOTH A_X and A_Y components: {both}; "
              f"A_X-line only: {planeX}; A_Y-line only: {planeY}")
    print("  A single fibre transport is supported on one line b = b(shadow), so no")
    print("  fibre can carry A_X (b=0) and A_Y (b=+-1) components simultaneously:")
    print("  the relative sign CANNOT come from single-fibre data, with or without z.")
    print("  The central register at lengths 2-3 therefore does NOT by itself fix the")
    print("  relative sign.  Retaining its values excludes RK only if an orientation")
    print("  of the central coordinate is separately declared fixed.")


def main():
    q, c = 29, 1
    AX, AY, TP, TM, F, Finv, R = part1_symmetry_table(q, c)
    part2_selector(q, c, AX, AY, TP, TM, F, Finv, R)
    part3_lengths23(q, c)
    print("Verdict: RK exchanges the two axes in strict D4 and remains a semilinear")
    print("automorphism after abstract cocycle enrichment.  The kill-switch therefore")
    print("still FIRES at both levels.  A coordinate-pinned O12/Korinman presentation")
    print("selects T_plus and is Z2-robust and c-odd, but that is conditional on an")
    print("extra orientation of the central/Weyl generators, not on the cocycle alone.")
    print("Re-entry legitimacy and axis naturality remain two gates unless one future")
    print("projection principle supplies both simultaneously.")


if __name__ == "__main__":
    main()
