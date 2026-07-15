#!/usr/bin/env python3
"""Erased-history-fibre front: is the missing ambient doubling carried by the
X-sign fibres of the Trajectory Branching canonical channel?

Archived reproducibility script for the no-sequential-re-entry note.

Pre-registered plan (Jerome, 15 July 2026):
  A. frame-stabiliser lemma: Stab_{SL(2,Fq)}({+-e1,+-e2}) = <w> ~ C4;
  B. exact fibre decomposition of the non-backtracking history space over the
     b-shadow channel (one sign per maximal zero-run; central register);
  C. C4 module structure of the history space (freeness);
  D. induced multiplicity profiles (full words from C4; shadows from C2):
     abundance test (kill-switch 1);
  E. fibre gauge group: per-run sign flips preserve the channel; abelian
     (kill-switch 3, direct form);
  F. operator transport: the flip-antisymmetric fibre line lands in the
     parity-odd sector D_minus.  Test separately whether C4 frame completion
     selects the odd Harper operator or leaves a two-axis ambiguity.

Pre-registered kill-switches (session order):
  KS-frame: a single-alphabet history space is not a canonical SL(2,Fq)
            carrier (fires if the stabiliser is exactly C4 -- expected);
  KS-1: inductions produce only abundant multiplicities, no privileged M2;
  KS-2: the only doublet is the Pell transfer vector (abelian commutant);
  KS-3: the fibre carries canonically only (Z2)^r, no forced non-abelian
        extension;
  KS-5: any lift that assumes SU(2) before finding it is circular.

Conventions follow O33 (SpectralO33.tex): rho_c(a,b,z) f(x) = psi_c(z + b(x-a)) f(x-a),
W_c(a,b) f(x) = psi_c(b(x - a/2)) f(x-a), and the Weyl covariance
rho(g) W(v) rho(g)^{-1} = W(phi(g) v), phi(g) = R0 g^{-T} R0.
Trajectory conventions follow TrajectoryBranching.tex: alphabet
S = {X, X^{-1}, Y, Y^{-1}}, non-backtracking, b-emissions (0, 0, +1, -1),
Heisenberg product (a,b,z)(a',b',z') = (a+a', b+b', z+z'+ab').
"""

from __future__ import annotations

import itertools
import math
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------------------
# Part A: frame-stabiliser lemma
# ---------------------------------------------------------------------------

def sl2_elements(q: int):
    for a in range(q):
        for b in range(q):
            for c in range(q):
                if a == 0:
                    if (-b * c) % q == 1:
                        for d in range(q):
                            yield (a, b, c, d)
                else:
                    d = (1 + b * c) * pow(a, -1, q) % q
                    yield (a, b, c, d)


def frame_stabiliser(q: int):
    """Exhaustive stabiliser of the set {+-e1, +-e2} in SL(2, Fq)."""
    frame = {(1, 0), (q - 1, 0), (0, 1), (0, q - 1)}
    stab = []
    for a, b, c, d in sl2_elements(q):
        image = {((a * x + b * y) % q, (c * x + d * y) % q) for x, y in frame}
        if image == frame:
            stab.append((a, b, c, d))
    return stab


def part_A(qs=(13, 29, 53)):
    print("== A. Frame stabiliser ==")
    for q in qs:
        stab = frame_stabiliser(q)
        w = (0, q - 1, 1, 0)
        expected = set()
        m = (1, 0, 0, 1)

        def mul(u, v):
            return ((u[0] * v[0] + u[1] * v[2]) % q, (u[0] * v[1] + u[1] * v[3]) % q,
                    (u[2] * v[0] + u[3] * v[2]) % q, (u[2] * v[1] + u[3] * v[3]) % q)

        g = m
        for _ in range(4):
            expected.add(g)
            g = mul(g, w)
        assert set(stab) == expected, (q, stab)
        print(f"  q={q}: |Stab| = {len(stab)} = C4 = <w>, exhaustive over {q**3 - q} elements: OK")
    print("  Analytic proof: a stabilising matrix is monomial with entries +-1;")
    print("  det=1 leaves exactly {I, -I, w, w^3}.  KS-frame CONFIRMED: the fixed")
    print("  Cayley alphabet is stabilised only by C4, so a single-frame history")
    print("  space is not a canonical SL(2,Fq) carrier; induction is mandatory.")


# ---------------------------------------------------------------------------
# Part B: history space, shadows, fibres
# ---------------------------------------------------------------------------

# letters: 0 = X, 1 = X^{-1}, 2 = Y, 3 = Y^{-1}
INV = {0: 1, 1: 0, 2: 3, 3: 2}
B_EMIT = {0: 0, 1: 0, 2: 1, 3: -1}
A_EMIT = {0: 1, 1: -1, 2: 0, 3: 0}


def words(n: int):
    if n == 0:
        yield ()
        return
    for first in range(4):
        stack = [(first,)]
        while stack:
            wrd = stack.pop()
            if len(wrd) == n:
                yield wrd
            else:
                for nxt in range(4):
                    if nxt != INV[wrd[-1]]:
                        stack.append(wrd + (nxt,))


def shadow_of(wrd):
    return tuple(B_EMIT[g] for g in wrd)


def endpoint(wrd):
    """Heisenberg endpoint (a, b, z) over Z (exact integers)."""
    a = b = z = 0
    for g in wrd:
        da, db = A_EMIT[g], B_EMIT[g]
        z += a * db
        a += da
        b += db
    return a, b, z


def zero_runs(shadow):
    runs = 0
    prev = None
    for s in shadow:
        if s == 0 and prev != 0:
            runs += 1
        prev = s
    return runs


def part_B(n_max=10):
    print("== B. Fibre decomposition of the history space ==")
    pell = {1: 3, 2: 7}
    for n in range(3, n_max + 1):
        pell[n] = 2 * pell[n - 1] + pell[n - 2]
    summary = {}
    for n in range(1, n_max + 1):
        fibres = defaultdict(list)
        total = 0
        for wrd in words(n):
            total += 1
            fibres[shadow_of(wrd)].append(wrd)
        assert total == 4 * 3 ** (n - 1)
        assert len(fibres) == pell[n], (n, len(fibres), pell[n])
        sign_ok = True
        sum_check = 0
        max_a_coll = 0
        z_splits = 0
        az_splits = 0
        for shd, fib in fibres.items():
            r = zero_runs(shd)
            sum_check += 2 ** r
            if len(fib) != 2 ** r:
                sign_ok = False
            eps = defaultdict(int)
            for wrd in fib:
                a, b, z = endpoint(wrd)
                eps[(a, z)] += 1
            avals = {a for a, _ in eps}
            max_a_coll = max(max_a_coll, len(fib) - len(avals))
            if len(eps) < len(fib):
                az_splits += 1
            zs = defaultdict(set)
            for (a, z), _ in eps.items():
                zs[a].add(z)
            if any(len(v) > 1 for v in zs.values()):
                z_splits += 1
        assert sum_check == total
        summary[n] = (pell[n], sign_ok, max_a_coll, az_splits)
        print(f"  n={n}: shadows={pell[n]} (Pell OK), words={total}, "
              f"fibre=2^runs everywhere: {sign_ok}, "
              f"fibres with (a,z)-collisions: {az_splits}, "
              f"max within-fibre a-collision count: {max_a_coll}")
    print("  Exact structure: word <-> (b-shadow, one sign per maximal zero-run);")
    print("  the central register z and displacement a are functions of the pair.")
    return summary


# ---------------------------------------------------------------------------
# Part C/D: C4 structure and induced multiplicity profiles
# ---------------------------------------------------------------------------

ROT = {0: 2, 2: 1, 1: 3, 3: 0}   # S: X -> Y -> X^{-1} -> Y^{-1} -> X


def part_CD(n_max=8):
    print("== C. C4 structure of the history space ==")
    for n in range(1, n_max + 1):
        fixed = {1: 0, 2: 0, 3: 0}
        for wrd in words(n):
            img = wrd
            for k in (1, 2, 3):
                img = tuple(ROT[g] for g in img)
                if img == wrd:
                    fixed[k] += 1
        assert fixed == {1: 0, 2: 0, 3: 0}, (n, fixed)
    print(f"  S^k has zero fixed words for k=1,2,3, all n <= {n_max}: the C4 action is FREE.")
    print("  Hence C[Words_n] = 3^(n-1) . C[C4] and")
    print("  Ind_C4^G C[Words_n] = 3^(n-1) . C[G]  (regular representation).")
    print("  End_G = (+)_pi M_(3^(n-1) dim pi): maximally abundant, no privileged factor.")
    print("  KS-1 FIRES for the raw frame-completed history space.")
    print("== D. Shadow induction from C2 = {+-I} ==")
    for n in range(1, n_max + 1):
        shadows = {shadow_of(wrd) for wrd in words(n)}
        neg_fixed = sum(1 for s in shadows if tuple(-x for x in s) == s)
        assert neg_fixed == 1, (n, neg_fixed)
    print("  -I acts on shadows by global negation with exactly ONE fixed shadow")
    print("  (the all-zero shadow), for every n tested.")
    print("  Hence m_pi(Ind_C2^G C[Shadows_n]) = (dim pi)(N_b(n) + eps_pi)/2 for every")
    print("  irreducible pi with central sign eps_pi: every irrep appears with")
    print("  multiplicity ~ N_b(n)/2 . dim pi.  KS-1 FIRES for the shadow induction too.")
    print("  KS-2 (analytic): the Pell transfer matrix [[2,1],[1,0]] has distinct")
    print("  eigenvalues 1+-sqrt2, so its commutant is C(+)C, abelian: the counting")
    print("  vector carries no doublet.  CONFIRMED without computation.")


# ---------------------------------------------------------------------------
# Part E: fibre gauge group
# ---------------------------------------------------------------------------

def flip_run(wrd, run_index):
    """Flip the sign of the run_index-th maximal zero-run (X <-> X^{-1})."""
    out = list(wrd)
    run = -1
    prev = None
    for i, g in enumerate(wrd):
        s = B_EMIT[g]
        if s == 0 and prev != 0:
            run += 1
        if s == 0 and run == run_index:
            out[i] = INV[g]
        prev = s
    return tuple(out)


def part_E(n_max=8):
    print("== E. Fibre gauge group ==")
    for n in range(1, n_max + 1):
        for wrd in words(n):
            shd = shadow_of(wrd)
            r = zero_runs(shd)
            for k in range(r):
                img = flip_run(wrd, k)
                ok = img != wrd
                ok = ok and shadow_of(img) == shd
                ok = ok and all(img[i + 1] != INV[img[i]] for i in range(n - 1))
                assert ok, (wrd, k, img)
    print(f"  Per-run sign flips are channel-preserving involutions for all n <= {n_max}:")
    print("  the canonical gauge group of the erased fibre over a shadow with r runs")
    print("  is exactly (Z2)^r, abelian.  The C4 frame rotation moves the shadow, so")
    print("  the fibre-preserving canonical group stays (Z2)^r.")
    print("  KS-3 FIRES in its direct form: no forced non-abelian action on a fibre.")


# ---------------------------------------------------------------------------
# Part F: operator transport of the antisymmetric fibre line
# ---------------------------------------------------------------------------

def rho_op(q, c, a, b, z):
    """rho_c(a,b,z) f(x) = psi_c(z + b(x-a)) f(x-a), as a q x q matrix."""
    mat = np.zeros((q, q), dtype=complex)
    for x in range(q):
        mat[x, (x - a) % q] = np.exp(2j * np.pi * c * ((z + b * (x - a)) % q) / q)
    return mat


def weyl_op(q, c, a, b):
    """W_c(a,b) f(x) = psi_c(b(x - a/2)) f(x-a)."""
    inv2 = pow(2, -1, q)
    mat = np.zeros((q, q), dtype=complex)
    for x in range(q):
        phase = (b * ((x - a * inv2) % q)) % q
        mat[x, (x - a) % q] = np.exp(2j * np.pi * c * phase / q)
    return mat


def parity_op(q):
    mat = np.zeros((q, q))
    for x in range(q):
        mat[x, (-x) % q] = 1.0
    return mat


def fourier_op(q, c):
    """Metaplectic operator for w, F[x,y] = q^{-1/2} psi_c(x y) (deposited model)."""
    mat = np.array([[np.exp(2j * np.pi * c * (x * y % q) / q) for y in range(q)]
                    for x in range(q)]) / math.sqrt(q)
    return mat


def odd_part(mat, R):
    return 0.5 * (mat - R @ mat @ R)


def part_F(q=29, c=1, n_max=6):
    print("== F. Operator transport of the erased fibre ==")
    R = parity_op(q)
    F = fourier_op(q, c)

    d1 = rho_op(q, c, 1, 0, 0) - rho_op(q, c, -1, 0, 0)
    w1 = weyl_op(q, c, 1, 0) - weyl_op(q, c, -1, 0)
    print(f"  q={q}: || (rho(X)-rho(X^-1)) - (W(X)-W(-X)) || = {np.abs(d1 - w1).max():.2e}")
    print(f"  R-odd check: || even part of fibre difference || = "
          f"{np.abs(d1 + R @ d1 @ R).max():.2e}  (0 = exactly odd)")

    Finv = np.linalg.inv(F)
    conj = F @ d1 @ Finv
    conj_inv = Finv @ d1 @ F
    wy = weyl_op(q, c, 0, 1) - weyl_op(q, c, 0, -1)
    err_plus = np.abs(conj - wy).max()
    err_minus = np.abs(conj + wy).max()
    print(f"  Frame completion: F (fibre diff) F^-1 = +-(W(Y)-W(-Y)): "
          f"errors {err_plus:.2e} / {err_minus:.2e} (one must vanish)")
    t_plus_ref = 1j * (d1 + wy)
    t_minus_ref = 1j * (d1 - wy)
    t_from_w = 1j * (d1 + conj)
    t_from_w_inv = 1j * (d1 + conj_inv)
    print(f"  Choosing w gives the + axis: error {np.abs(t_from_w - t_plus_ref).max():.2e}")
    print(f"  Choosing w^-1 gives the - axis: error "
          f"{np.abs(t_from_w_inv - t_minus_ref).max():.2e}")
    print("  Thus the C4 orbit spans the two-dimensional real plane <A_X,A_Y>;")
    print("  the unoriented frame set alone does not select T_plus over T_minus.")

    print("  General shadows: R-odd component of antisymmetrised fibre endpoints:")
    for n in range(1, n_max + 1):
        fibres = defaultdict(list)
        for wrd in words(n):
            fibres[shadow_of(wrd)].append(wrd)
        odd_norms = []
        for shd, fib in fibres.items():
            if zero_runs(shd) == 0:
                continue
            ops = np.zeros((q, q), dtype=complex)
            for wrd in fib:
                sign = 1
                run = -1
                prev = None
                for g in wrd:
                    s = B_EMIT[g]
                    if s == 0 and prev != 0:
                        run += 1
                        if g == 1:
                            sign = -sign
                    prev = s
                a, b, z = endpoint(wrd)
                ops += sign * rho_op(q, c, a % q, b % q, z % q)
            odd_norms.append(np.linalg.norm(odd_part(ops, R)))
        nz = sum(1 for v in odd_norms if v > 1e-9)
        print(f"    n={n}: porous-sign fibres = {len(odd_norms)}, "
              f"with non-zero D_minus component: {nz}")
    print("  The fully antisymmetrised fibre line is defined up to one global sign")
    print("  (relabelling every run flips it), so it is a canonical UNORIENTED line.")
    print("  POSITIVE RESIDUE: after freely linearising histories, the endpoint map")
    print("  sends the alternating length-1 fibre line exactly into D_minus.")
    print("  REMAINING GATE: D4 itself outputs residual norms and does not supply this")
    print("  linear re-entry; moreover C4 completion leaves two Hermitian axes.")
    print("  KS-5 does not fire: no SU(2) was assumed anywhere.")


def main():
    part_A()
    part_B()
    part_CD()
    part_E()
    part_F()
    print("Front verdict: KS-frame confirmed, KS-1 fires (twice), KS-2 confirmed")
    print("abelian, KS-3 fires in direct form; the M_L ~ C^2 target is NOT produced")
    print("at this stratum.  New exact positive residue: an alternating history-fibre")
    print("line maps algebraically to D_minus after a new linear re-entry.  C4 leaves")
    print("a two-axis ambiguity, so canonical odd activation is NOT yet derived.")


if __name__ == "__main__":
    main()
