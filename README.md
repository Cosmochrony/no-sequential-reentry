# Operator Doublets, Erased History Fibres, and the No-Sequential-Re-entry Theorem

Fermionic-matter sub-programme no-go note.

The note closes, with exact finite proofs and machine certification, the two deposited
candidates for the missing fermionic weak-isospin carrier of the Cosmochrony corpus.

## Main results

- Abelian carrier commutant: the commutant of the finite Weil module is spanned by the
  two parity projectors, so no non-abelian connected action commutes with the Weil
  action on `V_c`; the canonical Hecke `M_2(C)` doublets of `End(V_c)` have no
  fermionic lift to their own carrier.
- Erased-fibre structure: the oriented Cayley alphabet has stabiliser exactly `C4`;
  the fibre of the canonical projective channel over a `b`-shadow with `r` zero-runs is
  exactly `(Z2)^r`; the canonical inductions of the history space carry only
  near-regular multiplicities, with no privileged `M_2(C)` factor.
- Alternating transport and axis exchange: the flip-antisymmetric fibre line maps
  exactly to the parity-odd operator sector, but the two odd Hermitian axes
  `i(A_X +- A_Y)` are exchanged by a semilinear automorphism of the cocycle-enriched
  Heisenberg--Weil datum, implemented by the deposited antiunitary `RK`; no natural
  selector exists, and the coordinate formula that selects one is conditional on an
  oriented central/Weyl pinning.
- No-sequential-re-entry theorem: the channel acts on amplitudes as a conditional
  expectation `P`, and `Pv = 0 => QPv = 0` for every downstream linear `Q`; no tower of
  sequential projections recovers the erased alternating character.
- Architectural dichotomy: any realisation requires a parallel channel from the
  pre-projection space plus an orientation of `(Z, w)` --- two new axioms, hence an
  extension architecture, not a derivation.

## Reproducibility

The three archived scripts in `scripts/` certify every computational claim and rerun in
about one second in total:

- `front_hist_fibre_doubling.py`: frame stabiliser (exhaustive at `q = 13, 29, 53`),
  fibre bijection and Pell counts (`n <= 10`), `C4` freeness, flip involutions, exact
  alternating transport at `q = 29`;
- `front_axis_selector.py`: enriched-data symmetry table, `RK` intertwining and
  `RK F RK = F^{-1}`, conditional coordinate-pinned selector, single-line support
  obstruction;
- `front_joint_reentry_audit.py`: exact rational two-point fibre algebra, no-re-entry
  implication, exhaustive Heisenberg-automorphism check at `q = 5`.

## Status

Unpublished draft (branch `0.1`).
No release cascade has been triggered.

## Compilation

```bash
bash compile.sh
```

Output: `out/NoSequentialReentry.pdf`.
