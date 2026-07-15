# Operator Doublets, Erased History Fibres, and the No-Sequential-Re-entry Theorem

Fermionic-matter sub-programme no-go note.

The note closes, with exact finite proofs and machine certification, the two deposited
candidates for the missing fermionic weak-isospin carrier of the Cosmochrony corpus.

## Main results

- Abelian carrier commutant: the commutant of the finite Weil module is spanned by the
  two parity projectors, so no non-abelian connected action commutes with the Weil
  action on `V_c`; the canonical Hecke `M_2(C)` doublets of `End(V_c)` have no
  fermionic lift to their own carrier.
- Erased-fibre structure: the symmetric Cayley alphabet has stabiliser exactly `C4`;
  the fibre of the canonical projective channel over a `b`-shadow with `r` zero-runs is
  canonically a torsor under the abelian run-flip group `(Z2)^r`; the canonical
  inductions of the history space carry only near-regular multiplicities, with no
  privileged `M_2(C)` factor.
- Alternating transport and axis exchange: the flip-antisymmetric fibre line maps
  exactly to the parity-odd operator sector, but the two odd Hermitian axes
  `i(A_X +- A_Y)` are exchanged by a semilinear automorphism of the cocycle-enriched
  Heisenberg--Weil datum, implemented by the deposited antiunitary `RK`; no natural
  selector exists, and the coordinate formula that selects one is conditional on an
  oriented central/Weyl pinning.
- No-sequential-re-entry theorem: the canonical record quotient satisfies
  `pi v = 0 => Q pi v = 0` for every downstream linear `Q`, and the channel's
  conditional expectations annihilate the mu-centred fibre differences for every
  non-degenerate measure; no tower of sequential projections recovers the erased
  alternating character.
- Architectural dichotomy (history-fibre route): realisation requires a parallel
  channel from the pre-projection space plus an orientation of `(Z, w)` --- two new
  ingredients, or one new axiom supplying both; an independent new carrier outside the
  deposited corpus remains a distinct, unconstrained extension architecture.

## Reproducibility

The three archived scripts in `scripts/` certify every computational claim and rerun in
about one second in total:

- `front_hist_fibre_doubling.py`: frame stabiliser (exhaustive at `q = 13, 29, 53`),
  fibre bijection and Pell counts (`n <= 10`), `C4` freeness, flip involutions, exact
  alternating transport at `q = 29`;
- `front_axis_selector.py`: enriched-data symmetry table, `RK` intertwining and
  `RK F RK = F^{-1}`, conditional coordinate-pinned selector, single-line support
  obstruction;
- `front_joint_reentry_audit.py`: exact rational two-point fibre algebra (uniform and
  non-uniform conditional expectations, mu-centred kernel line), no-re-entry
  implication, exhaustive Heisenberg-automorphism check at `q = 5`.

## Status

Published working paper (v1.0), concept DOI [10.5281/zenodo.21381089](https://doi.org/10.5281/zenodo.21381089).


## Compilation

```bash
bash compile.sh
```

Output: `out/NoSequentialReentry.pdf`.
