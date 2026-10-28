# What the certificates contain and what the checkers do

## Model
The AC optimal power flow model of PGLib-OPF v23.07 (the `MODEL.tex` of that release), with the data of each case file
taken as exact decimals: complex bus voltages V_i, nodal active/reactive power balance with shunts and π-model branches
(tap ratio and phase shift), generator active/reactive limits, voltage-magnitude limits, apparent-power limits at both ends
of every rated branch, branch phase-angle-difference limits, and quadratic (convex) generation cost. OPT is the global
minimum of this nonconvex problem. All certificates bind the SHA-256 of the exact case file.

## Upper bound U (`certificates/upper/`, checker `checkers/independent_ub_krawczyk_r12.py`)
The certificate stores an approximate operating point (voltages and the free generator outputs) and a radius. The checker
rebuilds the power-flow equations F(z) = 0 from the case file, forms the Krawczyk operator
K(z̃, Z) = z̃ − C·F(z̃) + (I − C·F′(Z))(Z − z̃) with an approximate inverse C (floating point is used only to *propose* C)
and verifies in interval arithmetic (mpmath, outward rounding) that K(z̃, Z) lies in the interior of Z. By the Krawczyk
theorem the box Z then contains an exact solution of the power-flow equations (existence is all that the upper bound
needs). The checker then evaluates every operational limit and the cost over the whole box Z in interval arithmetic; if
all limits hold on Z, that exact solution is feasible and its cost is at most the upper end of the cost interval, which
must not exceed the declared U. Exit codes: 0 accept; 2 rejected (Krawczyk containment, operational limit or declared
value); 4 unsupported format; 5 input hash mismatch; 6 checker error.

## Lower bound L (`certificates/lower/`)
Two formats, both dual (sum-of-squares / Lagrangian) certificates over the same model:
- `*_lower_chordal.json`, checker `independent_chordal_r13.py`: multipliers for the power-balance, voltage, generator,
  thermal and angle constraints plus positive-semidefinite matrices on the cliques of a chordal extension of the network.
  The checker rebuilds the Lagrangian quadratic form from the case file in interval arithmetic, proves each clique matrix
  PSD after a small diagonal shift (interval LDLᵀ; floating point only proposes the shift), and subtracts the price of the
  shifts and of the residual between the form and the clique matrices over the voltage box. The result is a rigorous lower
  bound; the declared L must not exceed it.
- `*_lower_sos.json`, checker `independent_sparse.py`: a polynomial identity
  f(x) − t = Σ Gram blocks (order 1 or order 2 in the clique variables) + Σ multipliers × constraints + thermal terms
  + Σ polynomial multipliers × equalities + residual, in the real voltage coordinates x = (e, f). Gram and localizing
  matrices are proven PSD (interval LDLᵀ after a priced shift); every residual monomial is priced at its worst case over
  the voltage box. Bus costs with several generators (or, where marked, single generators) enter through affine minorants
  λ·P + ψ(λ) whose constants ψ are computed exactly.

## Independence
The three checkers were written separately from the programs that produced the certificates and import nothing from them;
each reads only the raw `.m` case file and one certificate. See `AI_DISCLOSURE.md` for who wrote and ran what.
