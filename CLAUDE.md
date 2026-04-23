# CausalTransPrediction — Project Brief

## What this project is

A simulation study for a **high-impact biostatistics paper** combining three methodological threads:

1. **Counterfactual (causal) risk prediction** — estimating E(Y^a | P, S=0), the mean counterfactual outcome given risk predictors P in a target population.
2. **Transportability** — using external cohort data (S=1,...,K) to inform inference about a target population (S=0) that lacks outcome/treatment data.
3. **Meta-analysis** — combining DR estimates across K heterogeneous databases via IVW fixed-effects pooling. Now implemented in `simulation_meta.ipynb`.

The paper is co-authored by at least Fuyu and Sinclair (sinclaircarr on git).

---

## Theoretical setup (`note.md`)

### Data structure
Each observation: O = (X, P, L, S, A, Y)
- X: baseline covariates; P ⊂ X are the risk predictors
- L: confounders (cohort-only)
- S: cohort indicator (S=0 target population, S=1,...,K external cohorts)
- A: treatment; Y: outcome

### Identification (Theorem 1)
Under A1–A5 (consistency, no unmeasured confounding, positivity of treatment, transportability, positivity of membership):

```
E(Y^a | P, S=0) = E[ E{ E(Y | A=a, X^(k), L^(k), S=k) | X^(k), S=0 } | P, S=0 ]
```

Three-step plug-in algorithm:
1. Fit outcome model m(X,L) = E(Y | A=a, X, L, S=k) in cohort k
2. Marginalise over L within cohort k → ψ(X)
3. Regress ψ(X) on P in target population → β̂

### Influence function / DR estimator (Theorem 2)
Linear model: E(Y^0 | P, S=0) = β⊤P

The influence function φ_β(O) has a **doubly-robust** structure with four nuisance functions:
- m(X, L): outcome model (E[Y | A=0, X, L, S=1])
- π(X, L): propensity score (Pr[A=1 | X, L, S=1])
- ψ(X): transportability bridge (E[m(X,L) | X, S=1])
- ρ₀(X), ρ₁(X): cohort membership probabilities

**Double robustness (Corollary):** β̂_DR is unbiased if:
- Either m or π is correctly specified, AND
- Either ψ or ρ₀ is correctly specified

---

## Simulation setup (`latex_code.tex` + `simulation.ipynb`)

### DGP (n=5,000 per replication)
- X₁, X₂ ~ N(0,1)
- S | X ~ Bernoulli(expit(-0.5 + 0.4·X₁ + 0.4·X₂))  →  ~57% cohort, ~43% target
- L = 0.5·X₁ + 0.5·X₂ + ε_L  (cohort only)
- A | X₁, L, S=1 ~ Bernoulli(expit(-1 + 0.5·X₁ + L))
- Y = 1 + 0.5·X₁ + 0.3·X₂ + 0.8·L + 1.5·A + ε_Y  (cohort only)

**True parameters:** β* = (β₀*, β₁*, β₂*) = (1.0, 0.9, 0.7)

True nuisance functions:
- m* = 1 + 0.5·X₁ + 0.3·X₂ + 0.8·L
- π* = expit(-1 + 0.5·X₁ + L)
- ψ* = 1 + 0.9·X₁ + 0.7·X₂
- ρ₁* = expit(-0.5 + 0.4·X₁ + 0.4·X₂)

### Monte Carlo design
- N_sim = 500 replications, n = 2,000 per replication (~1,150 cohort / ~850 target)
- All 2⁴ = 16 combinations of correct/misspecified (m, π, ψ, ρ)
- Misspecification strategy: omit key covariates (e.g., omit L from m, omit X₁ from ψ)
- Metrics: Bias and RMSE for β̂₁ and β̂₂

### Estimators compared
- **Plug-in** (3-step OLS chain)
- **DR** (doubly-robust, solves estimating equation P_n φ_β = 0)

---

## Key results (`table_results.tex`, `simulation_results.png`)

DR estimator confirms double robustness:
- Remains near-zero bias whenever (m OR π correct) AND (ψ OR ρ correct)
- Notably: m wrong + π correct + ψ correct → DR near zero; plug-in badly biased
- When both robustness conditions fail, DR also fails (as expected)
- Plug-in is only unbiased when m and ψ are both correct

---

## Meta-analysis extension (`simulation_meta.ipynb`)

Three fully heterogeneous databases, each pairing a target population with a unique external cohort. All identify β* = (1.0, 0.9, 0.7).

### DGP design principle

Each database has a unique covariate Z^(k) ~ N(0,1) **independent of (X₁,X₂)**. Z^(k) affects BOTH the selection mechanism S^(k) AND the outcome Y — it is a genuine member of X^(k), needed for transportability. Making Z^(k) independent of (X₁,X₂) ensures:
- Slopes β₁*, β₂* are analytically verified (Z^(k) contributes zero to slopes since E(Z^(k)|X₁,X₂)=0)
- Z^(k) in the S model is a genuine misspecification target for ρ
- β₀* = 1.0 requires one numerical calibration: α₀^(k) = 1 − c_Z^(k)·E[Z^(k)|S^(k)=0]

### DGP summary

| | DB1 | DB2 | DB3 |
|---|---|---|---|
| Extra covariate Z^(k) | Z₃ ~ N(0,1) | Z₄ ~ N(0,1) | Z₅ ~ N(0,1) |
| L structure | 0.8X₁+0.2X₂+0.6Z₃ | 0.5X₁+0.5X₂+0.4Z₄ | 0.4X₁+0.8X₂+0.5Z₅ |
| γ^(k) (L coef in Y) | 0.50 | 0.80 | 0.60 |
| α₁^(k), α₂^(k) (direct X₁,X₂ in Y) | 0.50, 0.60 | 0.50, 0.30 | 0.66, 0.22 |
| α_Z^(k) (direct Z in Y) | 0.30 | 0.20 | 0.40 |
| c_Z^(k) = α_Z+γl_Z | 0.60 | 0.52 | 0.70 |
| Treatment effect τ^(k) | 1.4 | 1.8 | 2.2 |
| A depends on | X₁, L | X₂, L | X₁, X₂, L |
| S^(k) model | expit(-0.5+0.3X₁+0.2X₂+0.03Z₃) | expit(-0.6+0.4X₁+0.3X₂+0.03Z₄) | expit(-0.4+0.2X₁+0.5X₂+0.03Z₅) |
| α₀^(k) | calibrated numerically | calibrated numerically | calibrated numerically |

True ψ^(k)* functions:
- ψ¹*(X₁,X₂,Z₃) = α₀¹ + 0.9X₁ + 0.7X₂ + 0.60·Z₃
- ψ²*(X₁,X₂,Z₄) = α₀² + 0.9X₁ + 0.7X₂ + 0.52·Z₄
- ψ³*(X₁,X₂,Z₅) = α₀³ + 0.9X₁ + 0.7X₂ + 0.70·Z₅

### Misspecification per database (uniform)
- m misspecified: omit L (use X₁,X₂,Z only)
- π misspecified: omit L
- ψ misspecified: omit Z (use X₁,X₂ only) — meaningful because c_Z^(k) ≠ 0
- ρ misspecified: omit Z (use X₁,X₂ only) — meaningful because Z^(k) is in S^(k)

### Sandwich variance + IVW pooling
Per-database sandwich: Var(β̂^(k)) = Â⁻¹BÂ⁻¹ where Â = Σ_{S=0} DDᵀ, B = Σᵢ sᵢsᵢᵀ.
IVW fixed-effects: β̂_MA = (Σ_k Σ̂^(k)⁻¹)⁻¹ · Σ_k Σ̂^(k)⁻¹ β̂^(k).

### Simulation parts
- **Part A**: 2⁴ = 16 misspecification scenarios × 3 databases (confirms per-DB double robustness)
- **Part B**: efficiency comparison — individual DB vs pooled (all nuisances correct)
- **Part C**: robustness — one DB misspecified, pooled vs individual

---

## Files

| File | Purpose |
|------|---------|
| `note.md` | Full theoretical write-up with proofs (LaTeX source) |
| `latex_code.tex` | Simulation section for the paper (single-DB DGP, estimators, MC design) |
| `simulation.ipynb` | Single-database simulation: 16 misspecification scenarios, plug-in vs DR |
| `simulation_meta.ipynb` | 3-database meta-analysis simulation: per-DB DR + IVW pooling |
| `table_results.tex` | LaTeX table of MC results for single-DB simulation |
| `simulation_results.png` | Visual summary of single-DB simulation results |
| `convergence.png` | Convergence diagnostic plot |
| `psi_rho_robustness.png` | Robustness plot for ψ and ρ nuisance misspecification |

---

## Open items / known issues

- Assumption A2 wording may need revision — Fuyu flagged it as potentially too strong (red comment in `note.md`)
- The proof of φ_β in `note.md` has a minor LaTeX artifact around line 272 (duplicated expression — copy-paste typo, not a theoretical error)
- Random-effects meta-analysis (allowing heterogeneous β* across databases) not yet implemented
- `latex_code.tex` covers single-DB simulation only; meta-analysis simulation section needs to be written for the paper
