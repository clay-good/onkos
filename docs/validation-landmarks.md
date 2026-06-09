# Scientific landmark validation

**NOT FOR CLINICAL USE.** This catalogue validates *model fidelity*, not clinical accuracy.

Onkos validates its reference kernels along **two independent axes**:

1. **Export round-trip** ([`test_reference.py`](../tests/test_reference.py)) — proves every
   exported artifact (NONMEM, SBML, PharmML, MathML, rxode2, Pumas, …) reproduces the
   reference kernel to numeric tolerance. This guards *projection consistency*: the
   exports faithfully encode the kernel.
2. **Scientific landmarks** ([`test_landmarks.py`](../tests/test_landmarks.py), this document)
   — proves each kernel reproduces the characteristic, analytically-derivable property of
   the *published model* it claims to implement. This guards *scientific fidelity*: the
   kernel is the model it says it is.

The two axes are complementary. A kernel can pass the round-trip (exports agree with it)
while still being the wrong model — landmark validation is what catches that. Each landmark
below is a quantitative property derived from the model's own equations; this is the honest
reading of the spec's §9 directive to compare against "published example simulations" — the
landmark *is* the published property, so no digitized data is invented.

## Growth laws

| Kernel | Landmark | Condition | Source |
|---|---|---|---|
| `growth_exponential` | Doubling time | `t = ln2/kg` ⇒ `V = 2·V₀` | definitional |
| `growth_logistic` | Inflection (max growth rate) | `V = Vmax/2` | Verhulst logistic |
| `growth_gompertz` | Inflection (max growth rate) | `V = Vmax/e` | Gompertz / Laird |

## Drug effect & TGI

| Kernel | Landmark | Condition | Source |
|---|---|---|---|
| `claret_tgi` | Initial perturbed log-slope | `d(ln V)/dt\|₀ = kL − kD·E` | Claret 2009 |
| `norton_simon` | Stationary tumor (all V) | `E = g/k` ⇒ `dV/dt = 0` | Norton 2005 |
| `biexp_tgi` | Nadir time | `t* = ln(ks·E/kg)/(kg+ks·E)`, a minimum | Stein 2008 / Wang 2009 |
| `simeoni_exp_linear` | Exp→linear transition | rate → `λ0·w` (small), → `λ1` (large) | Simeoni 2004 |
| `simeoni_tgi` | Tumor-static concentration | `c* = λ0/k2` ⇒ proliferating compartment stationary | Simeoni 2004 |

## Survival

| Kernel | Landmark | Condition | Source |
|---|---|---|---|
| `survival_weibull_ph` | Median survival (x=0) | `t = scale·(ln2)^(1/shape)` ⇒ `S = 0.5` | Weibull / Cox PH |
| `survival_weibull_ph` | Proportional hazards | `S_x = S₀^exp(β·x)` | Cox PH |
| `survival_cox_ph` | Baseline reproduction | `x=0` ⇒ `S(t) = S₀(t)`; median at `S₀ = 0.5` | Cox PH |

## Exposure-response

| Kernel | Landmark | Condition | Source |
|---|---|---|---|
| `er_emax` | Half-maximal effect | `C = EC50` ⇒ `E = Emax/2` | Holford 1981 |
| `er_sigmoid_emax` | Half-maximal (any Hill γ) | `C = EC50` ⇒ `E = Emax/2` | Hill / Holford 1981 |
| `er_power` | Scale-free | `E(2C)/E(C) = 2^θ` | power law |
| `ivive_power` | Scale-free | potency(2C)/potency(C) = `2^power` | Rocchetti 2007 |

## Immuno-oncology (hypothesis tier — DO NOT USE FOR PREDICTION)

| Kernel | Landmark | Condition | Source |
|---|---|---|---|
| `io_tumor_immune` | Immune homeostasis | tumor=0 ⇒ effectors → `s/δ` | Kuznetsov 1994 |

---

These checks run in CI on every push. Adding a kernel without a corresponding landmark is
discouraged: the landmark is what certifies the kernel is scientifically the model it names.
All parameter values used in the tests are illustrative — landmark validation certifies
*structural* fidelity (the model behaves as published), not that any specific parameter set
is clinically verified.
