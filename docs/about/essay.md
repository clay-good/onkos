# Why model-and-context-selection risk is the load-bearing idea

Oncology drug development runs on a chain of models. Drug exposure drives a
tumor-growth-inhibition (TGI) model; the TGI model produces a metric (a growth
rate constant, a week-8 change, a depth of response); the metric drives a
survival link that predicts overall survival. That prediction gates go/no-go
decisions worth hundreds of millions of dollars. The chain is standardized,
published, and reused constantly.

The reuse is where it fails. A TGI model is fit in one context: one drug, one
tumor type, one line of therapy, one trial. Then it is transported to another
context, often without comment, where its predictive validity is simply unknown.
The resistance and kill-rate terms that matter most for the survival tail are the
worst-identified parameters in the whole chain; coefficients of variation near
90% are routine because resistance is barely observable over a short trial. A
point estimate copied from a table hides all of that.

Onkos makes two normally-invisible facts first-class and machine-readable.

**Derivation context and transportability.** Every record states the exact
context it came from and the boundary within which it has actually been
validated. When a simulation crosses that boundary, the result is auto-tiered to
D and a warning is attached. You cannot get an A-looking survival forecast from a
model validated only on a different tumor type. This is the same discipline as
Nidus's per-parameter confidence tier and Hypnos's applicability envelope,
specialized for the failure mode oncology actually has.

**Propagated, worst-input-wins confidence.** A composed forecast is only as
trustworthy as its least-validated component or its furthest extrapolation, so
the composed tier is the worst of its parts, and that propagated tier travels
into every export as RDF.

The headline feature follows directly. The virtual-trial divergence view takes a
tumor context and a drug-effect size, overlays the simulated tumor-size and
population OS curves across every eligible TGI model, greys out the ones whose
transportability envelope the context violates, and reports how much the survival
prediction depends on the model choice. That number is model-selection risk made
measurable. Unquantified, it is exactly the risk that sends drugs into doomed
phase-3 trials.

Onkos does not predict a person's prognosis and does not rank therapies. It
reports what published models said about trials, with their uncertainty intact,
and it refuses to let a context-specific fit masquerade as a general truth. That
refusal is the contribution.
