"""Cox proportional-hazards survival link (nonparametric tabulated baseline)."""

import json

import numpy as np
import onkos
from onkos.export.registry import get_kernel
from onkos.simulate import _survival_vals

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
COX = "survival_link.nsclc_os_cox"
T = np.linspace(0.0, 208.0, 417)


def test_cox_kernel_registered_and_uses_baseline():
    spec = get_kernel(onkos.load()[COX])
    assert spec.name == "survival_cox_ph"
    assert spec.uses_baseline is True


def test_cox_is_not_auto_selected():
    """The non-default Cox link must not collide with the Weibull OS link."""
    ds = onkos.load()
    tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0)
    # auto-discovery yields exactly the default OS (Weibull) and PFS links
    assert set(tr.survival) == {"OS", "PFS"}


def test_cox_curve_is_valid_survival():
    ds = onkos.load()
    tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0,
                        t=T, survival_link=COX)
    s = tr.os_curve
    assert s is not None
    assert abs(s[0] - 1.0) < 1e-6        # starts at full survival
    assert np.all(np.diff(s) <= 1e-9)    # monotone non-increasing
    assert np.all((s >= 0) & (s <= 1))


def test_cox_covariate_direction():
    ds = onkos.load()
    cox = ds[COX]
    spec = get_kernel(cox)
    better = spec.analytic(T, _survival_vals(cox, -0.5))  # deep shrinkage
    worse = spec.analytic(T, _survival_vals(cox, +0.5))   # growth
    assert better[200] > worse[200]


def test_weibull_vs_cox_disagree_on_same_metric():
    """Survival-model choice is itself an uncertainty axis."""
    ds = onkos.load()
    w = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0, t=T)
    c = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0, t=T,
                       survival_link=COX)
    assert np.max(np.abs(w.os_curve - c.os_curve)) > 0.05


def test_cox_baseline_travels_in_exports():
    ds = onkos.load()
    from onkos.export import to_virtual_trial_json

    vt = json.loads(to_virtual_trial_json(ds[COX]))
    assert vt["baseline_survival"]["times"][0] == 0
    assert len(vt["baseline_survival"]["survival"]) == len(vt["baseline_survival"]["times"])
