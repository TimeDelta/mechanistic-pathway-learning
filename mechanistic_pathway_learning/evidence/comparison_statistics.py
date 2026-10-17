"""Each comparison's own statistic, carried beside its verdict (decision 9.48, option 3, the owner, 16 October 2026:
"Carrying stat to calibrated proper is the way").

Three sources state a comparison rather than a fact: a FAERS drug-reaction disproportion, a placebo-controlled trial
arm and an IMPC knockout test. Until now each was cut to a verdict (signal, absence or nothing) by a rule the owner
called "reasonable but unprincipled and slightly arbitrary"; option 3 keeps the comparison's statistic so that the
calibrated reliability (9.59) can learn from the training-fold labels what each level of it is worth. This module
computes the statistic in one shape for every source, the STATISTIC_COLUMNS:

  statistic_name                  what the estimate is
  statistic_estimate              the point estimate on a log scale, 0 for no difference, positive for an excess
  statistic_lower_95, _upper_95   its 95 percent interval on the same scale (credibility or confidence, by source)
  statistic_probability_of_excess P(true effect > 0 | the data), under the source's posterior or a normal approximation
  statistic_information           1 / variance of the estimate, the weight a pooled estimate would give it
  statistic_censored              "upper_bound" when the count behind it is only known to be at most the value used

FAERS: the shrunk observed-to-expected ratio of Noren, Hopstadius and Bate 2013 (Stat Methods Med Res 22(1):57-69,
doi 10.1177/0962280211403604, read in PMC6331976 on 17 October 2026), (O + alpha) / (E + alpha), which "can be viewed
as the Bayesian posterior mean of a parameter mu under the assumption that O is Poisson(mu E)-distributed with a Gamma
prior distribution for mu: Gamma(alpha, alpha)", the posterior then being "Gamma (but with parameters O + alpha and
E + alpha)"; with alpha = 1/2 it approximates the WHO programme's information component ("very well approximated by
() with alpha = alpha = 1/2"). E is the count expected under independence, drug reports times reaction reports over
all reports. The estimate is the natural logarithm of the shrunk ratio (the paper uses base 2), the interval and the
probability of excess are read from the Gamma posterior, and the information is 1 / trigamma(O + alpha), the inverse
posterior variance of log mu. The paper's own warning stands: its intervals "account for random but not systematic
variability" and "may exaggerate the precision of an OE ratio in the presence of clustered data points", as FAERS's
duplicates and co-reported drugs are. A count the API does not list is known only to be at most the smallest listed
one: the statistic is then computed at that bound and marked censored, so its estimate, its upper bound and its
probability of excess are upper bounds.

Trial arms: the log risk ratio of the pooled drug arms against the pooled placebo arms, with 1/2 added to every cell
when one of them is 0 (the Haldane-Anscombe correction), its Katz interval, a normal-approximation probability of
excess and the inverse of the Katz variance as its information. The risk difference with its Newcombe interval stays
in the table's arm columns as before.

What this is not: a verdict. Nothing here decides that a comparison is a signal or an absence; that remains the
builders' rule until the soft priors read the statistic (docs/open_decisions.md, 9.48).
"""
from __future__ import annotations

import math

from scipy import special, stats

STATISTIC_COLUMNS = ("statistic_name", "statistic_estimate", "statistic_lower_95", "statistic_upper_95", "statistic_probability_of_excess",
                     "statistic_information", "statistic_censored")
SHRINKAGE_ALPHA = 0.5
CENSORED_UPPER_BOUND = "upper_bound"
FAERS_STATISTIC = "log_shrunk_observed_to_expected"
TRIAL_STATISTIC = "log_risk_ratio"


def expected_count(drug_total: float, reaction_total: float, all_reports: float) -> float:
    """Reports expected to carry both the drug and the reaction if the two were independent."""
    if all_reports <= 0:
        raise ValueError("expected_count needs a positive number of reports")
    return float(drug_total) * float(reaction_total) / float(all_reports)


def shrunk_observed_to_expected(observed: float, expected: float, alpha: float = SHRINKAGE_ALPHA, censored: bool = False) -> dict:
    """The STATISTIC_COLUMNS of a disproportion: the log of (O + alpha) / (E + alpha) and its Gamma posterior."""
    if observed < 0 or expected < 0 or alpha <= 0:
        raise ValueError("observed and expected counts must be at least 0 and alpha positive")
    shape, rate = observed + alpha, expected + alpha
    lower, upper = stats.gamma.ppf([0.025, 0.975], a=shape, scale=1.0 / rate)
    return {"statistic_name": FAERS_STATISTIC, "statistic_estimate": math.log(shape / rate), "statistic_lower_95": math.log(lower),
            "statistic_upper_95": math.log(upper), "statistic_probability_of_excess": float(stats.gamma.sf(1.0, a=shape, scale=1.0 / rate)),
            "statistic_information": 1.0 / float(special.polygamma(1, shape)), "statistic_censored": CENSORED_UPPER_BOUND if censored else ""}


def log_risk_ratio(drug_affected: int, drug_at_risk: int, placebo_affected: int, placebo_at_risk: int) -> dict | None:
    """The STATISTIC_COLUMNS of a trial arm against placebo: the log risk ratio with the Katz interval; None when an arm
    has no one at risk or reports more affected than at risk."""
    if drug_at_risk <= 0 or placebo_at_risk <= 0 or not (0 <= drug_affected <= drug_at_risk) or not (0 <= placebo_affected <= placebo_at_risk):
        return None
    cells = (drug_affected, drug_at_risk - drug_affected, placebo_affected, placebo_at_risk - placebo_affected)
    correction = 0.5 if min(cells) == 0 else 0.0
    drug_cases, drug_total = drug_affected + correction, drug_at_risk + 2 * correction
    placebo_cases, placebo_total = placebo_affected + correction, placebo_at_risk + 2 * correction
    estimate = math.log((drug_cases / drug_total) / (placebo_cases / placebo_total))
    variance = 1.0 / drug_cases - 1.0 / drug_total + 1.0 / placebo_cases - 1.0 / placebo_total
    if variance <= 0:  # every participant affected in both arms: no information on a ratio
        return None
    standard_error = math.sqrt(variance)
    return {"statistic_name": TRIAL_STATISTIC, "statistic_estimate": estimate, "statistic_lower_95": estimate - 1.959963984540054 * standard_error,
            "statistic_upper_95": estimate + 1.959963984540054 * standard_error,
            "statistic_probability_of_excess": float(stats.norm.cdf(estimate / standard_error)), "statistic_information": 1.0 / variance,
            "statistic_censored": ""}


def empty_statistic() -> dict:
    """The STATISTIC_COLUMNS of a row that states no statistic."""
    return {"statistic_name": "", "statistic_estimate": float("nan"), "statistic_lower_95": float("nan"), "statistic_upper_95": float("nan"),
            "statistic_probability_of_excess": float("nan"), "statistic_information": float("nan"), "statistic_censored": ""}
