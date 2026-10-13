"""Resting pH per compartment, so a proton deviation can be read on the log scale it belongs on.

The owner, 13 October 2026: "for pH, ratio against HCO3 doesn't quite seem quite right because the ratio could be 10
and still only have 10 H+ ions, which would still limit a reaction that needs H+ despite the ratio being relatively
high ... why not use true pH on the log scale so it fits nicely into equations?" He is right, and it is the reason the
other cell-state quantities are ratios and this one is not. A ratio is dimensionless and so needs no scale, but it also
says nothing about how many protons there are, and a reaction consuming H+ is limited by the count and not by the
proportion. pH is the log-scale absolute, enzyme pH-activity curves are bell-shaped in pH rather than in [H+], and
having an absolute leaves Henderson-Hasselbalch available for the question these values will be needed for next, a
drug's ionisation against a compartment's pH and so the trapping of bases in the lysosome.

With these, compartment_ph.compartment_ph_weight reads

    pH = resting pH - log10(1 + sensitivity * deviation)

where sensitivity * deviation is the fractional change in [H+], so the sensitivity stays learned (the model's
deviations are in arbitrary units) while the resting value sets where the compartment sits.

**Provenance, per value, because they are not all of one kind.** The values marked `review` are from a curated
database entry read on 13 October 2026, whose primary source is the review it names; the ones marked `reasoned` are
not in that entry and the reasoning is given with each.

- Primary: Casey JR, Grinstein S, Orlowski J, "Sensors and regulators of intracellular pH", Nature Reviews Molecular
  Cell Biology 11(1):50-61, 2010, doi:10.1038/nrm2820. Not open access; the values were read from the BioNumbers entry
  below rather than from the paper, which is recorded here rather than implied.
- Also cited by that entry: Theillet FX et al., "Physicochemical properties of cells and their effects on
  intrinsically disordered proteins", Chemical Reviews 114(13):6661-714, 2014, doi:10.1021/cr400695p, figure 2B.
- Read at: `bionumbers.hms.harvard.edu/bionumber.aspx?id=115761` ("Average pH values of organelles and compartments in
  mammalian cells"), 13 October 2026.
- **Licence: unverified.** The entry page states none, and the site's own About page states none either: it gives
  citation guidance only ("When using a specific entry from the database it is recommended that you also specify the
  BioNumbers 6 digit ID"). Routes tried: the entry page and `bionumbers.hms.harvard.edu/aboutus.aspx`. A third route,
  the primary review itself, is behind a publisher paywall. The address is given here for the owner to open. What is
  recorded below is a dozen physical constants with their primary citation, not a copy of a database.
"""
from __future__ import annotations

# compartment letter -> (resting pH, how the value was arrived at)
RESTING_PH_BY_COMPARTMENT: dict[str, tuple[float, str]] = {
    "c": (7.2, "review: cytosol"),
    "n": (7.2, "review: nucleus, which the entry gives as the cytosol's value 'due to permeability'"),
    "r": (7.2, "review: endoplasmic reticulum"),
    "m": (8.0, "review: mitochondria"),
    "x": (7.0, "review: peroxisomes"),
    "l": (4.7, "review: lysosomes"),
    "g": (6.4, "review, midpoint: the entry separates cis-Golgi cisterns at 6.7 from the trans-Golgi network at 6.0 "
               "and the graph has one Golgi compartment, so the midpoint is taken and the endpoints recorded here"),
    "e": (7.4, "reasoned: arterial blood and interstitial fluid, which the entry does not cover; it is the value the "
               "extracellular compartment of a mammalian reconstruction is built at"),
    "i": (7.2, "reasoned: the mitochondrial intermembrane space is continuous with the cytosol through porins, so it "
               "is held at the cytosol's value. It is in fact somewhat more acidic, since the respiratory chain pumps "
               "protons into it, and no value for it is in the entry"),
    "v": (4.7, "reasoned: the vacuole is not a compartment of mammalian cells and reaches this graph only through "
               "Reactome; it is held at the lysosome's value, the acidic degradative compartment it corresponds to"),
}
# the Golgi endpoints, kept because a build that splits cis from trans should use them rather than the midpoint
CIS_GOLGI_PH, TRANS_GOLGI_NETWORK_PH = 6.7, 6.0
# values the entry gives from fluorescence imaging of HeLa cells, which differ from the review's and are kept for a
# reader comparing them: cytoplasm 7.4, mitochondria 8.0, endoplasmic reticulum 7.5, Golgi 6.6
HELA_IMAGED_PH = {"c": 7.4, "m": 8.0, "r": 7.5, "g": 6.6}


def resting_ph(compartment: str) -> float | None:
    """The resting pH of a compartment letter, or None where the graph names a compartment this does not cover (a
    transport reaction's "c;m", or a compartment a later source adds)."""
    entry = RESTING_PH_BY_COMPARTMENT.get(str(compartment))
    return None if entry is None else entry[0]


def resting_ph_per_compartment(compartments) -> list[float | None]:
    return [resting_ph(compartment) for compartment in compartments]
