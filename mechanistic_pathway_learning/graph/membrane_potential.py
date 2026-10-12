"""The neuronal membrane potential, computed from the ion pools rather than learned from the reactions.

Why. In the graph as built, `VM_c` received an edge from each of the 405 reactions that move net charge across the
plasma membrane, and the encoder averages a node's incoming messages within a relation, so each one arrived divided by
405 and nothing reached the node (docs/label_source_review.md, sections 22 and 23; the median field at `VM_c` was
2.9e-15). `VM_c` also carries none of the 138 node descriptor columns, so there was nothing else to read it by. The
owner, 12 October 2026: "Why not just change the VM_c to e.g. the Nernst equation", then "Fix the electrical layer".

What this computes. The potential of a membrane with several permeant ions is not the Nernst potential, which is one
ion's equilibrium potential. For a passive membrane at rest it is the parallel-conductance (Millman) equation

    V_m = sum_ion ( g_ion / g_total ) * E_ion,      E_ion = (RT / z_ion F) * ln( [ion]_out / [ion]_in )

and because the encoder's node state is a unitless deviation from rest and not a concentration, the law is used in its
linearised form, which is exact to first order around the resting state:

    change in V_m = sum_ion ( g_ion / g_total ) * (RT / z_ion F) * ( fractional change in [ion]_out
                                                                    - fractional change in [ion]_in )

So each ion pool feeds `VM_c` one edge whose weight is that ion's conductance share times RT/z F, with the sign of the
side it sits on. The weights are fixed numbers from a cited source and are not fitted, which is the point: the
aggregation stops being a learned mean over 405 edges and becomes a stated physical law over six.

Why the parallel-conductance equation and not Goldman-Hodgkin-Katz. GHK weights each ion by a permeability, and no
measured permeability ratio for a mammalian central neuron was found (docs/membrane_potential_constants.md names the
routes tried). The source used here states conductances, not permeabilities, and says so: "Leak currents were
calculated using Ohm's Law, I=g(Vm-E)". Putting its conductances into GHK as if they were permeabilities gives
-33.2 mV against the -72.6 mV it reports, a 39 mV error, so the two laws are not interchangeable here and the law has
to be the one the constants were measured under.

What this cannot claim. The 7 reactions `VM_c` gates are all also sources of charge movement, so the layer is seven
two-cycles through the potential, which is an excitable membrane: regenerative sodium current is positive feedback and
its fixed point is unstable. A model run to a fixed point under a contraction bound cannot show excitability. The
registered claim is therefore the subthreshold steady state, and the contraction bound is part of that claim rather
than a numerical convenience (docs/label_source_review.md 23.3).

The constants, their quotations and where each was read are in docs/membrane_potential_constants.md; the source rows
are in docs/data_sources.md.
"""
from __future__ import annotations

import math

# Düsterwald, Currin, Burman, Akerman, Kay and Raimondo 2018, eLife 7:e39575 (doi 10.7554/eLife.39575, PMID 30260315),
# Table 1, read from the Europe PMC full text on 12 October 2026 under CC BY 4.0. A pump-leak model of a rodent central
# neuron at 37 degrees Celsius. Its own reported steady state is V_m -72.6 mV, E_K -95.1 mV and E_Cl -83.8 mV.
GAS_CONSTANT_JOULES_PER_KELVIN_MOLE = 8.31446
FARADAY_CONSTANT_COULOMBS_PER_MOLE = 96485.33
TEMPERATURE_KELVIN = 310.15
RESTING_CONCENTRATION_MILLIMOLAR = {  # (ion, compartment letter) -> mM; "c" cytosolic, "e" extracellular
    ("Na+", "e"): 145.0, ("Na+", "c"): 14.0,
    ("K+", "e"): 3.5, ("K+", "c"): 122.9,
    ("Cl-", "e"): 119.0, ("Cl-", "c"): 5.2,
}
LEAK_CONDUCTANCE_MICROSIEMENS_PER_SQUARE_CENTIMETRE = {"Na+": 20.0, "K+": 70.0, "Cl-": 20.0}
ION_CHARGE = {"Na+": 1, "K+": 1, "Cl-": -1}
REPORTED_RESTING_POTENTIAL_MILLIVOLTS = -72.6  # Table 1, row Vm, of the same source
# Cytosolic Ca2+ is in the layer as a signal (transmitter release, calcium-activated channels) and not as a driving
# force: its equilibrium potential lies so far from rest that its relative change moves the potential little, and the
# voltage equation above is derived for monovalent ions (docs/label_source_review.md 23.2). It therefore gets no
# driving-force edge, and no extracellular neuronal calcium pool is added for one.
IONS_IN_THE_VOLTAGE_EQUATION: tuple[str, ...] = ("Na+", "K+", "Cl-")


def thermal_voltage_millivolts() -> float:
    """RT/F in millivolts at the source's temperature: 26.73 mV at 37 degrees Celsius."""
    return GAS_CONSTANT_JOULES_PER_KELVIN_MOLE * TEMPERATURE_KELVIN / FARADAY_CONSTANT_COULOMBS_PER_MOLE * 1000.0


def nernst_potential_millivolts(ion: str) -> float:
    """The ion's equilibrium potential at the resting concentrations, (RT/zF) ln([out]/[in])."""
    return (thermal_voltage_millivolts() / ION_CHARGE[ion]) * math.log(
        RESTING_CONCENTRATION_MILLIMOLAR[(ion, "e")] / RESTING_CONCENTRATION_MILLIMOLAR[(ion, "c")])


def conductance_share(ion: str) -> float:
    """The ion's share of the resting leak conductance: 0.636 for K+, 0.182 each for Na+ and Cl-."""
    return LEAK_CONDUCTANCE_MICROSIEMENS_PER_SQUARE_CENTIMETRE[ion] / sum(
        LEAK_CONDUCTANCE_MICROSIEMENS_PER_SQUARE_CENTIMETRE[name] for name in IONS_IN_THE_VOLTAGE_EQUATION)


def passive_resting_potential_millivolts() -> float:
    """The parallel-conductance potential of these constants, with no pump current: -64.4 mV.

    The source reports -72.6 mV, and the 8.2 mV between them is its electrogenic Na+/K+-ATPase, which its model carries
    as a current of its own. The difference is a baseline and not a response coefficient: with the conductances fixed,
    the pump reaches the potential only by moving the ion pools, and the encoder's state is a deviation from rest, so
    the offset cancels and the coefficients below are unaffected."""
    return sum(conductance_share(ion) * nernst_potential_millivolts(ion) for ion in IONS_IN_THE_VOLTAGE_EQUATION)


def driving_force_coefficients_millivolts() -> dict[tuple[str, str], float]:
    """(ion, compartment letter) -> dV_m per unit fractional change in that pool, in millivolts.

    The extracellular side of a cation raises the potential and its cytosolic side lowers it; chloride is the reverse,
    because its charge is negative. The magnitudes are the ion's conductance share times RT/zF, so potassium carries
    3.5 times sodium or chloride, which is the fact a weight of 1 on every edge threw away."""
    coefficients = {}
    for ion in IONS_IN_THE_VOLTAGE_EQUATION:
        outward = conductance_share(ion) * thermal_voltage_millivolts() / ION_CHARGE[ion]
        coefficients[(ion, "e")] = outward
        coefficients[(ion, "c")] = -outward
    return coefficients


def driving_force_edge_weights() -> dict[tuple[str, str], tuple[float, float]]:
    """(ion, compartment letter) -> (sign, weight) for the pool's edge into VM_c.

    The weights are the coefficients above rescaled to sum to 1 in absolute value, because the graph carries no units
    and the encoder reads an edge as sign x weight: only the ratios between the six are physical, and a total of 1
    keeps the node's summed input on the scale every other node's in-degree mean has. Multiplying all six by one
    constant is a change of the gain the relation learns, not of the physics."""
    coefficients = driving_force_coefficients_millivolts()
    total = sum(abs(value) for value in coefficients.values())
    return {key: (math.copysign(1.0, value), abs(value) / total) for key, value in coefficients.items()}
