import os

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.run import Run

from tests.diagnostic import all_timestamps

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from pyphare.pharesee import plotting
from pathlib import Path

mpl.use("Agg")


def densityMain(x):
    return 1.0

def densityBeam(x):
    return 0.2

def bx(x):
    return 1.0

def by(x):
    return 0.0

def bz(x):
    return 0.0

def vB(x):
    return 2.5

def v0(x):
    return 0.0

def vthMain(x):
    return np.sqrt(0.05)

def vthBeam(x):
    return np.sqrt(1.0)


vMain = {
    "vbulkx": v0,
    "vbulky": v0,
    "vbulkz": v0,
    "vthx": vthMain,
    "vthy": vthMain,
    "vthz": vthMain,
    }

vBulk = {
    "vbulkx": vB,
    "vbulky": v0,
    "vbulkz": v0,
    "vthx": vthBeam,
    "vthy": vthBeam,
    "vthz": vthBeam,
    }

def config(**kwargs):

    sim = ph.Simulation(
        final_time=100,
        time_step=0.01,
        boundary_types="periodic",
        cells=200,  # what is the k of the most unstable ES mode ? gary has the answer ?
        dl=0.2,
        hyper_resistivity=0.01,
        diag_options={
            "format": "phareh5",
            "options": {"dir": kwargs["diagdir"],
                        "mode": "overwrite"},
        },
    )


    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        main={"charge": 1, "density": densityMain, "nbr_part_per_cell": 200, **vMain},
        beam={"charge": 1, "density": densityBeam, "nbr_part_per_cell": 200, **vBulk},
    )

    ph.ElectronModel(closure="isothermal", Te=0.00025)

    timestamps = np.arange(0, sim.final_time, 0.1)

    for quantity in ["B", "E"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    for pop_name in ["main", "beam"]:
        ph.ParticleDiagnostics(
            quantity="domain",
            population_name=pop_name,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )
    return sim


def main():
    Simulator(config(diagdir="ib1")).run()


if __name__ == "__main__":
    main()
