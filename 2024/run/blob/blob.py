#!/usr/bin/env python3


import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
from pyphare.pharesee.run import Run


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')


from tests.diagnostic import all_timestamps



def config(**kwargs):

    def density(x, y):
        return 1.
    
    
    def bx(x, y):
        return kwargs.get("bx", 0.)
    
    
    def by(x, y):
        return kwargs.get("by", 1.)
    
    
    def bz(x, y):
        return 0.0
    
    
    def T(x,y):
        return kwargs.get("Ti", 0.1)
    
    def S(x, l=0.5):
        return 0.5*(1 + np.tanh(x/1))

    def square(x,y, xmin=8, xmax=16, ymin=28, ymax=36):
        sim = ph.global_vars.sim
        l=1
        s1 = S(x-xmin, l=l)
        s2 = 1-S(x-xmax, l=l)
        s3 = S(y-ymin, l=l)
        s4 = 1-S(y-ymax, l=l)
        return s1*s2*s3*s4


    def vx(x,y):
        V0 = 0.5
        return square(x,y)*V0
    
    
    
    def vy(x,y):
        return 0.
    
    
    def vz(x,y):
        return 0.
    
    
    def vthx(x,y):
        return np.sqrt(T(x,y))
    
    
    def vthy(x,y):
        return np.sqrt(T(x,y))
    
    
    def vthz(x,y):
        return np.sqrt(T(x,y))
    
    vvv = {"vbulkx": vx,
           "vbulky": vy,
           "vbulkz": vz,
           "vthx": vthx,
           "vthy": vthy,
           "vthz": vthz }
    
    Simulation(
        time_step=0.01,
        final_time=48.,
        hyper_resistivity=0.02 ,
        cells=(128,128),
        dl=(0.5,0.5),
        diag_options={"format": "phareh5",
                      "options": {"dir": kwargs["diagdir"],
                                  "mode":"overwrite"}
                     }
    )


    MaxwellianFluidModel(bx=bx,
                         by=by,
                         bz=bz,
                         protons={"charge": 1,
                                  "density": density,
                                  "nbr_part_per_cell": 100,
                                  **vvv}
                        )

    ElectronModel(closure="isothermal", Te=kwargs.get("Te", 0.005))


    sim = ph.global_vars.sim
    dt = sim.time_step*100
    timestamps = np.arange(0,sim.final_time, dt)


    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )


    for quantity in ["density", "bulkVelocity", "pressure_tensor"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            )


    for quantity in ['domain']:  # , 'levelGhost', 'patchGhost']:
        ParticleDiagnostics(quantity=quantity,
                            compute_timestamps=timestamps,
                            write_timestamps=timestamps,
                            population_name="protons")


def main():
    from pyphare.cpp import cpp_lib
    import sys
    cpp = cpp_lib()

    Te = float(sys.argv[1])
    Ti = float(sys.argv[2])
    bx = float(sys.argv[3])
    by = float(sys.argv[4])
    
    config(diagdir=f"blob_{Te}_{Ti}_{bx}_{by}", Te=Te, Ti=Ti, bx=bx, by=by)
    
    Simulator(gv.sim).run()
    gv.sim = None




if __name__=="__main__":
    main()
