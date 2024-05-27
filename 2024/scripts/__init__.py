import numpy as np
import matplotlib.pyplot as plt


def scatter_plot(particles, **kwargs):

    default_c = {pop: color  for pop, color in zip(particles, plt.rcParams['axes.prop_cycle'].by_key()['color'][:len(particles)])}
    colors = kwargs.get("colors", default_c)

    if "ax" not in kwargs:
        fig, ax = plt.subplots()
    else:
        ax = kwargs["ax"]
        fig = ax.figure

    for pop in particles:
        x_ = particles[pop].dl[:,0]*particles[pop].iCells[:,0]+particles[pop].deltas[:,0]
        v_ = particles[pop].v[:,0]
    
        im = ax.scatter(x_, v_, s=0.1, c=colors[pop])


def dist_plot(self, **kwargs):
    """
    plot phase space of a particle hierarchy
    """
    import copy

    from pyphare.pharesee.plotting import dist_plot as dp
    #### from .plotting import scatter_plot as sp
    sp=scatter_plot

    usr_lvls = kwargs.get("levels", (0,))
    finest = kwargs.get("finest", False)
    pops = kwargs.get("pop", [])
    time = kwargs.get("time", self.times()[0])
    axis = kwargs.get("axis", ("Vx", "Vy"))
    all_pops = list(self.level(0, time).patches[0].patch_datas.keys())

    plot_type=kwargs.get("plot_type", "pcolor")

    if finest:
        final = finest_part_data(self)
    else:
        final = {pop: None for pop in all_pops}
        for lvl_nbr, level in self.levels(time).items():
            if lvl_nbr not in usr_lvls:
                continue
            for ip, patch in enumerate(level.patches):
                if len(pops) == 0:
                    pops = list(patch.patch_datas.keys())

                for pop in pops:
                    tmp = copy.copy(patch.patch_datas[pop].dataset)

                    if final[pop] is None:
                        final[pop] = tmp
                    else:
                        final[pop].add(tmp)

    if plot_type == "pcolor":
        vmin = kwargs.get("vmin", -2)
        vmax = kwargs.get("vmax", 2)
        dv = kwargs.get("dv", 0.05)
        vbins = vmin + dv * np.arange(int((vmax - vmin) / dv))

        if finest:
            if axis[0] == "x":
                xbins = amr_grid(self, time)
                bins = (xbins, vbins)
            else:
                bins = (vbins, vbins)
            kwargs["bins"] = bins

    if "select" in kwargs:
        for pop, particles in final.items():
            final[pop] = kwargs["select"](particles)

    if plot_type == "pcolor":
        return final, dp(final, **kwargs)
    elif plot_type == "scatter":
        return final, sp(final, **kwargs)