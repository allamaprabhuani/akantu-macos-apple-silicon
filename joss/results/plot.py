#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

# needed to generate the plots on jed
import matplotlib

# matplotlib.use("TKAgg")
import matplotlib.pyplot as plt

# Same font as JOSS
plt.rcParams["font.sans-serif"] = "cmss10"

# Loading data
plots = {
    "elastic gcc v5": {
        "prefix": "timmings_",
        "material": "elastic",
        "compiler": "gcc",
        "suffix": "_jed_v5.0.4"
    },
    "cohesive gcc v5": {
        "prefix": "timmings_",
        "material": "cohesive",
        "compiler": "gcc",
        "suffix": "_jed_v5.0.4"
    },
    "elastic gcc v4": {
        "prefix": "timmings_",
        "material": "elastic",
        "compiler": "gcc",
        "suffix": "_jed_v4.0.1"
    },
    "cohesive gcc v4": {
        "prefix": "timmings_",
        "material": "cohesive",
        "compiler": "gcc",
        "suffix": "_jed_v4.0.1"
    },
}

fig, ax = plt.subplots(figsize=(4.5, 4))
# fig, ax = plt.subplots(1, 1)

plotting = "TTS"

handles = []
for plot_name, data in plots.items():
    data["df"] = pd.read_csv(
        f"""{data["prefix"]}{data["material"]}_{data["compiler"]}{data["suffix"]}.csv""",
        sep=",",
        skipinitialspace=True,
    )

    df = data["df"]
    step = df["solve_step"] * df["solve_step nb_rep"]
    if data["material"] == "cohesive":
        step = step + df["check_cohesive_stress"] * df["check_cohesive_stress nb_rep"]

    df["TTS"] = step / 1000
    df["speedup"] = step[0] / step
    df["mumps"] = df["static_solve"] * df["static_solve nb_rep"]


def plot_measure(ax, df, plotting, label, **kwargs):
    """Plot a given measure."""
    grouped = df.groupby("psize")  # compute stats grouped by number of procs
    med = grouped.median()
    min = grouped.min()
    max = grouped.max()
    min_psize = df["psize"][0]

    print(list(med[plotting]))

    (l,) = ax.plot(med.index, med[plotting], label=f"{label} (median)", **kwargs)

    ax.fill_between(
        med.index, min[plotting], max[plotting], color=l.get_color(), alpha=0.2
    )

    ax.plot(med.index, min_psize * med[plotting][min_psize] / med.index, ls="--", color=l.get_color())
    # ax.boxplot(
    #     data[plotting]["grouped"], positions=psize, widths=[0.1 * s for s in psize]
    # )


plot_measure(
    ax,
    plots["cohesive gcc v5"]["df"],
    plotting,
    "insertion",
    marker="o",
)

plot_measure(
    ax,
    plots["elastic gcc v5"]["df"],
    plotting,
    "no insertion v5",
    marker="o",
)


# Selecting appropriate tick values
psize = np.array(np.unique(plots[list(plots.keys())[0]]["df"]["psize"]))
labels = np.concatenate(
    [[psize[0]], psize[1:][psize[1:] >= 2 * psize[:-1]], [psize[-1]]]
)
# for name, ax in axes.items():
ax.set_xscale("log", base=2)
ax.set_yscale("log")

ylabel = plotting if plotting != "TTS" else "Time to solution"
yunit = "s" if plotting != "speedup" else "-"

ax.set_xlabel("Nb Cores [-]")
ax.set_ylabel(f"""{ylabel} [{yunit}]""")

ax.set_xticks(ticks=labels, labels=map(str, labels))

# Constructing legend with min/max and ideal labels
handles, labels = ax.get_legend_handles_labels()
handles += [
    matplotlib.lines.Line2D([], [], linestyle="--", color="k"),
    matplotlib.patches.Patch(color="k", alpha=0.2),
]
labels += ["ideal", "min/max"]
ax.legend(handles=handles, labels=labels)

fig.tight_layout()
fig.savefig(f"{plotting}.svg", transparent=True, bbox_inches="tight", pad_inches=0.1)
plt.show()
