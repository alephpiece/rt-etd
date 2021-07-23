#!/usr/bin/env python3

import math
import numpy as np  # multi-dimentional arrays and arithmetics
from matplotlib import pyplot, ticker   # library for making plots, similar to MatLab

def _plot_bars(axes, x, data, barwidth, baryscale, xlabel, lylabel, datalabel):
    """An auxiliary function for drawing bars.
    Parameters
    ==========
    x: x coordinates
    data: dictionary of name-list pairs
    barwidth: width of each bar
    baryscale: scale of the y-axis
    xlabel: label of the x-axis, e.g. "processes"
    lylabel: label of the left y-axis, e.g. "execution time/s"
    datalabel: boolean to indicate whether to show data labels or not
    """
    # X coordinates of center points
    xn = np.arange(1, len(x)+1)
    # Offset in the x direction
    offset = (len(data) - 1) * barwidth / 2

    # A counter used to shift each bar to the right position
    count = 0
    for yname, yvalue in data.items():
        # X coordinates of the current bar
        x1 = xn - offset + count * barwidth

        # Plotting the bar
        axes.bar(x=x1, height=yvalue, width=barwidth, label=yname)

        # Add data labels
        if datalabel:
            for xx, yy in zip(x1, yvalue):
                axes.text(xx, yy + 10, str(yy), ha="center", va="bottom", fontsize="x-small")

        count += 1

    # Set other attributes for bars
    # X-axis
    axes.set_xticks(xn)      # ticks of the x-axis, e.g. "1, 2, 3, ..."
    axes.set_xlabel(xlabel)  # labels along the x-axis, e.g. "8, 16, 32, ..."

    # Limits of axes
    barxlimit = max(xn + offset) + 0.5  # maximum x limit of bars
    barylimit = 0   # maximum y limit of bars
    for yvalue in data.values():
        barylimit = max(barylimit, max(yvalue))
    barylimit = barylimit * baryscale
    axes.set_xlim([0.5, barxlimit])

    # Format x tick labels as LaTeX formula
    formatted_x = list(map(lambda x: "$2^{{{0}}}$".format(int(math.log(x, 2))),
                            x)
                        )
    axes.set_xticklabels(formatted_x)
    #axes.set_xticklabels(x)

    # Setup the y-axis
    axes.set_ylabel(lylabel)
    axes.set_ylim([0, barylimit])
    axes.ticklabel_format(axis="y", style="scientific", scilimits=(0,0), useMathText=True)

    # Disable the background grid
    axes.grid(False)


def _plot_stacks(axes, x, data, barwidth, baryscale, xlabel, lylabel, datalabel):
    """An auxiliary function for drawing stacks.
    Parameters
    ==========
    x: x coordinates
    data: dictionary of name-list pairs
    barwidth: width of each bar
    baryscale: scale of the y-axis
    xlabel: label of the x-axis, e.g. "processes"
    lylabel: label of the left y-axis, e.g. "execution time/s"
    datalabel: boolean to indicate whether to show data labels or not
    """
    # X coordinates of center points
    xn = np.arange(1, len(x)+1)

    ybottom = [0] * len(x)
    for yname, yvalue in data.items():
        # Plotting the stack
        axes.bar(xn, height=yvalue, width=barwidth, bottom=ybottom, label=yname)

        # Add data labels
        if datalabel:
            for xx, yy in zip(xn, yvalue):
                axes.text(xx, yy + 10, str(yy), ha="center", va="bottom", fontsize="x-small")

        # Bottom of the current stack
        ybottom = [ybottom[i] + yvalue[i] for i in range(len(yvalue))]

    # Set other attributes for bars
    # X-axis
    barxlimit = max(xn) + 0.5  # maximum x limit of bars
    axes.set_xlim([0.5, barxlimit])
    axes.set_xticks(xn)      # ticks
    axes.set_xlabel(xlabel)

    # Format x tick labels as LaTeX formula
    formatted_x = list(map(lambda x: "$2^{{{0}}}$".format(int(math.log(x, 2))),
                            x)
                        )
    axes.set_xticklabels(formatted_x)
    #axes.set_xticklabels(x)

    # Setup the y-axis
    barylimit = max(ybottom) * baryscale   # maximum y limit of bars
    axes.set_ylim([0, barylimit])
    axes.set_ylabel(lylabel)
    axes.ticklabel_format(axis="y", style="scientific", scilimits=(0,0), useMathText=True)

    # Disable the background grid
    axes.grid(False)


def _plot_lines(axes, x, data, linefmt, lineylimit, rylabel):
    """An auxiliary function for drawing lines.
    Parameters
    ==========
    x: x coordinates
    data: dictionary of name-list pairs
    linefmt: format of lines
    lineylimit: limit of the y-axis
    rylabel: label of the right y-axis, e.g. "efficiency"
    """
    # X coordinates of center points
    xn = np.arange(1, len(x)+1)
    # Offsets used to adjust data labels
    xtextnudge = -0.03
    ytextnudge = 0.01 * lineylimit

    for yname, yvalue in data.items():
        axes.plot(xn, yvalue, linefmt, label=yname)
    
        # Add data labels
        for xx, yy in zip(xn, yvalue):
            axes.text(xx + xtextnudge, yy + ytextnudge, f"{yy:.1%}", ha="left", va="bottom", fontsize="x-small")
    
    # Other attributes
    axes.set_ylabel(rylabel)
    axes.set_ylim([0.0, lineylimit])
    axes.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=1))

    # Disable the background grid
    axes.grid(False)

def _locate_legend(loc="upper right"):
    """An auxiliary function to get the relative location of the legend.
    Parameters
    ==========
    loc: a string indicating the position.
        "upper right", "center left", ...

    Return
    ======
    loc, bbox
    """
    if loc == "upper right":
        bbox = (1, 1)
    elif loc == "upper left":
        bbox = (0, 1)
    elif loc == "lower right":
        bbox = (1, 0)
    elif loc == "center right":
        bbox = (1, 0.5)
    elif loc == "center left":
        bbox = (0, 0.5)

    return loc, bbox


def plot_nbars_nlines(x, bars, lines,
                      barwidth=0.31, baryscale=1.1, bardatalabel=True,
                      linefmt="^-", lineylimit=1.05,
                      xlabel="processes", lylabel="execution time/s", rylabel="percentage",
                      title="Scalability",
                      legendloc="upper right",
                      style="classic",
                      stacked=False,
                      out="figure.svg"):
    """Draw a figure with N bars and N lines.
    Parameters
    ==========
    x: list of x-coordinates
    bars: dictionary of "name-list" pairs
        {"total": [3, 4, 5],
         "communication":  [1, 1, 1]
         }
    lines: dictionary of "name-list" pairs
        {"efficiency": [0.5, 0.5, 0.5]
        }
    title: figure title
    out: output filename
    barwidth: width of each bar
    baryscale: scale of the y-axis of bars
    """
    # Style
    pyplot.style.use(style)

    # Create plot object
    fig, ax1 = pyplot.subplots()

    ## First plot: bars
    if stacked:
        bars_plotter = _plot_stacks
    else:
        bars_plotter = _plot_bars
    bars_plotter(axes=ax1, x=x, data=bars, barwidth=barwidth, baryscale=baryscale, xlabel=xlabel, lylabel=lylabel, datalabel=bardatalabel)

    ## Second plot: lines
    if lines:
        ax2 = ax1.twinx()
        _plot_lines(axes=ax2, x=x, data=lines, linefmt=linefmt, lineylimit=lineylimit, rylabel=rylabel)

    ## Common attributes
    fig.suptitle(title)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    loc, bbox = _locate_legend(legendloc)
    fig.legend(loc=loc, bbox_to_anchor=bbox, bbox_transform=ax1.transAxes, fontsize="small")

    ## Output the figure
    #pyplot.show()
    pyplot.savefig(out)


def plot_speedup_lines(x, efficiency,
                      linefmt="^-", 
                      xlabel="processes", ylabel="speedup",
                      title="Speedup",
                      legendloc="lower right",
                      style="classic",
                      out="figure.svg"):
    """Draw lines.
    Parameters
    ==========
    x: x coordinates
    efficiency: list, parallel efficiency
    """
    # Style
    pyplot.style.use(style)

    # Create plot object
    fig, ax = pyplot.subplots()

    # Compute speedup
    idealspeedup = [i/x[0] for i in x]
    actualspeedup = [i*j for i,j in zip(idealspeedup, efficiency)]

    # X coordinates of center points
    xn = np.arange(1, len(x)+1)
    ytextnudge = 1.05

    # Draw the ideal speedup
    yvalue = idealspeedup
    ax.plot(xn, yvalue, linefmt, label="ideal speedup")
    # Add data labels
    for xx, yy in zip(xn, yvalue):
        ax.text(xx, yy * ytextnudge, f"{yy:.0f}", ha="right", va="bottom", fontsize="x-small")

    # Draw the actual speedup
    yvalue = actualspeedup
    ax.plot(xn, yvalue, linefmt, label="actual speedup")

    # Add data labels
    for xx, yy in zip(xn, yvalue):
        ax.text(xx, yy / ytextnudge, f"{yy:.1f}", ha="left", va="top", fontsize="x-small")
    
    # Other attributes
    ax.set_xticks(xn)
    ax.set_xlabel(xlabel)
    ax.set_xticklabels(x)
    ax.set_yscale("log", basey=2)
    ax.set_ylabel(ylabel)
    ax.set_xlim([0.5, len(xn)+0.5])
    ax.set_ylim([0.75, 1.5*max(idealspeedup)])
    ax.grid(False)

    ## Common attributes
    fig.suptitle(title)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    loc, bbox = _locate_legend(legendloc)
    fig.legend(loc=loc, bbox_to_anchor=bbox, bbox_transform=ax.transAxes, fontsize="small")

    ## Output the figure
    #pyplot.show()
    pyplot.savefig(out)
