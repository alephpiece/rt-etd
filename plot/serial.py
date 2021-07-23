from fancyplot import plot_nbars_nlines

# Figure: performance of the serial code
size = [1024, 2048, 4096, 8192, 16384, 32768, 65536]
bars = {
    "original": [151.77, 329.39, 905.26, 2042.70, 5966.58, 13148.11, 26704.16],
    "optimized(v1)": [129.50, 260.18, 612.49, 1504.53, 3666.27, 8922.72, 20350.83],
    "optimized(v2)": [29.06, 59.43, 271.40, 843.52, 2759.06, 6587.04, 14397.25],
}
lines = {
    "time reduced(v2)": [0.8085, 0.8196, 0.7002, 0.5871, 0.5376, 0.4990, 0.4609]
}

plot_nbars_nlines(x=size, bars=bars, lines=lines,
    baryscale=2, lineylimit=0.90, bardatalabel=False,
    linefmt="k^-",
    xlabel="problem size",
    title="Performance of optimized code\n"
          "simulation time: $10^5$s",
    legendloc="center left",
    style="seaborn-deep",
    out="serial-optimization.svg")