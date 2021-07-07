from fancyplot import plot_nbars_nlines

# Figure: performance of the serial code
size = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576]
bars = {
    "original": [454, 479, 532, 646, 853, 1282, 2129, 4288, 10019, 22884, 50116],
    "optimized(v1)": [154, 166, 195, 258, 385, 630, 1108, 3229, 8429, 18574, 37558],
    "optimized(v2)": [263, 277, 308, 375, 497, 759, 1263, 2237, 6263, 15465, 34095]
}
lines = {
    "time reduced(v2)": [0.4213, 0.4218, 0.4215, 0.4202, 0.4175, 0.4079, 0.4070, 0.4783, 0.3749, 0.3242, 0.3197]
}
plot_nbars_nlines(x=size, bars=bars, lines=lines,
    baryscale=2, lineylimit=0.55, bardatalabel=False,
    linefmt="k^-",
    xlabel="problem size",
    title="Performance of optimized code\n"
          "simulation time: $10^5$s",
    legendloc="center left",
    style="seaborn-deep",
    out="perf-serial.svg")
