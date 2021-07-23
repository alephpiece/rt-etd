from fancyplot import plot_nbars_nlines, plot_speedup_lines

# Figure 1: strong scaling efficiency (on Sunway)
processes = [64, 128, 256, 512, 1024, 2048, 4096]

bars = {
    "time": [40434.10, 21715.28, 10710.07, 5245.50, 2717.61, 2062.84, 1887.04]

}
lines = {
    "efficiency": [1.0000, 0.9310, 0.9438, 0.9635, 0.9299, 0.6125, 0.3348]

}
plot_nbars_nlines(x=processes, bars=bars, lines=lines, baryscale=1.1,
    title="Strong scaling for large-scale problems (Sunway Taihulight)\n"
          "problem size: $2^{23}$, simulation time: $10^5$",
    out="strong-on-sunway.svg")

# Figure 2: weak scaling efficiency (on Sunway)
processes = [64, 128, 256, 512, 1024, 2048, 4096]
bars = {
    "time": [1915.92, 2060.68, 2144.15, 2334.73, 2652.48, 3134.78, 3785.52]
}
lines= {
    "efficiency": [1.0000, 0.9298, 0.8936, 0.8206, 0.7223, 0.6112, 0.5061]

}
plot_nbars_nlines(x=processes, bars=bars, lines=lines, baryscale=2.5,
    title="Weak scaling for large-scale problems (Sunway Taihulight)\n"
          "problem size (per process): 8192, simulation time: $10^5$",
    out="weak-on-sunway.svg")

# Figure 3: strong scaling speedup (on Sunway)
processes = [64, 128, 256, 512, 1024, 2048, 4096]
efficiency = [1.0000, 0.9310, 0.9438, 0.9635, 0.9299, 0.6125, 0.3348]

plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for the strong scaling test (Sunway Taihulight)\n"
          "problem size: $2^{23}$, simulation time: $10^5$",
    out="speedup-strong-on-sunway.svg")

# Figure 4: weak scaling speedup (on Sunway)
processes = [64, 128, 256, 512, 1024, 2048, 4096]
efficiency = [1.0000, 0.9298, 0.8936, 0.8206, 0.7223, 0.6112, 0.5061]
plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for the weak scaling test (Sunway Taihulight)\n"
          "problem size (per process): 8192, simulation time: $10^5$",
    out="speedup-weak-on-sunway.svg")