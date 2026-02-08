from fancyplot import plot_nbars_nlines, plot_speedup_lines

# Figure 9a: strong scaling
processes = [3, 6, 12, 24, 48, 96]
bars = {
    "total": [81046, 43457, 23236, 12769, 10485, 9515],
    "communication": [144, 176, 200, 1840, 3048, 4280]
}
lines = {
    "efficiency": [1.0000, 0.9325, 0.8720, 0.7934, 0.4831, 0.2662]
}
plot_nbars_nlines(x=processes, bars=bars, lines=lines,
    title="Comparison between total and communication time (strong scaling)\n"
          "problem size: 14400 vacancies + 14400 interstitials",
    out="9a-strong-IV.svg")

# Figure 9b: strong scaling speedup
processes = [3, 6, 12, 24, 48, 96]
efficiency = [1.0000, 0.9325, 0.8720, 0.7934, 0.4831, 0.2662]

plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for the strong scaling test\n"
          "problem size: 14400 vacancies + 14400 interstitials",
    out="9b-speedup-strong-IV.svg")

# Figure 10a: weak scaling
processes = [4, 8, 16, 32, 64, 128]
bars = {
    "total": [19057, 19025, 23744, 26248, 33407, 35638],
    "communication": [232, 256, 424, 2536, 3176, 5336]
}
lines= {
    "efficiency": [1.0000, 1.0017, 0.8026, 0.7260, 0.5704, 0.5347]
}
plot_nbars_nlines(x=processes, bars=bars, lines=lines, baryscale=2.5,
    title="Comparison between total and communication time (weak scaling)\n"
          "problem size (per process): 1000 vacancies + 1000 interstitials",
    out="10a-weak-IV.svg")

# Figure 10b: weak scaling speedup
processes = [4, 8, 16, 32, 64, 128]
efficiency = [1.0000, 1.0017, 0.8026, 0.7260, 0.5704, 0.5347]
plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for the weak scaling test\n"
          "problem size (per process): 1000 vacancies + 1000 interstitials",
    out="10b-speedup-weak-IV.svg")
