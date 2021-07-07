from fancyplot import plot_speedup_lines

# Figure 1: strong scaling speedup
processes = [3, 6, 12, 24, 48, 96]
efficiency = [1.0000, 0.9325, 0.8720, 0.7934, 0.4831, 0.2662]

plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for the strong scaling test\n"
          "problem size: 14400 vacancies + 14400 interstitials",
    out="speedup-strong.svg")

# Figure 2: weak scaling speedup
processes = [4, 8, 16, 32, 64, 128]
efficiency = [1.0000, 1.0017, 0.8026, 0.7260, 0.5704, 0.5347]
plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for the weak scaling test\n"
          "problem size (per process): 1000 vacancies + 1000 interstitials",
    out="speedup-weak.svg")

# Figure 3: strong scaling speedup on Era
processes = [24, 48, 96, 192, 384, 768, 1536]
efficiency = [1.0000, 0.9089, 0.7764, 0.6133, 0.4039, 0.2577, 0.1615]
plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for large-scale problems (Era)",
    out="speedup-strong-on-era.svg")

# Figure 4: strong scaling on Sunway
processes = [8, 16, 32, 64, 128, 256, 512, 1024]
efficiency = [1.0000, 0.9876, 0.9578, 0.9132, 0.8702, 0.7761, 0.6284, 0.4371]
plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for large-scale problems (Sunway Taihulight)",
    out="speedup-strong-on-sunway.svg")
