from fancyplot import plot_nbars_nlines

# Figure 1: strong scaling
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
    out="strong-comm-eff.svg")

# Figure 2: weak scaling
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
    out="weak-comm-eff.svg")

# Figure 3: strong scaling on Era
processes = [24, 48, 96, 192, 384, 768, 1536]
bars = {
    "total": [4239, 2332, 1365, 864, 656, 514, 410]
}
lines = {
    "efficiency": [1.0000, 0.9089, 0.7764, 0.6133, 0.4039, 0.2577, 0.1615]
}
plot_nbars_nlines(x=processes, bars=bars, lines=lines, baryscale=2.5,
    title="Strong scaling for large-scale problems (Era)",
    out="strong-on-era.svg")

# Figure 4: strong scaling on Sunway
processes = [8, 16, 32, 64, 128, 256, 512, 1024]
bars = {
    "total": [53492, 27082, 13962, 7322, 3842, 2154, 1330, 956]
}
lines = {
    "efficiency": [1.0000, 0.9876, 0.9578, 0.9132, 0.8702, 0.7761, 0.6284, 0.4371]
}
plot_nbars_nlines(x=processes, bars=bars, lines=lines, baryscale=1.1,
    title="Strong scaling for large-scale problems (Sunway Taihulight)",
    out="strong-on-sunway.svg")