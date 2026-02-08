from fancyplot import plot_nbars_nlines


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

# Figure 3: strong scaling speedup on Era
processes = [24, 48, 96, 192, 384, 768, 1536]
efficiency = [1.0000, 0.9089, 0.7764, 0.6133, 0.4039, 0.2577, 0.1615]
plot_speedup_lines(x=processes, efficiency=efficiency,
    title="Speedup for large-scale problems (Era)",
    out="speedup-strong-on-era.svg")

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

