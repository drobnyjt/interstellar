import numpy as np
import matplotlib.pyplot as plt
from figurefirst import FigureLayout

layout = FigureLayout('bubbles.svg')
layout.make_mplfigures()
ax = layout.axes['axis_name']['axis']
data = np.genfromtxt('RANGE_3D_H_200_05_c.txt', skip_header=18)
x = data[:,1]
mean = np.mean(x)
std = np.std(x)
num_bins = 50
bins = np.linspace(mean-6.*std, mean+6.*std, num_bins)
values, bins = np.histogram(x, bins=bins, density=True)
centers = bins[:-1] + (bins[1:] - bins[:-1])/2.
ax.plot(centers, values, linewidth=1, color='black')
ax.fill_between(centers, values, alpha=0.25, color='black')
ax.axis('off')
ax.set_xticks([])
ax.set_yticks([])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
layout.insert_figures('target_layer_name')
layout.write_svg('bubbles.svg')
