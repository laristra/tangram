# Simple script to read and plot the material polygons
import numpy, sys
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

# TODO - error checking on args
fname = sys.argv[1]
# THIS IS FRAGILE
# the first two lines are diagnostics
data = numpy.loadtxt(fname, skiprows=2)

# to color the polygons
colors = ['y', 'b', 'r']

# each entry is a material polygon
# since we are plotting in 2d, we are ignoring the third dimension
# this means that we should also only have four nodes
polys = []
for matpoly in data:
    matid = int(matpoly[0])
    # do this slowly and by hand for now
    xy = []
    xy.append([matpoly[1], matpoly[2]])
    xy.append([matpoly[4], matpoly[5]])
    xy.append([matpoly[10], matpoly[11]])
    xy.append([matpoly[7], matpoly[8]])
    polys.append(Polygon(xy, facecolor=colors[matid]))
p = PatchCollection(polys, match_original=True)

fig, ax = plt.subplots()
ax.add_collection(p)
ax.set_xticks(numpy.linspace(0, 1, 21), minor=True)
ax.set_yticks(numpy.linspace(0, 1, 21), minor=True)
plt.grid(which='minor', color='k', lw=1)
ax.set_aspect('equal')
plt.savefig("%s.png" % fname.split("_")[-1], bbox_inches="tight")
