me0 = "EnergyLandscape"

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import NullLocator
from Utils import fs, set_mplrc
set_mplrc(fs)
mpl.rc("lines", linewidth=4)

fig, ax = plt.subplots(1,1, figsize=fs["figsize"])

## xy
ax.arrow(-0.1, 0.0, 1.1, 0.0, head_width=0.03, head_length=0.07, length_includes_head=True, overhang=0.5, fc='k', ec='k')
ax.arrow(0.0, -0.1, 0.0, 1.1, head_width=0.03, head_length=0.07, length_includes_head=True, overhang=0.5, fc='k', ec='k')
ax.text(-0.08, 0.95, r"$E$", fontsize=fs["fsa"])

## Energy levels
EA1 = 0.15
EBC = 0.65
EA2 = 0.05
eab = EBC+0.10
ebc = EBC+0.07
## Width of energy levels and barriers
w = 0.17
wb = 0.08

## Unprimed
ax.plot([0*w+0*wb,1*w+0*wb],[EA1,EA1], "b-")
ax.plot([1*w+1*wb,2*w+1*wb],[EBC,EBC], "b-")
ax.plot([2*w+2*wb,3*w+2*wb],[EBC,EBC], "b-")
ax.plot([3*w+3*wb,4*w+3*wb],[EA2,EA2], "b-")
## Primed
ax.plot([0*w+0*wb,1*w+0*wb],[EA2,EA2], "r-")
ax.plot([1*w+1*wb,2*w+1*wb],[EBC-0.01,EBC-0.01], "r-")
ax.plot([2*w+2*wb,3*w+2*wb],[EBC-0.01,EBC-0.01], "r-")
ax.plot([3*w+3*wb,4*w+3*wb],[EA1,EA1], "r-")

## Barriers
ax.plot([1*w+0*wb,1*w+1*wb],[eab+EA1-EA2,eab+EA1-EA2], "b:")
ax.plot([2*w+1*wb,2*w+2*wb],[ebc,ebc], "b:")
ax.plot([3*w+2*wb,3*w+3*wb],[eab,eab], "b:")
ax.plot([1*w+0*wb,1*w+1*wb],[eab,eab], "r:")
ax.plot([2*w+1*wb,2*w+2*wb],[ebc-0.01,ebc-0.01], "r:")
ax.plot([3*w+2*wb,3*w+3*wb],[eab+EA1-EA2,eab+EA1-EA2], "r:")

## Scale labels
## TlnD
ax.arrow(w+0.02, EA2, 0.0, EA1-EA2, length_includes_head=True, overhang=0.0, fc='y', ec='y')
ax.arrow(w+0.02, EA1, 0.0, EA2-EA1, length_includes_head=True, overhang=0.0, fc='y', ec='y')
ax.text(w+0.04, 0.4*(EA1+EA2), r"$T\ln \Delta$", fontsize=fs["fsa"])
#
ax.arrow(w-0.02, eab+EA1-EA2, 0.0, EA2-EA1, length_includes_head=True, overhang=0.0, fc='y', ec='y')
ax.arrow(w-0.02, eab, 0.0, EA1-EA2, length_includes_head=True, overhang=0.0, fc='y', ec='y')
# ax.text(w-0.04, eab+0.4*(EA1-EA2), r"$T\ln \Delta$", ha="right", fontsize=fs["fsa"])
ax.arrow(3*w+3*wb+0.02, eab+EA1-EA2, 0.0, EA2-EA1, length_includes_head=True, overhang=0.0, fc='y', ec='y')
ax.arrow(3*w+3*wb+0.02, eab, 0.0, EA1-EA2, length_includes_head=True, overhang=0.0, fc='y', ec='y')
#
ax.arrow(3*w+3*wb-0.02, EA2, 0.0, EA1-EA2, length_includes_head=True, overhang=0.0, fc='y', ec='y')
ax.arrow(3*w+3*wb-0.02, EA1, 0.0, EA2-EA1, length_includes_head=True, overhang=0.0, fc='y', ec='y')

## Axis limits
ax.set_xlim(-0.1,1.0)
ax.set_ylim(-0.1,1.0)

## Axis labels
off1,off2 = 0.1,0.6
ax.text(w*off1, -0.08, r"$\mathcal A_1$", fontsize=fs["fsa"], color="b")
ax.text(w*off2, -0.087, r"$\mathcal A^\prime_1$", fontsize=fs["fsa"], color="r")
ax.text(w*off1+w+wb, -0.08, r"$\mathcal B$", fontsize=fs["fsa"], color="b")
ax.text(w*off2+w+wb, -0.08, r"$\mathcal B^\prime$", fontsize=fs["fsa"], color="r")
ax.text(w*off1+2*w+2*wb, -0.08, r"$\mathcal C$", fontsize=fs["fsa"], color="b")
ax.text(w*off2+2*w+2*wb, -0.08, r"$\mathcal C^\prime$", fontsize=fs["fsa"], color="r")
ax.text(w*off1+3*w+3*wb, -0.08, r"$\mathcal A_2$", fontsize=fs["fsa"], color="b")
ax.text(w*off2+3*w+3*wb, -0.087, r"$\mathcal A^\prime_2$", fontsize=fs["fsa"], color="r")

ax.yaxis.set_major_locator(NullLocator())

ax.axis("off")
ax.patch.set_visible(False)
plt.gcf().patch.set_visible(False)

plt.show()