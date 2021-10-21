#@file   holeparticle_correlation.py
#
#@date   Oct 12, 2021
#@author Anton Rikus
#@copyright \n
# This file is part of the program Serenity.\n\n
# Serenity is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.\n\n
# Serenity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.\n\n
# You should have received a copy of the GNU Lesser General
# Public License along with Serenity.
# If not, see <http://www.gnu.org/licenses/>.\n


import numpy as np
import matplotlib.pyplot as plt
import sys
filename = str(sys.argv[1])

correlationplot = np.loadtxt(filename)

fig, ax1 = plt.subplots()

im = ax1.imshow(correlationplot, "Greys", vmin=0, vmax=.8)

ax1.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

ax1.set_xticks(np.arange(len(correlationplot)))
ax1.set_yticks(np.arange(len(correlationplot)))
ax1.set_xticklabels(range(1, 7), fontsize=9)
ax1.set_yticklabels(range(6, 0, -1), fontsize=9)
ax1.set_xlabel("Hole")
ax1.set_ylabel("Particle")
ax1.set_title("Electron correlation")

cbar = ax1.figure.colorbar(im)

fig.tight_layout()
plt.savefig(filename + ".pdf", format="pdf", bbox_inches='tight')
print("Plotted the hole-particle correlation plot to "+ filename + ".pdf")
