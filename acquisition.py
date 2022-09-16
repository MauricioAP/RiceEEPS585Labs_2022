'''
=
=    Author: Mauricio Araya Polo
=    Date: 09/2021 - present
=
'''

import numpy as np
import matplotlib.pyplot as plt
import sys

assert len(sys.argv) == 7, "Missing params: X[grid] Y[grid] dx[m] dy[m] padx[m] pady[m]"

# meters
dx = float(sys.argv[3])
dy = float(sys.argv[4])
sdx = float(sys.argv[5])
sdy = float(sys.argv[6])
# grid points
nx = int(sys.argv[1])
ny = int(sys.argv[2])
padx =int(sdx/2) 
pady =int(sdy/2) 

nshotsx = int(nx*dx/sdx)
nshotsy = int(ny*dy/sdy)
print( 'Shots on X', nshotsx )
print( 'Shots on Y', nshotsy )
nshotsx = nshotsx - 1
nshotsy = nshotsy - 1
shotsx = []
shotsx.append( padx )
current = padx
for i in range( 1, nshotsx+1 ):
    current += sdx
    shotsx.append(current)
shotsy = []
shotsy.append( pady )
current = pady
for i in range( 1, nshotsy+1 ):
    current += sdy
    shotsy.append(current)
print( 'Shots locations X and Y', shotsx, shotsy )
shots = []
for i in range( len(shotsx) ):
    for j in range( len(shotsy)  ):
        shots.append( shotsx[i] )
        shots.append( shotsy[j] )
        plt.plot( shotsx[i],shotsy[j], marker='o', markerfacecolor='blue', markeredgecolor='blue' )
print( 'All shots location', int(len(shots)/2), shots )
plt.grid()
plt.ylabel("Y")
plt.xlabel("X")
plt.show()

print( "Saving acquisition plan")
np.save("acquisition_plan.npy", np.array(shotsx)/dx )
