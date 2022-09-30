'''
=
=    Author: Mauricio Araya Polo
=    Date: 09/2020 - present
=
'''

import os, sys
from utils import *

# CLI processing
if len(sys.argv) < 2:
    print( '\n ERROR!. Only '+str(len(sys.argv)-1)+' parameters read, need 1' )
    print(' Usage: > python3 modeling2D.py synthmodel2D or SEGsaltmodel')
    print( ' Example: > python3 modeling2D.py synthmodel2D' )
    print()
    exit()

# Loading models
print( "\nLoading "+sys.argv[1]+".\n" )
model = load_model( sys.argv[1], "npy", "." )

# model dimensions (nz is depth and nx is the horizontal axis)
(nz, nx) = model.vel.shape
print( "Model size: z %d, x %d\n" % (nz, nx) )

# peak source frequency
freq = 18.0

# Setting number of time steps
nt = 2000

# 15% padding from edges, for shot and receivers location
start = int(0.15*nx)
stop = int(nx - 0.15*nx)

multishot = True
# ns number of shots
if multishot:
    ns = 10
    # shots locations
    sloc = np.linspace( start, stop, ns ).astype(np.int32)
else:
    # just one shot in the middle of the surface
    ns = 1
    sloc = [ int(nx/2) ]

# number of receivers and locations
nr = 100
rloc = np.linspace( start, stop, nr ).astype(np.int32)

# data structure with experimental septup
case = caseStruct( model.name, freq, nt, ns, nr, sloc, rloc, 4, 10 )

print("Experiment main parameters:\n", case)
# running simulation
traces, waves = modeling( model, case )

# saving unsorted  seismic data to disk
print('Saving seismic data (traces_'+sys.argv[1]+'.npy) to disk', traces.shape)
np.save( 'traces_'+sys.argv[1]+'.npy', traces )
print('Saving wavefields (waves_'+sys.argv[1]+'.npy) to disk', waves.shape)
np.save( 'waves_'+sys.argv[1]+'.npy', waves )
