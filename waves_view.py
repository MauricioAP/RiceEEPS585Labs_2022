'''
=
=    Author: Mauricio Araya Polo
=    Date: 09/2020 - present
=
'''

import os, sys
import numpy as np
import matplotlib.pyplot as plt

def input_file_desc( fname ):
    print( '\n Waves file: '+fname )
    waves = np.load( fname )
    print( ' Shots: '+str(waves.shape[0])+', time steps: '+str(waves.shape[1])+', Depth: '+str(waves.shape[3])+', X dimension: '+str(waves.shape[2]) )
    return waves

# CLI processing
if len(sys.argv) < 4:
    print( '\n ERROR!. Only '+str(len(sys.argv)-1)+' parameters read, need 3' )
    print(' Usage: > python3 waves_view.py file.npy shot# step')
    print( ' Example: > python3 waves_view.py file.npy 1 100' )
    if len(sys.argv) >= 2:
        input_file_desc( sys.argv[1] )
    print()
    exit()
shot = int(sys.argv[2])
step = int(sys.argv[3])
# Loading models
waves = input_file_desc( sys.argv[1] )

fig, axs = plt.subplots( nrows=1, ncols=1, figsize=(6,4) )
ax = axs.matshow( waves[shot, step, :, :].T, cmap='Greys', aspect=2 )
plt.xlabel( 'X' )
plt.ylabel( 'Depth' )
plt.title( 'Wavefield: '+str(shot) )
plt.grid()
fig.colorbar( ax, ax=axs )
plt.show()

print( '\n Done.' )
