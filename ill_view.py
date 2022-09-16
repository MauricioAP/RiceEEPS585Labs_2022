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
    print( '\n Trace file: '+fname )
    imaps = np.load( fname )
    print( ' Shots: 0 - '+str(imaps.shape[0]-1)+', Depth: '+str(imaps.shape[2])+', X dimension: '+str(imaps.shape[1]) )
    print( ' Shot: '+str(imaps.shape[0])+' is an aggregation of the shots illumination maps' )

    return imaps

# CLI processing
if len(sys.argv) < 3:
    print( '\n ERROR!. Only '+str(len(sys.argv)-1)+' parameters read, need 2' )
    print(' Usage: > python3 ill_view.py file.npy shot#')
    print( ' Example: > python3 ill_view.py file.npy 1' )
    if len(sys.argv) >= 2:
        input_file_desc( sys.argv[1] )
    print()
    exit()
shot = int(sys.argv[2])
# Loading models
imaps = input_file_desc( sys.argv[1] )
tmp = np.zeros( (imaps.shape[1], imaps.shape[2]), dtype=np.float32)
if shot == imaps.shape[0]:
    print( "Showing aggregated illumination map" )
    for i in range(imaps.shape[0]):
        tmp += imaps[i,:,:]
    vmin = tmp.min() * .15
    vmax = tmp.max() * .15
    print( "Original data range: ", tmp.max(), tmp.min() )
else:
    vmin = imaps[ shot ].min() * .15
    vmax = imaps[ shot ].max() * .15
    print( "Original data range: ", imaps[shot].max(), imaps[shot].min() )

fig, axs = plt.subplots( nrows=1, ncols=1, figsize=(6,4) )
print( "Reduced data range: ",vmin, vmax )
if shot == imaps.shape[0]:
    ax = axs.matshow( tmp.T, cmap='Greys', aspect=2, vmin=vmin, vmax=vmax )
else:
    ax = axs.matshow( imaps[ shot ].T, cmap='Greys', aspect=2, vmin=vmin, vmax=vmax )
plt.xlabel( 'X' )
plt.ylabel( 'Depth' )
plt.title( 'Illumination map: '+str(shot) )
plt.grid()
fig.colorbar( ax, ax=axs )
plt.show()

print( '\n Done.' )
