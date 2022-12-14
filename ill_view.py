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
    print( '\n Illumination/Image file: '+fname )
    imaps = np.load( fname )
    print( ' Shots: 0 - '+str(imaps.shape[0]-2)+', Depth: '+str(imaps.shape[2])+', X dimension: '+str(imaps.shape[1]) )
    print( ' Shot: '+str(imaps.shape[0]-1)+' is an aggregation of the shots illumination maps' )

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

fig, axs = plt.subplots( nrows=1, ncols=1, figsize=(6,4) )
vmin = imaps[ shot ].min() * .15
vmax = imaps[ shot ].max() * .15
print( "Original data range: ", imaps[shot].max(), imaps[shot].min() )
print( "Reduced data range: ",vmin, vmax )
ax = axs.matshow( imaps[ shot ].T, cmap='Greys', aspect=2, vmin=vmin, vmax=vmax )
plt.xlabel( 'X' )
plt.ylabel( 'Depth' )
plt.title( 'Illumination/Image map: '+str(shot) )
plt.grid()
fig.colorbar( ax, ax=axs )
plt.show()

print( '\n Done.' )
