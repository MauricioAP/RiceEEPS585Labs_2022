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
    model = np.load( fname )
    print( ' Depth: '+str(model.shape[0])+', X dimension: '+str(model.shape[1]) )
    return model

# CLI processing
if len(sys.argv) < 2:
    print( '\n ERROR!. Only '+str(len(sys.argv)-1)+' parameters read, need 1' )
    print( ' Usage: > python3 model_view.py synthmodel2D.npy or SEGsaltmodel.npy' )
    print( ' Example: > python3 model_view.py synthmodel2D.npy' )
    print()
    exit()

# Loading models
model = input_file_desc( sys.argv[1] )

fig, axs = plt.subplots( nrows=1, ncols=1, figsize=(4,4) )
ax = axs.matshow( model.T, cmap='seismic' )
plt.xlabel('X')
plt.ylabel('Depth')
fig.colorbar( ax, ax=axs )
plt.show()

print('\n Done.')
