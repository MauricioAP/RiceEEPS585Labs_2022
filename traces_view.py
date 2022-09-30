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
    traces = np.load( fname )
    print( ' Shots: '+str(traces.shape[0])+', receivers: '+str(traces.shape[1])+', timesteps: '+str(traces.shape[2]) )
    return traces

# CLI processing
if len(sys.argv) < 5:
    print( '\n ERROR!. Only '+str(len(sys.argv)-1)+' parameters read, need 4' )
    print( ' Usage: > python3 traces_view.py file.npy shot# time_start time_end' )
    print( ' Example: > python3 traces_view.py file.npy 1 100 500' )
    if len(sys.argv) >= 2:
        input_file_desc( sys.argv[1] )
    print()
    exit()
shotS = int(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])
if start > end:
    print( '\n ERROR!. time_end has to be larger than time_start' )
    print( ' Usage: > python3 traces_view.py file.npy shot# time_start time_end' )
    print()
    exit()

# Loading models
traces = input_file_desc( sys.argv[1] )

#fig, axs = plt.subplots( nrows=1, ncols=traces.shape[0], figsize=(6,4) )
#if traces.shape[0] > 1:
#    for shot in range(traces.shape[0]):
#        ax = axs[ shot ]
#        tmp = traces[shot, :, :].flatten()
#        pcm = ax.hist(tmp, bins=10 )
#        ax.set_title('Shot '+str(shot))
#        if shot > 0:
#            axs[shot].set_yticks([])
#    axs[0].set_ylabel('Frequency')
#    axs[0].set_xlabel('Amplitude range')
#else:
#    tmp = traces[0, :, :].flatten()
#    axs.hist(tmp, bins=10 )
#    axs.set_title('Shot '+str(shotS))
#    axs.set_ylabel('Frequency')
#    axs.set_xlabel('Amplitude range')
#plt.show()

fig, axs = plt.subplots( nrows=1, ncols=traces.shape[0], figsize=(7,4) )
if traces.shape[0] > 1:
    for shot in range(traces.shape[0]):
        ax = axs[ shot ]
        pcm = ax.matshow(traces[shot, :, :].T, cmap='Greys', vmin=-3, vmax=3, aspect='auto' )
        ax.set_xticks([ 1, 50, 100 ])
        ax.set_title('Shot '+str(shot))
        if shot > 0:
            axs[shot].set_yticks([])
    fig.colorbar( pcm, ax=axs )
    axs[0].set_xlabel('X')
    axs[0].set_ylabel('Time')
else:
    pcm = axs.matshow(traces[0, :, :].T, cmap='Greys', vmin=-3, vmax=3, aspect='auto' )
    axs.set_xticks([ 1, 50, 100 ])
    axs.set_title('Shot '+str(shotS))
    fig.colorbar( pcm, ax=axs )
    axs.set_xlabel('X')
    axs.set_ylabel('Time')
plt.show()

fig, axs = plt.subplots( nrows=1, ncols=1, figsize=(4,4) )
window = start + ( traces.shape[2] - end)
axs.plot( traces[shotS,50,start:end], np.linspace(1, traces.shape[2]-window, traces.shape[2]-window) )
plt.xlabel('Amplitude')
plt.ylabel('Time')
axs.invert_yaxis()
plt.title("Shot "+str(shotS)+", receiver 50")
plt.show()

print('\n Done.')
