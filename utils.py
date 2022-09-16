'''
=
=    Author: Mauricio Araya Polo
=    Date: 09/2020 - present
=
'''

from typing import NamedTuple
import numpy as np
import matplotlib.pyplot as plt
import ctypes
from numpy.ctypeslib import ndpointer

bcoeffs = [ 0.0, -205.0/72.0, 8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.10 ]

class modelStruct(NamedTuple):
  name: str
  path: str
  ftype: str
  vel: float
  dx: float
  dz: float
  nz: int
  nx: int

class caseStruct(NamedTuple):
  name: str
  freq: float
  nt: int
  ns: int
  nr: int
  sloc: int
  rloc: int
  stencil: int
  skip: int

def load_model( name, ftype, path ):
  vel = np.load( path+"/"+name+"."+ftype )
  vel = np.transpose( vel )
  vel = np.ascontiguousarray(vel, dtype=np.float32)
  model = modelStruct( name, path, ftype, vel, 10.0, 10.0, vel.shape[0], vel.shape[1] )

  return model

def check_num( model, dt, freq ):
  cond = False

  # CFL check
  v_max = model.vel.max()
  dt_max = min(model.dx, model.dz) * np.sqrt(3/8) / v_max
  if dt_max < dt:
      print("WARNING: chosen dt violates CFL condition (dt_max = %f). The simulation can be unstable.\n" % dt_max)
      cond = True

  # Wavelength check
  v_min = model.vel.min()
  np_min = v_min / freq / max(model.dx, model.dz)
  if np_min < 8:
      print("WARNING: only %f grid points per wavelength in at least one dimension. This can lead to dispersion.\n" % np_min)
      cond = True

  return cond

def compute_coeffs( dz, dx, coeffs ):
  dxi2 = 1.0 / dx
  dxi2 *= dxi2
  b0x = bcoeffs[1] * dxi2
  coeffs[1] = bcoeffs[2] * dxi2
  coeffs[2] = bcoeffs[3] * dxi2
  coeffs[3] = bcoeffs[4] * dxi2
  coeffs[4] = bcoeffs[5] * dxi2

  dzi2 = 1.0 / dz
  dzi2 *= dzi2
  b0z = bcoeffs[1] * dzi2
  coeffs[5] = bcoeffs[2] * dzi2
  coeffs[6] = bcoeffs[3] * dzi2
  coeffs[7] = bcoeffs[4] * dzi2
  coeffs[8] = bcoeffs[5] * dzi2

  coeffs[0] = b0x + b0z

def prep( case, model ):
  coeffs = np.zeros( (16), dtype=np.float32 )
  params = np.zeros( (14), dtype=np.int32 )

  vmax = model.vel.max()
  print( "Vel max: %f" % vmax )
  dt = 8.0 * (0.09 * model.dz/vmax/np.sqrt(2))
  print( "dt: %f" % dt )

  coeffs[ 13 ] = model.dx
  coeffs[ 14 ] = model.dz
  coeffs[ 15 ] = dt
  if check_num( model, dt, case.freq ) == True:
      exit()
  print("Numerical check passed")

  compute_coeffs( model.dz, model.dx, coeffs )

  stencil = case.stencil
  nz = model.nz
  nx = model.nx
  #params[0] = nx + 2*stencil
  params[0] = stencil
  params[1] = nz
  params[2] = nx
  col = nz + 2*stencil
  plane = col * params[0]
  params[4] = plane
  params[5] = 2*plane
  params[6] = 3*plane
  params[7] = 4*plane
  params[8] = nz * nx
  params[9] = col
  params[10] = 2*col
  params[11] = 3*col
  params[12] = 4*col
  params[13] = stencil + stencil*nz

  return params, coeffs

def make_wavelet( case, dt ):
  #wavelet = np.zeros( (case.nt), dtype=np.float32 )

  # Ricker wavelet
  t_peak = 1.0 / case.freq
  t = np.linspace(0, (case.nt - 1) * dt, case.nt) - t_peak
  t = t**2 * case.freq**2 * np.pi**2
  return (1 - 2*t) * np.exp(-t)

def shot_injection_transp( field, vel, sample, shotloc, dt, stencil ):
  if sample != 0.0:
    imin = max(stencil, shotloc - stencil)
    imax = min(shotloc + stencil, vel.shape[0] + stencil)
    #print( imin, imax)
    jmin = stencil
    jmax = 2 * stencil
    for i in range(imin, imax):
      for j in range(jmin, jmax):
        dist = (j - jmin)**2 + (i - imin)**2
        dist = np.exp(-dist)
        tmp = vel[i - stencil, j - stencil] * dt
        tmp2 = sample * dist * tmp * tmp
        field[i, j] = field[i, j] + tmp2

def data_injection_transp( field, vel, samples, recloc, dt, stencil ):
  if np.sum(samples) != 0.0:
    for k in range( 0, recloc.shape[0] ):
        imin = max(stencil, recloc[k] - stencil)
        imax = min(recloc[k] + stencil, vel.shape[0] + stencil)
        #print( imin, imax)
        jmin = stencil
        jmax = 2 * stencil
        for i in range(imin, imax):
          for j in range(jmin, jmax):
            dist = (j - jmin)**2 + (i - imin)**2
            dist = np.exp(-dist)
            tmp = vel[i - stencil, j - stencil] * dt
            tmp2 = samples[k] * dist * tmp * tmp
            field[i, j] = field[i, j] + tmp2

def fwd_step2D( p1, p2, vel, params, cfs ):
  stencil = params[0]
  tmp = np.zeros((params[2] + 2*stencil), dtype=np.float32)
  for x in range( stencil, params[2] + stencil ):
    for z in range( stencil, params[1] + stencil ):
      tmp[ z ] = cfs[0] * p2[ z, x ] + \
                 cfs[1] * ( p2[ z, x + 1 ] + p2[ z, x - 1 ] ) + \
                 cfs[2] * ( p2[ z, x + 2 ] + p2[ z, x - 2 ] ) + \
                 cfs[3] * ( p2[ z, x + 3 ] + p2[ z, x - 3 ] ) + \
                 cfs[4] * ( p2[ z, x + 4 ] + p2[ z, x - 4 ] )
    for z in range( stencil, params[1] + stencil ):
      tmp2 = vel[ z - stencil, x - stencil ] * cfs[15]
      tmp2 *= tmp2
      tmp[ z ] = tmp2 * ( tmp[ z ] + \
                  cfs[5] * ( p2[ z + 1, x ] + p2[ z - 1, x ] ) + \
                  cfs[6] * ( p2[ z + 2, x ] + p2[ z - 2, x ] ) + \
                  cfs[7] * ( p2[ z + 3, x ] + p2[ z - 3, x ] ) + \
                  cfs[8] * ( p2[ z + 4, x ] + p2[ z - 4, x ] ) )
      p1[ z, x ] = p2[ z, x ] + p2[ z, x ] + tmp[ z ] - p1[ z, x ]

def forward( case, shot, vel, params, coeffs, wavelet ):
  stencil = params[0]
  data = np.zeros( (case.nr, case.nt), dtype=np.float32 )
  nxe = params[2] + 2*stencil
  pin = np.zeros( ( nxe, params[9] ), dtype=np.float32 )
  pout = np.zeros( ( nxe, params[9] ), dtype=np.float32 )
  waves = np.zeros( ( int(case.nt/case.skip), nxe, params[9] ), dtype=np.float32 )

  mylibc = ctypes.cdll.LoadLibrary("./libprops.so")
  fpointer = ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")
  ipointer = ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")
  mylibc.fwd_step2D.argtypes = [ fpointer, fpointer, fpointer, ipointer, fpointer ]
  mylibc.fwd_step2D.restype = None
  mylibc.absmult.argtypes = [ fpointer, fpointer, fpointer, ipointer ]
  mylibc.absmult.restype = None
  vel = np.transpose( vel )
  vel = np.ascontiguousarray(vel, dtype=np.float32)

  # Computing absorbing boundary conditions
  abc = absorb( nxe, params[9], False )

  # Wave propagation
  for it in range(case.nt):
    if it%100 == 0:
      print( "It: %d" % it )
      #plt.matshow( pin.T )
      #plt.matshow( vel.T )
      #plt.show()
    if it%case.skip == 0:
        waves[ int(it/10), :, : ] = pin[:,:]

    # calling C code, stencil
    mylibc.fwd_step2D( pin, pout, vel, params, coeffs )
    #fwd_step2D( pin, pout, vel, params, coeffs )
    # calling boundary condition
    mylibc.absmult( pin, pout, abc, params )
    #pin *= abc
    #pout *= abc
    # Source injection
    shot_injection_transp( pin, vel, wavelet[it], case.sloc[shot], coeffs[15], stencil )
    # recording
    data[ :, it ] = pin[ case.rloc, stencil ]

    tmp = pin
    pin = pout
    pout = tmp

  return data, waves

def modeling( model, case ):
  traces = np.zeros( (case.ns, case.nr, case.nt), dtype=np.float32 )

  # Simulation set up
  ( params, coeffs ) = prep( case, model )
  nxe = params[2] + 2*params[0]
  waves = np.zeros( ( case.ns, int(case.nt/case.skip), nxe, params[9] ), dtype=np.float32 )
  # source computation
  wavelet = make_wavelet( case, coeffs[15] )

  # Loop over shots
  for shot in range(case.ns):
    print( "Processing shot %d" % shot )
    traces[ shot, :, : ], waves[ shot, :, :, : ] = forward( case, shot, model.vel, params, coeffs, wavelet )

  return traces, waves

def illuminating( model, case ):

  # Simulation set up
  ( params, coeffs ) = prep( case, model )
  nxe = params[2] + 2*params[0]
  imaps = np.zeros( ( case.ns + 1, nxe, params[9] ), dtype=np.float32 )
  # source computation
  wavelet = make_wavelet( case, coeffs[15] )

  # Loop over shots
  for shot in range(case.ns):
    print( "Processing shot %d" % shot )
    imaps[ shot, :, : ] = forward_ill( case, shot, model.vel, params, coeffs, wavelet )
    imaps[ -1, :, : ] += imaps[ shot, :, : ]

  return imaps

def forward_ill( case, shot, vel, params, coeffs, wavelet ):
  stencil = params[0]
  nxe = params[2] + 2*stencil
  pin = np.zeros( ( nxe, params[9] ), dtype=np.float32 )
  pout = np.zeros( ( nxe, params[9] ), dtype=np.float32 )
  imap = np.zeros( ( nxe, params[9] ), dtype=np.float32 )

  mylibc = ctypes.cdll.LoadLibrary("./libprops.so")
  fpointer = ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")
  ipointer = ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")
  mylibc.fwd_step2D.argtypes = [ fpointer, fpointer, fpointer, ipointer, fpointer ]
  mylibc.fwd_step2D.restype = None
  vel = np.transpose( vel )
  vel = np.ascontiguousarray(vel, dtype=np.float32)

  # Computing absorbing boundary conditions
  abc = absorb( nxe, params[9], False )

  # Wave propagation
  for it in range(case.nt):
    if it%100 == 0:
      print( "It: %d" % it )
    if it%case.skip == 0:
        imap += pin[:,:] * pin[:,:]

    # calling C code, stencil
    mylibc.fwd_step2D( pin, pout, vel, params, coeffs )
    # calling boundary condition
    pin *= abc
    pout *= abc
    # Source injection
    shot_injection_transp( pin, vel, wavelet[it], case.sloc[shot], coeffs[15], stencil )

    tmp = pin
    pin = pout
    pout = tmp

  return imap

# Simple absorbing boundary condition
# Ideas taken from Daniel Koehn github repo
# fabs = exp(-a^2 (FW -i)^2), from Cerjan et al. 1985
def absorb( x, z, top ):
    FW = 50
    a = 0.0053

    coeff = np.zeros( FW, dtype=np.float32 )

    for i in range(FW):
        coeff[i] = np.exp(-(a**2 * (FW-i)**2))

    acoeffs = np.ones( (x, z), dtype=np.float32 )

    # letf and right of 2D model
    for i in range(FW):
        ze = z - i - 1
        ii = x - i - 1
        for j in range(0, ze):
            acoeffs[ii,j] = coeff[i]
            acoeffs[i,j] = coeff[i]

    # top and bottom of 2D model
    if top == False:
         xb = 0
         for j in range(FW):
             jj = z - j - 1
             xe = x - j 
             xb = j
             for i in range(xb, xe):
                 acoeffs[i,jj] = coeff[j]
    else:
         xb = 0
         for j in range(FW):
             jj = z - j - 1
             xe = x - j 
             xb = j
             for i in range(xb, xe):
                 acoeffs[i,jj] = coeff[j]
                 acoeffs[i,j] = coeff[j]

    return acoeffs

def migrating( model, case ):

  # Simulation set up
  ( params, coeffs ) = prep( case, model )
  nxe = params[2] + 2*params[0]
  simaps = np.zeros( ( case.ns + 1, nxe, params[9] ), dtype=np.float32 )
  rimaps = np.zeros( ( case.ns + 1, nxe, params[9] ), dtype=np.float32 )
  images = np.zeros( ( case.ns + 1, nxe, params[9] ), dtype=np.float32 )
  # source computation
  wavelet = make_wavelet( case, coeffs[15] )

  # Loop over shots
  for shot in range(case.ns):
    print( "Processing shot %d" % shot )
    (simaps[ shot, :, : ], rimaps[ shot, :, : ], images[ shot, :, : ]) = mig_shot( case, shot, model.vel, params, coeffs, wavelet )
    simaps[ -1, :, : ] += simaps[ shot, :, : ]
    rimaps[ -1, :, : ] += rimaps[ shot, :, : ]
    images[ -1, :, : ] += images[ shot, :, : ]

  return simaps, rimaps, images

def mig_shot( case, shot, vel, params, coeffs, wavelet ):
  # calling forward
  ( data, waves ) = forward( case, shot, vel, params, coeffs, wavelet )

  stencil = params[0]
  nxe = params[2] + 2*stencil
  pin = np.zeros( ( nxe, params[9] ), dtype=np.float32 )
  pout = np.zeros( ( nxe, params[9] ), dtype=np.float32 )
  image = np.zeros( ( nxe, params[9] ), dtype=np.float32 )
  rimap = np.zeros( ( nxe, params[9] ), dtype=np.float32 )
  simap = np.zeros( ( nxe, params[9] ), dtype=np.float32 )

  mylibc = ctypes.cdll.LoadLibrary("./libprops.so")
  fpointer = ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")
  ipointer = ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")
  mylibc.fwd_step2D.argtypes = [ fpointer, fpointer, fpointer, ipointer, fpointer ]
  mylibc.fwd_step2D.restype = None
  mylibc.mmult.argtypes = [ fpointer, fpointer, fpointer, ipointer, fpointer, fpointer ]
  mylibc.mmult.restype = None
  mylibc.absmult.argtypes = [ fpointer, fpointer, fpointer, ipointer ]
  mylibc.absmult.restype = None
  vel = np.transpose( vel )
  vel = np.ascontiguousarray(vel, dtype=np.float32)

  # Computing absorbing boundary conditions
  abc = absorb( nxe, params[9], False )

  # Wave propagation
  for it in range(case.nt):
    if it%100 == 0:
      print( "It: %d" % it )
    if it%case.skip == 0:
        mylibc.mmult( pin, waves[ waves.shape[0] - int(it/case.skip) - 1, :, : ], image, params, simap, rimap )

    # calling C code, stencil
    mylibc.fwd_step2D( pin, pout, vel, params, coeffs )
    #fwd_step2D( pin, pout, vel, params, coeffs )
    # calling boundary condition
    mylibc.absmult( pin, pout, abc, params )
    #pin *= abc
    #pout *= abc
    # Source injection
    data_injection_transp( pin, vel, data[ :, case.nt - it - 1], case.rloc, coeffs[15], stencil )

    tmp = pin
    pin = pout
    pout = tmp

  return simap, rimap, image

