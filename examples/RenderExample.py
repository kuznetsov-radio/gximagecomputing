import ctypes
import numpy as np
import matplotlib.pyplot as plt

f=open('model.bin', mode='rb')
model=f.read()
f.close()

f=open('ebtel.bin', mode='rb')
ebtel=f.read()
f.close()

Nx=150
Ny=150
xc=-50.0
yc=390.0
dx=2.0
dy=2.0
freqlist=[5.8, 6.2, 6.6, 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 9.8, 10.2, 10.6, 11.0, 11.4, 11.8]

Tbase=1e6
nbase=1e8
Q0=4.5e-3
a=1.5
b=2.5
force_isothermal=1

#------------------------------------------------

dt_s=np.dtype([('Nx', np.int32),
               ('Ny', np.int32),
               ('Nf', np.int32),
               ('_r1', np.int32),
               ('xc', np.float64),
               ('yc', np.float64),
               ('dx', np.float64),
               ('dy', np.float64),
               ('freqlist', np.float64, len(freqlist))])
simbox=np.zeros(1, dtype=dt_s)
simbox['Nx']=Nx
simbox['Ny']=Ny
simbox['Nf']=len(freqlist)
simbox['xc']=xc
simbox['yc']=yc
simbox['dx']=dx
simbox['dy']=dy
simbox['freqlist']=freqlist

dt_c=np.dtype([('Tbase', np.float64),
               ('nbase', np.float64),
               ('Q0', np.float64),
               ('a', np.float64),
               ('b', np.float64),
               ('iso', np.int32)])
cparms=np.zeros(1, dtype=dt_c)
cparms['Tbase']=Tbase
cparms['nbase']=nbase
cparms['Q0']=Q0
cparms['a']=a
cparms['b']=b
cparms['iso']=force_isothermal

dt_o=np.dtype([('flagsAll', np.int32, 6),
               ('flagsCorona', np.int32, 6),
               ('TI', np.float64, (Nx, Ny, len(freqlist))),
               ('TV', np.float64, (Nx, Ny, len(freqlist)))])
outspace=np.zeros(1, dtype=dt_o)

libname='.\RenderGRFF_64.dll'
libc_mw=ctypes.CDLL(libname)
mwfunc=libc_mw.pyComputeMW
_dt_s=np.ctypeslib.ndpointer(dtype=dt_s)
_dt_c=np.ctypeslib.ndpointer(dtype=dt_c)
_dt_o=np.ctypeslib.ndpointer(dtype=dt_o)
mwfunc.argtypes=[ctypes.c_void_p, ctypes.c_void_p, _dt_s, _dt_c, _dt_o]
mwfunc.restype=ctypes.c_int

r=mwfunc(model, ebtel, simbox, cparms, outspace)

o=outspace[0]
TI=np.reshape(np.ravel(o['TI']), (Nx, Ny, len(freqlist)), order='F')
TV=np.reshape(np.ravel(o['TV']), (Nx, Ny, len(freqlist)), order='F')

#-----------------------------------------

plt.figure(1)
img=plt.imshow(np.transpose(TI[:, :, 0]), origin='lower')

plt.figure(2)
img=plt.imshow(np.transpose(TV[:, :, 0]), origin='lower')

plt.show()
