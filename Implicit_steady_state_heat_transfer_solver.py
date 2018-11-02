import Ellipsoid as ellipse
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ua_matrix_generator as ua_m
#import xlsxwriter as xl
#from numpy import linalg as la
import scipy.sparse.linalg as spla
#import scipy.sparse as sp
import sys
import faulthandler; faulthandler.enable()


def Heatsolver(Heat_rate,Th_cond,HTC,Tamb):
    
    
    voxel_db = ellipse.voxel_db
    
    K = np.zeros((ellipse.voxel_n,24),dtype = float)

    for vc in range(ellipse.voxel_n):
        for i in range(24):
            if(voxel_db[vc,(i+2)].mat != 0):
                K[vc,i] = Th_cond
    
    nx = ellipse.nx
    ny = ellipse.ny
    nz = ellipse.nz
    
    n_tetra = nx*ny*nz*24
    
    tetra_vol = ellipse.dx*ellipse.dy*ellipse.dz/24
    
    q = tetra_vol*Heat_rate

    Q = np.zeros((n_tetra,1),dtype = float)
    
    # Creating the UA matrix
    
    UA,Q = ua_m.matrixgenerator(voxel_db,K,HTC,ellipse.dx,ellipse.dy,ellipse.dz,ellipse.voxel_n,ellipse.G,nx,ny,nz,q,Q,Tamb)
    

    return(UA,Q)



Heat_rate = 1000
HTC = 2.0
k = 0.3
Tamb = 20.0

UA,Q = Heatsolver(Heat_rate,k,HTC,Tamb)

print("Matrix generated, now solving it")

print("memory = ",sys.getsizeof(UA))



T = spla.spsolve(UA,Q)

Tdomtetra = np.zeros((ellipse.voxel_n,24),dtype = float)
count = 0
for vc in range(ellipse.voxel_n):
    for i in range(24):
        Tdomtetra[vc,i] = T[count]
        count = count + 1

