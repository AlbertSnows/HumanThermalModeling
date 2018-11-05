import numpy as np
# import math as mth
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import Ellipsoid as ellip
# from numpy import linalg as la
# import matplotlib.pyplot as plt
# import timeit as time
import Implicit_steady_state_heat_transfer_solver as isht

# start = time.default_timer()
i_cor_vox = []
j_cor_vox = []
UA_vox_data = []

R = ellip.a  # mth.sqrt(ellip.a**2 + ellip.b**2 + ellip.c**2)
dx = ellip.dx
n = ellip.nx
dom = ellip.voxel
# dom = np.zeros((n,n,n),dtype = int)
# cx = n/2
# cy = n/2
# cz = n/2
# for k in range(n):
#    for j in range(n):
#        for i in range(n):   
#            if(i>2 and i<n-2):
#                if(j>2 and j<n-2):
#                    if(k>2 and k<n-2):
#                        dom[i,j,k] = 1
#            if((i-cx)**2 + (j-cy)**2 + (k-cz)**2 <= (R/dx)**2):
#                    dom[i,j,k] = 1

k = isht.k
q = isht.Heat_rate
h = isht.HTC
Ta = isht.Tamb
total_n = ellip.nx * ellip.ny * ellip.nz

side_exposed = ellip.side_exposed
# side_exposed = np.zeros((n,n,n),dtype = int)
# for r in range(n):
#    for q in range(n):
#        for p in range(n):
#            count = 0
#            if(dom[p,q,r] != 0):
#                if(dom[p+1,q,r] == 0):
#                    count = count + 1
#                if(dom[p-1,q,r] == 0):
#                    count = count + 1
#                if(dom[p,q+1,r] == 0):
#                    count = count + 1
#                if(dom[p,q-1,r] == 0):
#                    count = count + 1
#                if(dom[p,q,r+1] == 0):
#                    count = count + 1
#                if(dom[p,q,r-1] == 0):
#                    count = count + 1
#                side_exposed[p,q,r] = count

UA_v = np.zeros((total_n, total_n), dtype=float)

Q = np.zeros((total_n, 1), dtype=float)

a = 1
b = ellip.ny
c = ellip.nx * ellip.ny

for g in range(ellip.nz):
    for j in range(ellip.ny):
        for i in range(ellip.nx):
            p = a * i + b * j + c * g
            if (dom[i, j, g] == 0):

                i_cor_vox.append(p);
                j_cor_vox.append(p);
                UA_vox_data.append(1.0)
                #                UA_v[p,p] = 1.0
                Q[p, 0] = Ta

            # if(dom[i,j,g] != 0):
            else:
                if (dom[i + 1, j, g] == 0):
                    uar = 1 / (1 / h + dx * 0.5 / k) * dx * dx
                    pr = a * (i + 1) + b * j + c * g
                if (dom[i - 1, j, g] == 0):
                    ual = 1 / (1 / h + dx * 0.5 / k) * dx * dx
                    pl = a * (i - 1) + b * j + c * g
                if (dom[i, j + 1, g] == 0):
                    uat = 1 / (1 / h + dx * 0.5 / k) * dx * dx
                    pt = a * i + b * (j + 1) + c * g
                if (dom[i, j - 1, g] == 0):
                    uab = 1 / (1 / h + dx * 0.5 / k) * dx * dx
                    pb = a * i + b * (j - 1) + c * g
                if (dom[i, j, g + 1] == 0):
                    uaf = 1 / (1 / h + dx * 0.5 / k) * dx * dx
                    pf = a * i + b * j + c * (g + 1)
                if (dom[i, j, g - 1] == 0):
                    uabk = 1 / (1 / h + dx * 0.5 / k) * dx * dx
                    pbk = i + b * j + c * (g - 1)

                if (dom[i + 1, j, g] != 0):
                    uar = 1 / (dx * 0.5 / k + dx * 0.5 / k) * dx * dx
                    pr = a * (i + 1) + b * j + c * g
                if (dom[i - 1, j, g] != 0):
                    ual = 1 / (dx * 0.5 / k + dx * 0.5 / k) * dx * dx
                    pl = a * (i - 1) + b * j + c * g
                if (dom[i, j + 1, g] != 0):
                    uat = 1 / (dx * 0.5 / k + dx * 0.5 / k) * dx * dx
                    pt = a * i + b * (j + 1) + c * g
                if (dom[i, j - 1, g] != 0):
                    uab = 1 / (dx * 0.5 / k + dx * 0.5 / k) * dx * dx
                    pb = a * i + b * (j - 1) + c * g
                if (dom[i, j, g + 1] != 0):
                    uaf = 1 / (dx * 0.5 / k + dx * 0.5 / k) * dx * dx
                    pf = a * i + b * j + c * (g + 1)
                if (dom[i, j, g - 1] != 0):
                    uabk = 1 / (dx * 0.5 / k + dx * 0.5 / k) * dx * dx
                    pbk = a * i + b * j + c * (g - 1)

                i_cor_vox.append(p);
                j_cor_vox.append(pr);
                UA_vox_data.append(uar)
                i_cor_vox.append(p);
                j_cor_vox.append(pl);
                UA_vox_data.append(ual)
                i_cor_vox.append(p);
                j_cor_vox.append(pt);
                UA_vox_data.append(uat)
                i_cor_vox.append(p);
                j_cor_vox.append(pb);
                UA_vox_data.append(uab)
                i_cor_vox.append(p);
                j_cor_vox.append(pf);
                UA_vox_data.append(uaf)
                i_cor_vox.append(p);
                j_cor_vox.append(pbk);
                UA_vox_data.append(uabk)
                i_cor_vox.append(p);
                j_cor_vox.append(p);
                UA_vox_data.append(-(uar + ual + uat + uab + uaf + uabk))

                #                UA_v[p,pr] = uar
                #                UA_v[p,pl] = ual
                #                UA_v[p,pt] = uat
                #                UA_v[p,pb] = uab
                #                UA_v[p,pf] = uaf
                #                UA_v[p,pbk] = uabk
                #                UA_v[p,p] = -(uar + ual + uat + uab + uaf + uabk)
                Q[p, 0] = -q * dx * dx * dx
            # print(i,j,g,p)

# M = int(len(Q))
UA_v = sp.csc_matrix((UA_vox_data, (i_cor_vox, j_cor_vox)))  # ,[shape = (M,M)])#

# det = (la.det(UA_v))
# print(det)

# UA_v = UA_csc_vox.todens()
Tv = spla.spsolve(UA_v, Q)

Tdom = np.zeros((ellip.nx, ellip.ny, ellip.nz), dtype=float)
count = 0
for g in range(ellip.nz):
    for j in range(ellip.ny):
        for i in range(ellip.nx):
            Tdom[i, j, g] = Tv[count]
            count = count + 1
