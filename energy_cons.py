# Energy Conservation Crosscheck
import faulthandler
import numpy as np
import matplotlib.pyplot as plt
import math as mt
import Ellipsoid as Ellipse
import Test3D as Vox
import Implicit_steady_state_heat_transfer_solver as steadyStateHeatTransfer
# import sys
# import matplotlib
faulthandler.enable()
# matplotlib.use('Agg')
# t_count = 0
voxel_mat = Vox.dom
G = Ellipse.G
ua_test_voxel = Vox.UA_v
ua_test_tetra = steadyStateHeatTransfer.UA  # .todense()
Q_test_tetra = steadyStateHeatTransfer.Q
Q_test_voxel = Vox.Q


def perp_diffmat(vc1, t1, vc2, t2, dx, G, k1, h):
    L = mt.sqrt((G[vc1, t1].x - G[vc2, t2].x) ** 2 + (G[vc1, t1].y - G[vc2, t2].y) ** 2 + (
            G[vc1, t1].z - G[vc2, t2].z) ** 2) / 2
    ua_perp = 1 / (L / k1 + 1 / h) * (0.25 * dx ** 2)
    #    global t_count
    #    t_count = t_count + 1

    # print("perp")
    # print(L)
    return ua_perp, L


def side_diffmat(vc1, t1, vc2, t2, dx, G, k1, h):
    L = mt.sqrt((G[vc1, t1].x - G[vc2, t2].x) ** 2 + (G[vc1, t1].y - G[vc2, t2].y) ** 2 + (
            G[vc1, t1].z - G[vc2, t2].z) ** 2) / 2
    us_side = 1 / (L / k1 + 1 / h) * (0.17677 * dx ** 2)
    #    global t_count
    #    t_count = t_count + 1
    # print("side")
    return us_side, L


def diag_diffmat(vc1, t1, vc2, t2, dx, G, k1, h):
    L = mt.sqrt((G[vc1, t1].x - G[vc2, t2].x) ** 2 + (G[vc1, t1].y - G[vc2, t2].y) ** 2 + (
            G[vc1, t1].z - G[vc2, t2].z) ** 2) / 2
    ua_diag = 1 / (L / k1 + 1 / h) * (0.3535 * dx ** 2)
    #    global t_count
    #    t_count = t_count + 1
    # print("diag")
    return ua_diag, L


Temp = steadyStateHeatTransfer.Tdomtetra
Tvox = Vox.Tdom

Ts = []
Tsv = []

q_tetra = []
q_vox = []

k = steadyStateHeatTransfer.k
HTC = steadyStateHeatTransfer.HTC


def N1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("N1")
    if voxel_db[p + nx * ny, 6].mat != mat:  # S1
        ua, L = perp_diffmat(p, 0, p + nx * ny, 4, dx, G, k, HTC)
        Tt = (Temp[p, 0] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 0] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 3].mat != mat:  # N2
        #        # Ts.append(Temp[p,0])
        ua, L = side_diffmat(p, 0, p, 1, dx, G, k, HTC)
        Tt = (Temp[p, 0] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 0] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 5].mat != mat:  # N4
        #        # Ts.append(Temp[p,0])
        ua, L = side_diffmat(p, 0, p, 3, dx, G, k, HTC)
        Tt = (Temp[p, 0] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 0] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 22].mat != mat:  # B4
        #        # Ts.append(Temp[p,0])
        ua, L = side_diffmat(p, 0, p, 20, dx, G, k, HTC)
        Tt = (Temp[p, 0] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 0] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def N2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("N2")
    if voxel_db[p + nx * ny, 7].mat != mat:  # S2
        #        # Ts.append(Temp[p,1])
        ua, L = perp_diffmat(p, 1, p + nx * ny, 5, dx, G, k, HTC)
        Tt = (Temp[p, 1] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 1] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 4].mat != mat:  # N3
        #        # Ts.append(Temp[p,1])
        ua, L = side_diffmat(p, 1, p, 2, dx, G, k, HTC)
        Tt = (Temp[p, 1] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 1] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 2].mat != mat:  # N1
        #        # Ts.append(Temp[p,1])
        ua, L = side_diffmat(p, 1, p, 0, dx, G, k, HTC)
        Tt = (Temp[p, 1] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 1] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 10].mat != mat:  # W1
        #        # Ts.append(Temp[p,1])
        ua, L = diag_diffmat(p, 1, p, 8, dx, G, k, HTC)
        Tt = (Temp[p, 1] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 1] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def N3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("N3")
    if voxel_db[p + nx * ny, 8].mat != mat:  # S3
        #        # Ts.append(Temp[p,2])
        ua, L = perp_diffmat(p, 2, p + nx * ny, 6, dx, G, k, HTC)
        Tt = (Temp[p, 2] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 2] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 5].mat != mat:  # N4
        #        # Ts.append(Temp[p,2])
        ua, L = side_diffmat(p, 2, p, 3, dx, G, k, HTC)
        Tt = (Temp[p, 2] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 2] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 3].mat != mat:  # N2
        #        # Ts.append(Temp[p,2])
        ua, L = side_diffmat(p, 2, p, 1, dx, G, k, HTC)
        Tt = (Temp[p, 2] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 2] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 18].mat != mat:  # F1
        #        # Ts.append(Temp[p,2])
        ua, L = diag_diffmat(p, 2, p, 16, dx, G, k, HTC)
        Tt = (Temp[p, 2] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 2] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def N4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("N4")
    if voxel_db[p + nx * ny, 9].mat != mat:  # S4
        # # Ts.append(Temp[p,3])
        ua, L = perp_diffmat(p, 3, p + nx * ny, 7, dx, G, k, HTC)
        Tt = (Temp[p, 3] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 3] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 2].mat != mat):  # N1
        # # Ts.append(Temp[p,3])
        ua, L = side_diffmat(p, 3, p, 0, dx, G, k, HTC)
        Tt = (Temp[p, 3] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 3] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 4].mat != mat):  # N3
        # Ts.append(Temp[p,3])
        ua, L = side_diffmat(p, 3, p, 2, dx, G, k, HTC)
        Tt = (Temp[p, 3] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 3] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 14].mat != mat):  # E1
        # Ts.append(Temp[p,3])
        ua, L = diag_diffmat(p, 3, p, 12, dx, G, k, HTC)
        Tt = (Temp[p, 3] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 3] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def S1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("S1")
    if (voxel_db[p - nx * ny, 2].mat != mat):  # N1
        # Ts.append(Temp[p,4])
        ua, L = perp_diffmat(p, 4, p - nx * ny, 0, dx, G, k, HTC)
        Tt = (Temp[p, 4] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 4] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 9].mat != mat):  # S4
        # Ts.append(Temp[p,4])
        ua, L = side_diffmat(p, 4, p, 7, dx, G, k, HTC)
        Tt = (Temp[p, 4] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 4] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 7].mat != mat):  # S2
        # Ts.append(Temp[p,4])
        ua, L = side_diffmat(p, 4, p, 5, dx, G, k, HTC)
        Tt = (Temp[p, 4] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 4] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 24].mat != mat):  # B3
        # Ts.append(Temp[p,4])
        ua, L = diag_diffmat(p, 4, p, 22, dx, G, k, HTC)
        Tt = (Temp[p, 4] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 4] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def S2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("S2")
    if (voxel_db[p - nx * ny, 3].mat != mat):  # N2
        # Ts.append(Temp[p,5])
        ua, L = perp_diffmat(p, 5, p - nx * ny, 1, dx, G, k, HTC)
        Tt = (Temp[p, 5] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 5] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 9].mat != mat):  # S4
        # Ts.append(Temp[p,5])
        ua, L = side_diffmat(p, 5, p, 7, dx, G, k, HTC)
        Tt = (Temp[p, 5] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 5] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 8].mat != mat):  # S3
        # Ts.append(Temp[p,5])
        ua, L = side_diffmat(p, 5, p, 6, dx, G, k, HTC)
        Tt = (Temp[p, 5] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 5] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 12].mat != mat):  # W3
        # Ts.append(Temp[p,5])
        ua, L = diag_diffmat(p, 5, p, 10, dx, G, k, HTC)
        Tt = (Temp[p, 5] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 5] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def S3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("S3")
    if (voxel_db[p - nx * ny, 4].mat != mat):  # N3
        # Ts.append(Temp[p,6])
        ua, L = perp_diffmat(p, 6, p - nx * ny, 2, dx, G, k, HTC)
        Tt = (Temp[p, 6] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 6] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 9].mat != mat):  # S4
        # Ts.append(Temp[p,6])
        ua, L = side_diffmat(p, 6, p, 7, dx, G, k, HTC)
        Tt = (Temp[p, 6] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 6] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 7].mat != mat):  # S2
        # Ts.append(Temp[p,6])
        ua, L = side_diffmat(p, 6, p, 5, dx, G, k, HTC)
        Tt = (Temp[p, 6] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 6] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 20].mat != mat):  # F3
        # Ts.append(Temp[p,6])
        ua, L = diag_diffmat(p, 6, p, 18, dx, G, k, HTC)
        q = (Temp[p, 6] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def S4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("S4")
    if (voxel_db[p - nx * ny, 5].mat != mat):  # N4
        # Ts.append(Temp[p,7])
        ua, L = perp_diffmat(p, 7, p - nx * ny, 3, dx, G, k, HTC)
        Tt = (Temp[p, 7] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 7] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 6].mat != mat):  # S1
        # Ts.append(Temp[p,7])
        ua, L = side_diffmat(p, 7, p, 4, dx, G, k, HTC)
        Tt = (Temp[p, 7] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 7] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 8].mat != mat):  # S3
        # Ts.append(Temp[p,7])
        ua, L = side_diffmat(p, 7, p, 6, dx, G, k, HTC)
        Tt = (Temp[p, 7] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 7] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 16].mat != mat):  # E3
        # Ts.append(Temp[p,7])
        ua, L = diag_diffmat(p, 7, p, 14, dx, G, k, HTC)
        Tt = (Temp[p, 7] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 7] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def W1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("W1")
    if (voxel_db[p - ny, 14].mat != mat):  # E1
        # Ts.append(Temp[p,8])
        ua, L = perp_diffmat(p, 8, p - ny, 12, dx, G, k, HTC)
        Tt = (Temp[p, 8] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 8] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 11].mat != mat):  # W2
        # Ts.append(Temp[p,8])
        ua, L = side_diffmat(p, 8, p, 9, dx, G, k, HTC)
        Tt = (Temp[p, 8] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 8] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 13].mat != mat):  # W4
        # Ts.append(Temp[p,8])
        ua, L = side_diffmat(p, 8, p, 11, dx, G, k, HTC)
        Tt = (Temp[p, 8] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 8] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 3].mat != mat):  # N2
        # Ts.append(Temp[p,8])
        ua, L = diag_diffmat(p, 8, p, 1, dx, G, k, HTC)
        Tt = (Temp[p, 8] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 8] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def W2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("W2")
    if (voxel_db[p - ny, 17].mat != mat):  # E4
        # Ts.append(Temp[p,9])
        ua, L = perp_diffmat(p, 9, p - ny, 15, dx, G, k, HTC)
        Tt = (Temp[p, 9] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 9] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 12].mat != mat):  # W3
        # Ts.append(Temp[p,9])
        ua, L = side_diffmat(p, 9, p, 10, dx, G, k, HTC)
        Tt = (Temp[p, 9] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 9] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 10].mat != mat):  # W1
        # Ts.append(Temp[p,9])
        ua, L = side_diffmat(p, 9, p, 8, dx, G, k, HTC)
        Tt = (Temp[p, 9] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 9] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 25].mat != mat):  # B4
        # Ts.append(Temp[p,9])
        ua, L = diag_diffmat(p, 9, p, 23, dx, G, k, HTC)
        Tt = (Temp[p, 9] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 9] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def W3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("W3")
    if (voxel_db[p - ny, 16].mat != mat):  # E3
        # Ts.append(Temp[p,10])
        ua, L = perp_diffmat(p, 10, p - ny, 14, dx, G, k, HTC)
        Tt = (Temp[p, 10] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 10] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 13].mat != mat):  # W4
        # Ts.append(Temp[p,10])
        ua, L = side_diffmat(p, 10, p, 11, dx, G, k, HTC)
        Tt = (Temp[p, 10] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 10] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 11].mat != mat):  # W2
        # Ts.append(Temp[p,10])
        ua, L = side_diffmat(p, 10, p, 9, dx, G, k, HTC)
        Tt = (Temp[p, 10] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 10] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 7].mat != mat):  # S2
        # Ts.append(Temp[p,10])
        ua, L = diag_diffmat(p, 10, p, 5, dx, G, k, HTC)
        Tt = (Temp[p, 10] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 10] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def W4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("W4")
    if (voxel_db[p - ny, 15].mat != mat):  # E2
        # Ts.append(Temp[p,11])
        ua, L = perp_diffmat(p, 11, p - ny, 13, dx, G, k, HTC)
        Tt = (Temp[p, 11] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 11] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 10].mat != mat):  # W1
        # Ts.append(Temp[p,11])
        ua, L = side_diffmat(p, 11, p, 8, dx, G, k, HTC)
        Tt = (Temp[p, 11] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 11] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 12].mat != mat):  # W3
        # Ts.append(Temp[p,11])
        ua, L = side_diffmat(p, 11, p, 10, dx, G, k, HTC)
        Tt = (Temp[p, 11] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 11] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 19].mat != mat):  # F2
        # Ts.append(Temp[p,11])
        ua, L = diag_diffmat(p, 11, p, 17, dx, G, k, HTC)
        Tt = (Temp[p, 11] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 11] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def E1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("E1")
    if (voxel_db[p + ny, 10].mat != mat):  # W1
        # Ts.append(Temp[p,12])
        ua, L = perp_diffmat(p, 12, p + ny, 8, dx, G, k, HTC)
        Tt = (Temp[p, 12] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 12] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 15].mat != mat):  # E2
        # Ts.append(Temp[p,12])
        ua, L = side_diffmat(p, 12, p, 13, dx, G, k, HTC)
        Tt = (Temp[p, 12] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 12] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 17].mat != mat):  # E4
        # Ts.append(Temp[p,12])
        ua, L = side_diffmat(p, 12, p, 15, dx, G, k, HTC)
        Tt = (Temp[p, 12] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 12] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 5].mat != mat):  # N4
        # Ts.append(Temp[p,12])
        ua, L = diag_diffmat(p, 12, p, 3, dx, G, k, HTC)
        Tt = (Temp[p, 12] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 12] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def E2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("E2")
    if (voxel_db[p + ny, 13].mat != mat):  # W4
        # Ts.append(Temp[p,13])
        ua, L = perp_diffmat(p, 13, p + ny, 11, dx, G, k, HTC)
        Tt = (Temp[p, 13] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 13] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 16].mat != mat):  # E3
        # Ts.append(Temp[p,13])
        ua, L = side_diffmat(p, 13, p, 14, dx, G, k, HTC)
        Tt = (Temp[p, 13] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 13] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 14].mat != mat):  # E1
        # Ts.append(Temp[p,13])
        ua, L = side_diffmat(p, 13, p, 12, dx, G, k, HTC)
        Tt = (Temp[p, 13] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 13] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 21].mat != mat):  # F4
        # Ts.append(Temp[p,13])
        ua, L = diag_diffmat(p, 13, p, 19, dx, G, k, HTC)
        Tt = (Temp[p, 13] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 13] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def E3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("E3")
    if (voxel_db[p + ny, 12].mat != mat):  # W3
        # Ts.append(Temp[p,14])
        ua, L = perp_diffmat(p, 14, p + ny, 10, dx, G, k, HTC)
        Tt = (Temp[p, 14] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 14] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 17].mat != mat):  # E4
        # Ts.append(Temp[p,14])
        ua, L = side_diffmat(p, 14, p, 15, dx, G, k, HTC)
        Tt = (Temp[p, 14] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 14] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 15].mat != mat):  # E2
        # Ts.append(Temp[p,14])
        ua, L = side_diffmat(p, 14, p, 13, dx, G, k, HTC)
        Tt = (Temp[p, 14] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 14] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 9].mat != mat):  # S4
        # Ts.append(Temp[p,14])
        ua, L = diag_diffmat(p, 14, p, 7, dx, G, k, HTC)
        Tt = (Temp[p, 14] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 14] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def E4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("E4")
    if (voxel_db[p + ny, 11].mat != mat):  # W2
        # Ts.append(Temp[p,15])
        ua, L = perp_diffmat(p, 15, p + ny, 9, dx, G, k, HTC)
        Tt = (Temp[p, 15] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 15] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 14].mat != mat):  # E1
        # Ts.append(Temp[p,15])
        ua, L = side_diffmat(p, 15, p, 12, dx, G, k, HTC)
        Tt = (Temp[p, 15] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 15] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 16].mat != mat):  # E3
        # Ts.append(Temp[p,15])
        ua, L = side_diffmat(p, 15, p, 14, dx, G, k, HTC)
        Tt = (Temp[p, 15] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 15] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 23].mat != mat):  # B2
        # Ts.append(Temp[p,15])
        ua, L = diag_diffmat(p, 15, p, 21, dx, G, k, HTC)
        Tt = (Temp[p, 15] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 15] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def F1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("F1")
    if (voxel_db[p + 1, 22].mat != mat):  # B1
        # Ts.append(Temp[p,16])
        ua, L = perp_diffmat(p, 16, p + 1, 20, dx, G, k, HTC)
        Tt = (Temp[p, 16] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 16] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 19].mat != mat):  # F2
        # Ts.append(Temp[p,16])
        ua, L = side_diffmat(p, 16, p, 17, dx, G, k, HTC)
        Tt = (Temp[p, 16] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 16] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 21].mat != mat):  # F4
        # Ts.append(Temp[p,16])
        ua, L = side_diffmat(p, 16, p, 19, dx, G, k, HTC)
        Tt = (Temp[p, 16] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 16] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 4].mat != mat):  # N3
        # Ts.append(Temp[p,16])
        ua, L = diag_diffmat(p, 16, p, 2, dx, G, k, HTC)
        Tt = (Temp[p, 16] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 16] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def F2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("F2")
    if (voxel_db[p + 1, 25].mat != mat):  # B4
        # Ts.append(Temp[p,17])
        ua, L = perp_diffmat(p, 17, p + 1, 23, dx, G, k, HTC)
        Tt = (Temp[p, 17] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 17] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 20].mat != mat):  # F3
        # Ts.append(Temp[p,17])
        ua, L = side_diffmat(p, 17, p, 18, dx, G, k, HTC)
        Tt = (Temp[p, 17] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 17] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 18].mat != mat):  # F1
        # Ts.append(Temp[p,17])
        ua, L = side_diffmat(p, 17, p, 16, dx, G, k, HTC)
        Tt = (Temp[p, 17] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 17] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 13].mat != mat):  # W4
        # Ts.append(Temp[p,17])
        ua, L = diag_diffmat(p, 17, p, 11, dx, G, k, HTC)
        Tt = (Temp[p, 17] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 17] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def F3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("F3")
    if (voxel_db[p + 1, 24].mat != mat):  # B3
        # Ts.append(Temp[p,18])
        ua, L = perp_diffmat(p, 18, p + 1, 22, dx, G, k, HTC)
        Tt = (Temp[p, 18] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 18] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 21].mat != mat):  # F4
        # Ts.append(Temp[p,18])
        ua, L = side_diffmat(p, 18, p, 19, dx, G, k, HTC)
        Tt = (Temp[p, 18] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 18] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 19].mat != mat):  # F2
        # Ts.append(Temp[p,18])
        ua, L = side_diffmat(p, 18, p, 17, dx, G, k, HTC)
        Tt = (Temp[p, 18] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 18] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 8].mat != mat):  # S3
        # Ts.append(Temp[p,18])
        ua, L = diag_diffmat(p, 18, p, 6, dx, G, k, HTC)
        Tt = (Temp[p, 18] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 18] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def F4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("F4")
    if (voxel_db[p + 1, 23].mat != mat):  # B2
        # Ts.append(Temp[p,19])
        ua, L = perp_diffmat(p, 19, p + 1, 21, dx, G, k, HTC)
        Tt = (Temp[p, 19] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 19] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 18].mat != mat):  # F1
        # Ts.append(Temp[p,19])
        ua, L = side_diffmat(p, 19, p, 16, dx, G, k, HTC)
        Tt = (Temp[p, 19] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 19] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 20].mat != mat):  # F3
        # Ts.append(Temp[p,19])
        ua, L = side_diffmat(p, 19, p, 18, dx, G, k, HTC)
        Tt = (Temp[p, 19] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 19] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 15].mat != mat):  # E2
        # Ts.append(Temp[p,19])
        ua, L = diag_diffmat(p, 19, p, 13, dx, G, k, HTC)
        Tt = (Temp[p, 19] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 19] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def B1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("B1")
    if (voxel_db[p - 1, 18].mat != mat):  # F1
        # Ts.append(Temp[p,20])
        ua, L = perp_diffmat(p, 20, p - 1, 16, dx, G, k, HTC)
        Tt = (Temp[p, 20] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 20] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 23].mat != mat):  # B2
        # Ts.append(Temp[p,20])
        ua, L = side_diffmat(p, 20, p, 21, dx, G, k, HTC)
        Tt = (Temp[p, 20] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 20] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 25].mat != mat):  # B4
        # Ts.append(Temp[p,20])
        ua, L = side_diffmat(p, 20, p, 23, dx, G, k, HTC)
        Tt = (Temp[p, 20] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 20] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if (voxel_db[p, 2].mat != mat):  # N1
        # Ts.append(Temp[p,20])
        ua, L = diag_diffmat(p, 20, p, 0, dx, G, k, HTC)
        Tt = (Temp[p, 20] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 20] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def B2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("B2")
    if voxel_db[p - 1, 21].mat != mat:  # F4
        # Ts.append(Temp[p,21])
        ua, L = perp_diffmat(p, 21, p - 1, 19, dx, G, k, HTC)
        Tt = (Temp[p, 21] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 21] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 24].mat != mat:  # B3
        # Ts.append(Temp[p,21])
        ua, L = side_diffmat(p, 21, p, 22, dx, G, k, HTC)
        Tt = (Temp[p, 21] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 21] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 22].mat != mat:  # B1
        # Ts.append(Temp[p,21])
        ua, L = side_diffmat(p, 21, p, 20, dx, G, k, HTC)
        Tt = (Temp[p, 21] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 21] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 17].mat != mat:  # E4
        # Ts.append(Temp[p,21])
        ua, L = diag_diffmat(p, 21, p, 15, dx, G, k, HTC)
        Tt = (Temp[p, 21] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 21] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def B3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("B3")
    if voxel_db[p - 1, 20].mat != mat:  # F3
        # Ts.append(Temp[p,22])
        ua, L = perp_diffmat(p, 22, p - 1, 18, dx, G, k, HTC)
        Tt = (Temp[p, 22] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 22] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 25].mat != mat:  # B4
        # Ts.append(Temp[p,22])
        ua, L = side_diffmat(p, 22, p, 23, dx, G, k, HTC)
        Tt = (Temp[p, 22] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 22] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 23].mat != mat:  # B2
        # Ts.append(Temp[p,22])
        ua, L = side_diffmat(p, 22, p, 21, dx, G, k, HTC)
        Tt = (Temp[p, 22] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 22] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 6].mat != mat:  # S3
        # Ts.append(Temp[p,22])
        ua, L = diag_diffmat(p, 22, p, 4, dx, G, k, HTC)
        Tt = (Temp[p, 22] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 22] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def B4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    #    print("B4")
    if voxel_db[p - 1, 19].mat != mat:  # F2
        # Ts.append(Temp[p,23])
        ua, L = perp_diffmat(p, 23, p - 1, 17, dx, G, k, HTC)
        Tt = (Temp[p, 23] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 23] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 22].mat != mat:  # B1
        # Ts.append(Temp[p,23])
        ua, L = side_diffmat(p, 23, p, 20, dx, G, k, HTC)
        Tt = (Temp[p, 23] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 23] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 24].mat != mat:  # B3
        # Ts.append(Temp[p,23])
        ua, L = side_diffmat(p, 23, p, 22, dx, G, k, HTC)
        Tt = (Temp[p, 23] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 23] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)
    if voxel_db[p, 11].mat != mat:  # W2
        # Ts.append(Temp[p,23])
        ua, L = diag_diffmat(p, 23, p, 9, dx, G, k, HTC)
        Tt = (Temp[p, 23] + HTC * L / k * steadyStateHeatTransfer.Tamb) / (1 + HTC * L / k)
        Ts.append(Tt)
        q = (Temp[p, 23] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(q)


def func(tag, voxel_db, p, nx, ny, nz, dx, mat, k, HTC):
    if tag == 'N1':
        N1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'N2':
        N2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'N3':
        N3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'N4':
        N4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)

    elif tag == 'S1':
        S1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'S2':
        S2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'S3':
        S3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'S4':
        S4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)

    elif tag == 'W1':
        W1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'W2':
        W2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'W3':
        W3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'W4':
        W4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)

    elif tag == 'E1':
        E1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'E2':
        E2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'E3':
        E3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'E4':
        E4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)

    elif tag == 'F1':
        F1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'F2':
        F2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'F3':
        F3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'F4':
        F4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)

    elif tag == 'B1':
        B1(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'B2':
        B2(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'B3':
        B3(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)
    elif tag == 'B4':
        B4(voxel_db, p, nx, ny, nz, dx, mat, k, HTC)


for vc in range(Ellipse.voxel_n):
    for i in range(24):
        if Ellipse.voxel_db[vc, i + 2].mat != 0:
            func(Ellipse.voxel_db[vc, i + 2].pos, Ellipse.voxel_db, vc, Ellipse.nx, Ellipse.ny, Ellipse.nz, Ellipse.dx,
                 Ellipse.voxel_db[vc, i + 2].mat, steadyStateHeatTransfer.k, steadyStateHeatTransfer.HTC)

Tsa = steadyStateHeatTransfer.Heat_rate * Ellipse.a / (3 * steadyStateHeatTransfer.HTC) + steadyStateHeatTransfer.Tamb
# plt.plot(Ts)
# plt.axhline(Tsa,color='red')
# plt.show()

error = []
for i in range(len(Ts)):
    error.append(Tsa - Ts[i])

# plt.boxplot(error)
# plt.show()


# voxelplot = 1
# if(voxelplot == 1):

errv = []
for k in range(Vox.n):
    for j in range(Vox.n):
        for i in range(Vox.n):
            if Vox.dom[i, j, k] == 1:
                if Vox.side_exposed[i, j, k] > 0:
                    Ttv = (Vox.Tdom[i, j, k] + (steadyStateHeatTransfer.HTC * (Ellipse.dx / 2) / steadyStateHeatTransfer.k) * steadyStateHeatTransfer.Tamb) / (
                            1 + steadyStateHeatTransfer.HTC * (Ellipse.dx / 2) / steadyStateHeatTransfer.k)
                    Tsv.append(Ttv)
                    errv.append(Tsa - Ttv)
                    ua = 1 / (0.5 * Ellipse.dx / steadyStateHeatTransfer.k + 1 / steadyStateHeatTransfer.HTC) * Ellipse.dx * Ellipse.dx
                    q = Vox.side_exposed[i, j, k] * ua * (Vox.Tdom[i, j, k] - steadyStateHeatTransfer.Tamb)
                    q_vox.append(q)

data = [error, errv]
labels = ["Tetrahedrons", "Voxels"]

plt.figure(figsize=(10, 6))
plt.subplot(121)
bplot = plt.boxplot(data, vert=True, patch_artist=True, labels=labels)
colors = ['pink', 'lightblue']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
    plt.grid(color='grey')
    plt.ylabel("Surface Temperature Error in degrees Celsius")

x_t = []
x_v = []
t = 1
for i in range(len(Ts)):
    x_t.append((4 * mt.pi * 0.1 ** 2) / t)
    t = t + 1
t = 1
for i in range(len(Tsv)):
    x_v.append((4 * mt.pi * 0.1 ** 2) / t)
    t = t + 1

# plt.figure(figsize = (10,6))
plt.subplot(122)
# plt.plot(x_t,Ts,label = "Tetra")
# plt.plot(x_v,Tsv,label = "Voxel")
plt.axhline(np.mean(Ts), label="Mean Tetra", color="Green")
plt.axhline(np.mean(Tsv), label="Mean Voxel", color="Purple")
plt.axhline(Tsa, color="red", label="Analytical")
plt.yticks(np.arange(steadyStateHeatTransfer.Tamb - 2, Tsa + 2, 2))
plt.ylabel("Surface Temperature")
plt.legend()
# plt.savefig("Results_1mm_R10cm.png")
plt.show()
plt.show()

print("\nAnalytical Ts =", Tsa)
print("\nMean Ts tetra", np.mean(Ts))
print("Mean Ts Voxel", np.mean(Tsv), "\n")

print("Tetrahedron error =", steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume - np.sum(q_tetra))
print("Heart rate = ", steadyStateHeatTransfer.Heat_rate)
print("tetra_volume = ", Ellipse.tetra_volume)
print("percentage error = ",
      (steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume - np.sum(q_tetra)) / (steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume) * 100, "%")
print("Tetra T_surf for given", steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume / (Ellipse.areasum1 * steadyStateHeatTransfer.HTC) + steadyStateHeatTransfer.Tamb)

print("voxel error =", steadyStateHeatTransfer.Heat_rate * Ellipse.voxel_volume - np.sum(q_vox), "W")
print("Voxel heat error percentage =",
      (steadyStateHeatTransfer.Heat_rate * Ellipse.voxel_volume - np.sum(q_vox)) / steadyStateHeatTransfer.Heat_rate * Ellipse.voxel_volume * 100, "%")

print("\nActual tetra Q = ", steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume)
print("Tetra Q =  ", np.sum(q_tetra))
print("voxel Q =  ", np.sum(q_vox))

test_tetra = 0
for i in range(len(Q_test_tetra)):
    if Q_test_tetra[i] != steadyStateHeatTransfer.Tamb:
        test_tetra = test_tetra + Q_test_tetra[i]
print("Q_tetra_test =", test_tetra)

test_tetra = 0
for i in range(len(Q_test_voxel)):
    if Q_test_voxel[i] != steadyStateHeatTransfer.Tamb:
        test_tetra = test_tetra + Q_test_voxel[i]
print("Q_voxel_test =", test_tetra)
