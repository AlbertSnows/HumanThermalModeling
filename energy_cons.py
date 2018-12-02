# Energy Conservation Crosscheck
import faulthandler
import numpy
import matplotlib.pyplot as pyplot
import math
import ellipsoid_optimized_full as Ellipse
import Test3D as Vox
import Implicit_steady_state_heat_transfer_solver as steadyStateHeatTransfer
# import sys
# import matplotlib
faulthandler.enable()
voxel_mat = Vox.dom
outer_g = Ellipse.g_coord
ua_test_voxel = Vox.UA_v
ua_test_tetra = steadyStateHeatTransfer.UA  # .todense()
Q_test_tetra = steadyStateHeatTransfer.Q
Q_test_voxel = Vox.Q
# matplotlib.use('Agg')
# t_count = 0


def perp_diffmat(vc1, t1, vc2, t2, dx, g, k1, h):
    l_val = math.sqrt((g[vc1, t1].x - g[vc2, t2].x) ** 2 + (g[vc1, t1].y - g[vc2, t2].y) ** 2 + (
            g[vc1, t1].z - g[vc2, t2].z) ** 2) / 2
    ua_perp = 1 / (l_val / k1 + 1 / h) * (0.25 * dx ** 2)
    #    global t_count
    #    t_count = t_count + 1

    # print("perp")
    # print(l_val)
    return ua_perp, l_val


def side_diffmat(vc1, t1, vc2, t2, dx, inner_g, k1, h):
    l_val = \
        math.sqrt((inner_g[vc1, t1].x - inner_g[vc2, t2].x) ** 2
                  + (inner_g[vc1, t1].y - inner_g[vc2, t2].y) ** 2
                  + (inner_g[vc1, t1].z - inner_g[vc2, t2].z) ** 2) / 2
    us_side = 1 / (l_val / k1 + 1 / h) * (0.17677 * dx ** 2)
    #    global t_count
    #    t_count = t_count + 1
    # print("side")
    return us_side, l_val


def diag_diffmat(vc1, t1, vc2, t2, dx, inner_g, k1, h):
    l_val = \
        math.sqrt((inner_g[vc1, t1].x - inner_g[vc2, t2].x) ** 2
                  + (inner_g[vc1, t1].y - inner_g[vc2, t2].y) ** 2
                  + (inner_g[vc1, t1].z - inner_g[vc2, t2].z) ** 2) / 2
    ua_diag = 1 / (l_val / k1 + 1 / h) * (0.3535 * dx ** 2)
    #    global t_count
    #    t_count = t_count + 1
    # print("diag")
    return ua_diag, l_val


Temp = steadyStateHeatTransfer.Tdomtetra

Tvox = Vox.Tdom

Ts = []
Tsv = []

q_tetra = []
q_vox = []

outer_k = steadyStateHeatTransfer.k_val
outer_htc = steadyStateHeatTransfer.HTC


def n1(voxel_db, p, nx, ny, dx, mat, k_val, htc):
    #    print("N1")
    if voxel_db[p + nx * ny, 6].mat != mat:  # S1
        ua, l_val = perp_diffmat(p, 0, p + nx * ny, 4, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 0] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 0] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 3].mat != mat:  # n2
        #        # Ts.append(Temp[p,0])
        ua, l_val = side_diffmat(p, 0, p, 1, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 0] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 0] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 5].mat != mat:  # n4
        #        # Ts.append(Temp[p,0])
        ua, l_val = side_diffmat(p, 0, p, 3, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 0] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 0] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 22].mat != mat:  # b4
        #        # Ts.append(Temp[p,0])
        ua, l_val = side_diffmat(p, 0, p, 20, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 0] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 0] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def n2(voxel_db, p, nx, ny, dx, mat, k_val, htc):
    #    print("N2")
    if voxel_db[p + nx * ny, 7].mat != mat:  # S2
        #        # Ts.append(Temp[p,1])
        ua, l_val = perp_diffmat(p, 1, p + nx * ny, 5, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 1] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 1] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 4].mat != mat:  # n3
        #        # Ts.append(Temp[p,1])
        ua, l_val = side_diffmat(p, 1, p, 2, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 1] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 1] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 2].mat != mat:  # n1
        #        # Ts.append(Temp[p,1])
        ua, l_val = side_diffmat(p, 1, p, 0, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 1] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 1] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 10].mat != mat:  # w1
        #        # Ts.append(Temp[p,1])
        ua, l_val = diag_diffmat(p, 1, p, 8, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 1] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 1] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def n3(voxel_db, p, nx, ny, dx, mat, k_val, htc):
    #    print("N3")
    if voxel_db[p + nx * ny, 8].mat != mat:  # S3
        #        # Ts.append(Temp[p,2])
        ua, l_val = perp_diffmat(p, 2, p + nx * ny, 6, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 2] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 2] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 5].mat != mat:  # n4
        #        # Ts.append(Temp[p,2])
        ua, l_val = side_diffmat(p, 2, p, 3, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 2] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 2] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 3].mat != mat:  # n2
        #        # Ts.append(Temp[p,2])
        ua, l_val = side_diffmat(p, 2, p, 1, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 2] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 2] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 18].mat != mat:  # f1
        #        # Ts.append(Temp[p,2])
        ua, l_val = diag_diffmat(p, 2, p, 16, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 2] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 2] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def n4(voxel_db, p, nx, ny, dx, mat, k_val, htc):
    #    print("N4")
    if voxel_db[p + nx * ny, 9].mat != mat:  # S4
        # # Ts.append(Temp[p,3])
        ua, l_val = perp_diffmat(p, 3, p + nx * ny, 7, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 3] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 3] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 2].mat != mat:  # n1
        # # Ts.append(Temp[p,3])
        ua, l_val = side_diffmat(p, 3, p, 0, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 3] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 3] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 4].mat != mat:  # n3
        # Ts.append(Temp[p,3])
        ua, l_val = side_diffmat(p, 3, p, 2, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 3] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 3] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 14].mat != mat:  # E1
        # Ts.append(Temp[p,3])
        ua, l_val = diag_diffmat(p, 3, p, 12, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 3] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 3] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def s1(voxel_db, p, nx, ny, dx, mat, k_val, htc):
    #    print("S1")
    if voxel_db[p - nx * ny, 2].mat != mat:  # n1
        # Ts.append(Temp[p,4])
        ua, l_val = perp_diffmat(p, 4, p - nx * ny, 0, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 4] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 4] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 9].mat != mat:  # s4
        # Ts.append(Temp[p,4])
        ua, l_val = side_diffmat(p, 4, p, 7, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 4] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 4] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 7].mat != mat:  # s2
        # Ts.append(Temp[p,4])
        ua, l_val = side_diffmat(p, 4, p, 5, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 4] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 4] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 24].mat != mat:  # b3
        # Ts.append(Temp[p,4])
        ua, l_val = diag_diffmat(p, 4, p, 22, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 4] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 4] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def s2(voxel_db, p, nx, ny, dx, mat, k_val, htc):
    #    print("S2")
    if voxel_db[p - nx * ny, 3].mat != mat:  # n2
        # Ts.append(Temp[p,5])
        ua, l_val = perp_diffmat(p, 5, p - nx * ny, 1, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 5] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 5] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 9].mat != mat:  # s4
        # Ts.append(Temp[p,5])
        ua, l_val = side_diffmat(p, 5, p, 7, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 5] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 5] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 8].mat != mat:  # s3
        # Ts.append(Temp[p,5])
        ua, l_val = side_diffmat(p, 5, p, 6, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 5] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 5] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 12].mat != mat:  # w3
        # Ts.append(Temp[p,5])
        ua, l_val = diag_diffmat(p, 5, p, 10, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 5] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 5] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def s3(voxel_db, p, nx, ny, dx, mat, k_val, htc):
    #    print("S3")
    if voxel_db[p - nx * ny, 4].mat != mat:  # n3
        # Ts.append(Temp[p,6])
        ua, l_val = perp_diffmat(p, 6, p - nx * ny, 2, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 6] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 6] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 9].mat != mat:  # s4
        # Ts.append(Temp[p,6])
        ua, l_val = side_diffmat(p, 6, p, 7, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 6] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 6] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 7].mat != mat:  # s2
        # Ts.append(Temp[p,6])
        ua, l_val = side_diffmat(p, 6, p, 5, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 6] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 6] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 20].mat != mat:  # f3
        # Ts.append(Temp[p,6])
        ua, l_val = diag_diffmat(p, 6, p, 18, dx, outer_g, k_val, htc)
        inner_q = (Temp[p, 6] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def s4(voxel_db, p, nx, ny, dx, mat, k_val, htc):
    #    print("S4")
    if voxel_db[p - nx * ny, 5].mat != mat:  # n4
        # Ts.append(Temp[p,7])
        ua, l_val = perp_diffmat(p, 7, p - nx * ny, 3, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 7] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 7] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 6].mat != mat:  # s1
        # Ts.append(Temp[p,7])
        ua, l_val = side_diffmat(p, 7, p, 4, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 7] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 7] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 8].mat != mat:  # s3
        # Ts.append(Temp[p,7])
        ua, l_val = side_diffmat(p, 7, p, 6, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 7] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 7] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 16].mat != mat:  # E3
        # Ts.append(Temp[p,7])
        ua, l_val = diag_diffmat(p, 7, p, 14, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 7] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 7] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def w1(voxel_db, p, ny, dx, mat, k_val, htc):
    #    print("W1")
    if voxel_db[p - ny, 14].mat != mat:  # E1
        # Ts.append(Temp[p,8])
        ua, l_val = perp_diffmat(p, 8, p - ny, 12, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 8] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 8] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 11].mat != mat:  # w2
        # Ts.append(Temp[p,8])
        ua, l_val = side_diffmat(p, 8, p, 9, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 8] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 8] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 13].mat != mat:  # w4
        # Ts.append(Temp[p,8])
        ua, l_val = side_diffmat(p, 8, p, 11, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 8] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 8] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 3].mat != mat:  # n2
        # Ts.append(Temp[p,8])
        ua, l_val = diag_diffmat(p, 8, p, 1, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 8] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 8] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def w2(voxel_db, p, ny, dx, mat, k_val, htc):
    #    print("W2")
    if voxel_db[p - ny, 17].mat != mat:  # E4
        # Ts.append(Temp[p,9])
        ua, l_val = perp_diffmat(p, 9, p - ny, 15, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 9] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 9] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 12].mat != mat:  # w3
        # Ts.append(Temp[p,9])
        ua, l_val = side_diffmat(p, 9, p, 10, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 9] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 9] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 10].mat != mat:  # w1
        # Ts.append(Temp[p,9])
        ua, l_val = side_diffmat(p, 9, p, 8, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 9] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 9] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 25].mat != mat:  # b4
        # Ts.append(Temp[p,9])
        ua, l_val = diag_diffmat(p, 9, p, 23, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 9] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 9] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def w3(voxel_db, p, ny, dx, mat, k_val, htc):
    #    print("W3")
    if voxel_db[p - ny, 16].mat != mat:  # E3
        # Ts.append(Temp[p,10])
        ua, l_val = perp_diffmat(p, 10, p - ny, 14, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 10] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 10] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 13].mat != mat:  # w4
        # Ts.append(Temp[p,10])
        ua, l_val = side_diffmat(p, 10, p, 11, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 10] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 10] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 11].mat != mat:  # w2
        # Ts.append(Temp[p,10])
        ua, l_val = side_diffmat(p, 10, p, 9, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 10] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 10] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 7].mat != mat:  # s2
        # Ts.append(Temp[p,10])
        ua, l_val = diag_diffmat(p, 10, p, 5, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 10] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 10] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def w4(voxel_db, p, ny, dx, mat, k_val, htc):
    #    print("W4")
    if voxel_db[p - ny, 15].mat != mat:  # E2
        # Ts.append(Temp[p,11])
        ua, l_val = perp_diffmat(p, 11, p - ny, 13, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 11] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 11] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 10].mat != mat:  # w1
        # Ts.append(Temp[p,11])
        ua, l_val = side_diffmat(p, 11, p, 8, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 11] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 11] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 12].mat != mat:  # w3
        # Ts.append(Temp[p,11])
        ua, l_val = side_diffmat(p, 11, p, 10, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 11] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 11] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 19].mat != mat:  # f2
        # Ts.append(Temp[p,11])
        ua, l_val = diag_diffmat(p, 11, p, 17, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 11] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 11] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def e1(voxel_db, p, ny, dx, mat, k_val, htc):
    #    print("E1")
    if voxel_db[p + ny, 10].mat != mat:  # w1
        # Ts.append(Temp[p,12])
        ua, l_val = perp_diffmat(p, 12, p + ny, 8, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 12] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 12] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 15].mat != mat:  # e2
        # Ts.append(Temp[p,12])
        ua, l_val = side_diffmat(p, 12, p, 13, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 12] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 12] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 17].mat != mat:  # e4
        # Ts.append(Temp[p,12])
        ua, l_val = side_diffmat(p, 12, p, 15, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 12] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 12] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 5].mat != mat:  # n4
        # Ts.append(Temp[p,12])
        ua, l_val = diag_diffmat(p, 12, p, 3, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 12] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 12] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def e2(voxel_db, p, ny, dx, mat, k_val, htc):
    #    print("E2")
    if voxel_db[p + ny, 13].mat != mat:  # w4
        # Ts.append(Temp[p,13])
        ua, l_val = perp_diffmat(p, 13, p + ny, 11, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 13] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 13] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 16].mat != mat:  # e3
        # Ts.append(Temp[p,13])
        ua, l_val = side_diffmat(p, 13, p, 14, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 13] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 13] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 14].mat != mat:  # e1
        # Ts.append(Temp[p,13])
        ua, l_val = side_diffmat(p, 13, p, 12, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 13] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 13] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 21].mat != mat:  # f4
        # Ts.append(Temp[p,13])
        ua, l_val = diag_diffmat(p, 13, p, 19, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 13] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 13] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def e3(voxel_db, p, ny, dx, mat, k_val, htc):
    #    print("E3")
    if voxel_db[p + ny, 12].mat != mat:  # w3
        # Ts.append(Temp[p,14])
        ua, l_val = perp_diffmat(p, 14, p + ny, 10, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 14] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 14] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 17].mat != mat:  # e4
        # Ts.append(Temp[p,14])
        ua, l_val = side_diffmat(p, 14, p, 15, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 14] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 14] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 15].mat != mat:  # e2
        # Ts.append(Temp[p,14])
        ua, l_val = side_diffmat(p, 14, p, 13, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 14] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 14] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 9].mat != mat:  # s4
        # Ts.append(Temp[p,14])
        ua, l_val = diag_diffmat(p, 14, p, 7, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 14] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 14] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def e4(voxel_db, p, ny, dx, mat, k_val, htc):
    #    print("E4")
    if voxel_db[p + ny, 11].mat != mat:  # w2
        # Ts.append(Temp[p,15])
        ua, l_val = perp_diffmat(p, 15, p + ny, 9, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 15] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 15] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 14].mat != mat:  # e1
        # Ts.append(Temp[p,15])
        ua, l_val = side_diffmat(p, 15, p, 12, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 15] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 15] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 16].mat != mat:  # e3
        # Ts.append(Temp[p,15])
        ua, l_val = side_diffmat(p, 15, p, 14, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 15] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 15] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 23].mat != mat:  # b2
        # Ts.append(Temp[p,15])
        ua, l_val = diag_diffmat(p, 15, p, 21, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 15] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 15] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def f1(voxel_db, p, dx, mat, k_val, htc):
    #    print("F1")
    if voxel_db[p + 1, 22].mat != mat:  # b1
        # Ts.append(Temp[p,16])
        ua, l_val = perp_diffmat(p, 16, p + 1, 20, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 16] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 16] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 19].mat != mat:  # f2
        # Ts.append(Temp[p,16])
        ua, l_val = side_diffmat(p, 16, p, 17, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 16] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 16] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 21].mat != mat:  # f4
        # Ts.append(Temp[p,16])
        ua, l_val = side_diffmat(p, 16, p, 19, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 16] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 16] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 4].mat != mat:  # n3
        # Ts.append(Temp[p,16])
        ua, l_val = diag_diffmat(p, 16, p, 2, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 16] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 16] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def f2(voxel_db, p, dx, mat, k_val, htc):
    #    print("F2")
    if voxel_db[p + 1, 25].mat != mat:  # b4
        # Ts.append(Temp[p,17])
        ua, l_val = perp_diffmat(p, 17, p + 1, 23, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 17] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 17] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 20].mat != mat:  # f3
        # Ts.append(Temp[p,17])
        ua, l_val = side_diffmat(p, 17, p, 18, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 17] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 17] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 18].mat != mat:  # f1
        # Ts.append(Temp[p,17])
        ua, l_val = side_diffmat(p, 17, p, 16, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 17] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 17] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 13].mat != mat:  # w4
        # Ts.append(Temp[p,17])
        ua, l_val = diag_diffmat(p, 17, p, 11, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 17] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 17] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def f3(voxel_db, p, dx, mat, k_val, htc):
    #    print("F3")
    if voxel_db[p + 1, 24].mat != mat:  # b3
        # Ts.append(Temp[p,18])
        ua, l_val = perp_diffmat(p, 18, p + 1, 22, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 18] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 18] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 21].mat != mat:  # f4
        # Ts.append(Temp[p,18])
        ua, l_val = side_diffmat(p, 18, p, 19, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 18] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 18] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 19].mat != mat:  # f2
        # Ts.append(Temp[p,18])
        ua, l_val = side_diffmat(p, 18, p, 17, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 18] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 18] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 8].mat != mat:  # s3
        # Ts.append(Temp[p,18])
        ua, l_val = diag_diffmat(p, 18, p, 6, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 18] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 18] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def f4(voxel_db, p, dx, mat, k_val, htc):
    #    print("F4")
    if voxel_db[p + 1, 23].mat != mat:  # b2
        # Ts.append(Temp[p,19])
        ua, l_val = perp_diffmat(p, 19, p + 1, 21, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 19] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 19] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 18].mat != mat:  # f1
        # Ts.append(Temp[p,19])
        ua, l_val = side_diffmat(p, 19, p, 16, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 19] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 19] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 20].mat != mat:  # f3
        # Ts.append(Temp[p,19])
        ua, l_val = side_diffmat(p, 19, p, 18, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 19] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 19] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 15].mat != mat:  # e2
        # Ts.append(Temp[p,19])
        ua, l_val = diag_diffmat(p, 19, p, 13, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 19] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 19] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def b1(voxel_db, p, dx, mat, k_val, htc):
    #    print("B1")
    if voxel_db[p - 1, 18].mat != mat:  # f1
        # Ts.append(Temp[p,20])
        ua, l_val = perp_diffmat(p, 20, p - 1, 16, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 20] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 20] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 23].mat != mat:  # b2
        # Ts.append(Temp[p,20])
        ua, l_val = side_diffmat(p, 20, p, 21, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 20] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 20] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 25].mat != mat:  # b4
        # Ts.append(Temp[p,20])
        ua, l_val = side_diffmat(p, 20, p, 23, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 20] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 20] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 2].mat != mat:  # n1
        # Ts.append(Temp[p,20])
        ua, l_val = diag_diffmat(p, 20, p, 0, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 20] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 20] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def b2(voxel_db, p, dx, mat, k_val, htc):
    #    print("B2")
    if voxel_db[p - 1, 21].mat != mat:  # f4
        # Ts.append(Temp[p,21])
        ua, l_val = perp_diffmat(p, 21, p - 1, 19, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 21] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 21] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 24].mat != mat:  # b3
        # Ts.append(Temp[p,21])
        ua, l_val = side_diffmat(p, 21, p, 22, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 21] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 21] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 22].mat != mat:  # b1
        # Ts.append(Temp[p,21])
        ua, l_val = side_diffmat(p, 21, p, 20, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 21] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 21] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 17].mat != mat:  # E4
        # Ts.append(Temp[p,21])
        ua, l_val = diag_diffmat(p, 21, p, 15, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 21] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 21] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def b3(voxel_db, p, dx, mat, k_val, htc):
    #    print("B3")
    if voxel_db[p - 1, 20].mat != mat:  # f3
        # Ts.append(Temp[p,22])
        ua, l_val = perp_diffmat(p, 22, p - 1, 18, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 22] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 22] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 25].mat != mat:  # b4
        # Ts.append(Temp[p,22])
        ua, l_val = side_diffmat(p, 22, p, 23, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 22] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 22] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 23].mat != mat:  # b2
        # Ts.append(Temp[p,22])
        ua, l_val = side_diffmat(p, 22, p, 21, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 22] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 22] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 6].mat != mat:  # s3
        # Ts.append(Temp[p,22])
        ua, l_val = diag_diffmat(p, 22, p, 4, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 22] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 22] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def b4(voxel_db, p, dx, mat, k_val, htc):
    #    print("B4")
    if voxel_db[p - 1, 19].mat != mat:  # f2
        # Ts.append(Temp[p,23])
        ua, l_val = perp_diffmat(p, 23, p - 1, 17, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 23] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 23] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 22].mat != mat:  # b1
        # Ts.append(Temp[p,23])
        ua, l_val = side_diffmat(p, 23, p, 20, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 23] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 23] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 24].mat != mat:  # b3
        # Ts.append(Temp[p,23])
        ua, l_val = side_diffmat(p, 23, p, 22, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 23] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 23] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)
    if voxel_db[p, 11].mat != mat:  # w2
        # Ts.append(Temp[p,23])
        ua, l_val = diag_diffmat(p, 23, p, 9, dx, outer_g, k_val, htc)
        tt_val = (Temp[p, 23] + htc * l_val / k * steadyStateHeatTransfer.Tamb) / (1 + htc * l_val / k)
        Ts.append(tt_val)
        inner_q = (Temp[p, 23] - steadyStateHeatTransfer.Tamb) * ua
        q_tetra.append(inner_q)


def func(tag, voxel_db, p, nx, ny, dx, mat, k_val, htc):
    if tag == 'N1':
        n1(voxel_db, p, nx, ny, dx, mat, k_val, htc)
    elif tag == 'N2':
        n2(voxel_db, p, nx, ny, dx, mat, k_val, htc)
    elif tag == 'N3':
        n3(voxel_db, p, nx, ny, dx, mat, k_val, htc)
    elif tag == 'N4':
        n4(voxel_db, p, nx, ny, dx, mat, k_val, htc)

    elif tag == 'S1':
        s1(voxel_db, p, nx, ny, dx, mat, k_val, htc)
    elif tag == 'S2':
        s2(voxel_db, p, nx, ny, dx, mat, k_val, htc)
    elif tag == 'S3':
        s3(voxel_db, p, nx, ny, dx, mat, k_val, htc)
    elif tag == 'S4':
        s4(voxel_db, p, nx, ny, dx, mat, k_val, htc)

    elif tag == 'W1':
        w1(voxel_db, p, ny, dx, mat, k_val, htc)
    elif tag == 'W2':
        w2(voxel_db, p, ny, dx, mat, k_val, htc)
    elif tag == 'W3':
        w3(voxel_db, p, ny, dx, mat, k_val, htc)
    elif tag == 'W4':
        w4(voxel_db, p, ny, dx, mat, k_val, htc)

    elif tag == 'E1':
        e1(voxel_db, p, ny, dx, mat, k_val, htc)
    elif tag == 'E2':
        e2(voxel_db, p, ny, dx, mat, k_val, htc)
    elif tag == 'E3':
        e3(voxel_db, p, ny, dx, mat, k_val, htc)
    elif tag == 'E4':
        e4(voxel_db, p, ny, dx, mat, k_val, htc)

    elif tag == 'F1':
        f1(voxel_db, p, dx, mat, k_val, htc)
    elif tag == 'F2':
        f2(voxel_db, p, dx, mat, k_val, htc)
    elif tag == 'F3':
        f3(voxel_db, p, dx, mat, k_val, htc)
    elif tag == 'F4':
        f4(voxel_db, p, dx, mat, k_val, htc)

    elif tag == 'B1':
        b1(voxel_db, p, dx, mat, k_val, htc)
    elif tag == 'B2':
        b2(voxel_db, p, dx, mat, k_val, htc)
    elif tag == 'B3':
        b3(voxel_db, p, dx, mat, k_val, htc)
    elif tag == 'B4':
        b4(voxel_db, p, dx, mat, k_val, htc)


for vc in range(Ellipse.voxel_n):
    for i in range(24):
        if Ellipse.voxel_db[vc, i + 2].mat != 0:
            func(Ellipse.voxel_db[vc, i + 2].pos, Ellipse.voxel_db, vc, Ellipse.nx, Ellipse.ny, Ellipse.dx,
                 Ellipse.voxel_db[vc, i + 2].mat, steadyStateHeatTransfer.k_val, steadyStateHeatTransfer.HTC)


Tsa = steadyStateHeatTransfer.Heat_rate * Ellipse.a / (3 * steadyStateHeatTransfer.HTC) + steadyStateHeatTransfer.Tamb
# pyplot.plot(Ts)
# pyplot.axhline(Tsa,color='red')
# pyplot.show()

error = []
for i in range(len(Ts)):
    error.append(Tsa - Ts[i])

# pyplot.boxplot(error)
# pyplot.show()


# voxelplot = 1
# if(voxelplot == 1):


errv = []
for k in range(Vox.n):
    for j in range(Vox.n):
        for i in range(Vox.n):
            if Vox.dom[i, j, k] == 1:
                if Vox.side_exposed[i, j, k] > 0:
                    tt_valv = \
                        (Vox.Tdom[i, j, k]
                         + (steadyStateHeatTransfer.HTC
                            * (Ellipse.dx / 2) / steadyStateHeatTransfer.k_val)
                            * steadyStateHeatTransfer.Tamb) / (1 + steadyStateHeatTransfer.HTC
                                                               * (Ellipse.dx / 2) / steadyStateHeatTransfer.k_val)
                    Tsv.append(tt_valv)
                    errv.append(Tsa - tt_valv)
                    outer_ua = \
                        1 / (0.5 * Ellipse.dx / steadyStateHeatTransfer.k_val + 1 / steadyStateHeatTransfer.HTC) \
                        * Ellipse.dx * Ellipse.dx
                    outer_q = Vox.side_exposed[i, j, k] * outer_ua * (Vox.Tdom[i, j, k] - steadyStateHeatTransfer.Tamb)
                    q_vox.append(outer_q)

data = [error, errv]
labels = ["Tetrahedrons", "Voxels"]

pyplot.figure(figsize=(10, 6))
pyplot.subplot(121)
bplot = pyplot.boxplot(data, vert=True, patch_artist=True, labels=labels)
colors = ['pink', 'lightblue']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
    pyplot.grid(color='grey')
    pyplot.ylabel("Surface Temperature Error in degrees Celsius")

x_t = []
x_v = []
t = 1
for i in range(len(Ts)):
    x_t.append((4 * math.pi * 0.1 ** 2) / t)
    t = t + 1
t = 1
for i in range(len(Tsv)):
    x_v.append((4 * math.pi * 0.1 ** 2) / t)
    t = t + 1

# pyplot.figure(figsize = (10,6))
pyplot.subplot(122)
# pyplot.plot(x_t,Ts,label = "Tetra")
# pyplot.plot(x_v,Tsv,label = "Voxel")
pyplot.axhline(numpy.mean(Ts), label="Mean Tetra", color="g_valreen")
pyplot.axhline(numpy.mean(Tsv), label="Mean Voxel", color="Purple")
pyplot.axhline(Tsa, color="red", label="Analytical")
pyplot.yticks(numpy.arange(steadyStateHeatTransfer.Tamb - 2, Tsa + 2, 2))
pyplot.ylabel("Surface Temperature")
pyplot.legend()
# pyplot.savefig("Results_1mm_R10cm.png")
pyplot.show()
pyplot.show()

print("\nAnalytical Ts =", Tsa)
print("\nMean Ts tetra", numpy.mean(Ts))
print("Mean Ts Voxel", numpy.mean(Tsv), "\n")

print("Tetrahedron error =", steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume - numpy.sum(q_tetra))
print("Heart rate = ", steadyStateHeatTransfer.Heat_rate)
print("tetra_volume = ", Ellipse.tetra_volume)
print("percentage error = ",
      (steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume - numpy.sum(q_tetra)) / (
                  steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume) * 100, "%")
print("Tetra T_surf for given", steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume / (
            Ellipse.area_sum * steadyStateHeatTransfer.HTC) + steadyStateHeatTransfer.Tamb)

print("voxel error =", steadyStateHeatTransfer.Heat_rate * Ellipse.voxel_volume - numpy.sum(q_vox), "W")
print("Voxel heat error percentage =",
      (steadyStateHeatTransfer.Heat_rate * Ellipse.voxel_volume - numpy.sum(
          q_vox)) / steadyStateHeatTransfer.Heat_rate * Ellipse.voxel_volume * 100, "%")

print("\nActual tetra Q = ", steadyStateHeatTransfer.Heat_rate * Ellipse.tetra_volume)
print("Tetra Q =  ", numpy.sum(q_tetra))
print("voxel Q =  ", numpy.sum(q_vox))

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
