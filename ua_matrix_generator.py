import faulthandler
import math as mt
import sys

import numpy as np  # from numpy import linalg as la
import scipy.sparse as sp

faulthandler.enable()
#####################
deci = 5


def perp_samemat(vc1, t1, vc2, t2, dx, g_val, k1, k2):
    l_val = round(mt.sqrt((g_val[vc1, t1].x - g_val[vc2, t2].x) ** 2 + (g_val[vc1, t1].y - g_val[vc2, t2].y) ** 2 + (
            g_val[vc1, t1].z - g_val[vc2, t2].z) ** 2) / 2, deci)
    ua_perp = 1 / (l_val / k1 + l_val / k2) * (0.25 * dx ** 2)
    # print("perp distance", l_val)
    #    print("U same perp",ua_perp/(0.25*dx**2))
    return ua_perp


def perp_diffmat(vc1, t1, vc2, t2, dx, g_val, k1, h):
    l_val = round(mt.sqrt((g_val[vc1, t1].x - g_val[vc2, t2].x) ** 2 + (g_val[vc1, t1].y - g_val[vc2, t2].y) ** 2 + (
            g_val[vc1, t1].z - g_val[vc2, t2].z) ** 2) / 2, deci)
    ua_perp = 1 / (l_val / k1 + 1 / h) * (0.25 * dx ** 2)
    # print("perp distance", l_val)
    #    print("U diff perp",ua_perp/(0.25*dx**2))
    return ua_perp


def side_samemat(vc1, t1, vc2, t2, dx, g_val, k1, k2):
    l_val = round(mt.sqrt((g_val[vc1, t1].x - g_val[vc2, t2].x) ** 2 + (g_val[vc1, t1].y - g_val[vc2, t2].y) ** 2 + (
            g_val[vc1, t1].z - g_val[vc2, t2].z) ** 2) / 2, deci)
    ua_side = 1 / (l_val / k1 + l_val / k2) * (0.17677 * dx ** 2)
    # print("side distance", l_val)
    #    print("U same side",ua_side/(0.17677*dx**2))
    return ua_side


def side_diffmat(vc1, t1, vc2, t2, dx, g_val, k1, h):
    l_val = round(mt.sqrt((g_val[vc1, t1].x - g_val[vc2, t2].x) ** 2 + (g_val[vc1, t1].y - g_val[vc2, t2].y) ** 2 + (
            g_val[vc1, t1].z - g_val[vc2, t2].z) ** 2) / 2, deci)
    ua_side = 1 / (l_val / k1 + 1 / h) * (0.17677 * dx ** 2)
    # print("side distance", l_val)
    #    print("U diff side",ua_side/(0.17677*dx**2))
    return ua_side


def diag_samemat(vc1, t1, vc2, t2, dx, g_val, k1, k2):
    l_val = round(mt.sqrt((g_val[vc1, t1].x - g_val[vc2, t2].x) ** 2 + (g_val[vc1, t1].y - g_val[vc2, t2].y) ** 2 + (
            g_val[vc1, t1].z - g_val[vc2, t2].z) ** 2) / 2, deci)
    ua_diag = 1 / (l_val / k1 + l_val / k2) * (0.3535 * dx ** 2)
    # print("diag distance", l_val)
    #    print("U same diag",ua_diag/(0.3535*dx**2))
    return ua_diag


def diag_diffmat(vc1, t1, vc2, t2, dx, g_val, k1, h):
    l_val = round(mt.sqrt((g_val[vc1, t1].x - g_val[vc2, t2].x) ** 2 + (g_val[vc1, t1].y - g_val[vc2, t2].y) ** 2 + (
            g_val[vc1, t1].z - g_val[vc2, t2].z) ** 2) / 2, deci)
    ua_diag = 1 / (l_val / k1 + 1 / h) * (0.3535 * dx ** 2)
    # print("diag distance", l_val)
    #    print("U diff diag",ua_diag/(0.3535*dx**2))
    return ua_diag


# north Tetras ####################################################

def n1(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("n1")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + nx * ny, 6].mat == 0:  # s1

        ua_perp = perp_diffmat(p, 0, p + nx * ny, 4, dx, g_val, k_val[p, 0], htc)

    if voxel_db[p + nx * ny, 6].mat == mat:
        ua_perp = perp_samemat(p, 0, p + nx * ny, 4, dx, g_val, k_val[p, 0], k_val[p + nx * ny, 4])

    if voxel_db[p, 3].mat == 0:  # n2

        ua_s1 = side_diffmat(p, 0, p, 1, dx, g_val, k_val[p, 0], htc)

    if voxel_db[p, 3].mat == mat:
        ua_s1 = side_samemat(p, 0, p, 1, dx, g_val, k_val[p, 0], k_val[p, 1])

    if voxel_db[p, 5].mat == 0:  # n4

        ua_s2 = side_diffmat(p, 0, p, 3, dx, g_val, k_val[p, 0], htc)

    if voxel_db[p, 5].mat == mat:
        ua_s2 = side_samemat(p, 0, p, 3, dx, g_val, k_val[p, 0], k_val[p, 3])

    if voxel_db[p, 22].mat == 0:  # b1

        ua_diag = diag_diffmat(p, 0, p, 20, dx, g_val, k_val[p, 0], htc)

    if voxel_db[p, 22].mat == mat:
        ua_diag = diag_samemat(p, 0, p, 20, dx, g_val, k_val[p, 0], k_val[p, 20])

    i_cor.append(p * 24 + 0)
    j_cor.append(((p + nx * ny) * 24 + 4))
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 0)
    j_cor.append(p * 24 + 1)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 0)
    j_cor.append(p * 24 + 3)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 0)
    j_cor.append(p * 24 + 20)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 0)
    j_cor.append(p * 24 + 0)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def n2(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("n2")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + nx * ny, 7].mat == 0:  # s2

        ua_perp = perp_diffmat(p, 1, p + nx * ny, 5, dx, g_val, k_val[p, 1], htc)

    if voxel_db[p + nx * ny, 7].mat == mat:
        ua_perp = perp_samemat(p, 1, p + nx * ny, 5, dx, g_val, k_val[p, 1], k_val[p + nx * ny, 5])

    if voxel_db[p, 4].mat == 0:  # n3

        ua_s1 = side_diffmat(p, 1, p, 2, dx, g_val, k_val[p, 1], htc)

    if voxel_db[p, 4].mat == mat:
        ua_s1 = side_samemat(p, 1, p, 2, dx, g_val, k_val[p, 1], k_val[p, 2])

    if voxel_db[p, 2].mat == 0:  # n1

        ua_s2 = side_diffmat(p, 1, p, 0, dx, g_val, k_val[p, 1], htc)

    if voxel_db[p, 2].mat == mat:
        ua_s2 = side_samemat(p, 1, p, 0, dx, g_val, k_val[p, 1], k_val[p, 0])

    if voxel_db[p, 10].mat == 0:  # w1

        ua_diag = diag_diffmat(p, 1, p, 8, dx, g_val, k_val[p, 1], htc)

    if voxel_db[p, 10].mat == mat:
        ua_diag = diag_samemat(p, 1, p, 8, dx, g_val, k_val[p, 1], k_val[p, 8])

    i_cor.append(p * 24 + 1)
    j_cor.append((p + nx * ny) * 24 + 5)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 1)
    j_cor.append(p * 24 + 2)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 1)
    j_cor.append(p * 24 + 0)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 1)
    j_cor.append(p * 24 + 8)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 1)
    j_cor.append(p * 24 + 1)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def n3(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("n3")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + nx * ny, 8].mat == 0:  # s3

        ua_perp = perp_diffmat(p, 2, (p + nx * ny), 6, dx, g_val, k_val[p, 2], htc)

    if voxel_db[p + nx * ny, 8].mat == mat:
        ua_perp = perp_samemat(p, 2, (p + nx * ny), 6, dx, g_val, k_val[p, 2], k_val[p + nx * ny, 6])

    if voxel_db[p, 5].mat == 0:  # n4

        ua_s1 = side_diffmat(p, 2, p, 3, dx, g_val, k_val[p, 2], htc)

    if voxel_db[p, 5].mat == mat:
        ua_s1 = side_samemat(p, 2, p, 3, dx, g_val, k_val[p, 2], k_val[p, 3])

    if voxel_db[p, 3].mat == 0:  # n2

        ua_s2 = side_diffmat(p, 2, p, 1, dx, g_val, k_val[p, 2], htc)

    if voxel_db[p, 3].mat == mat:
        ua_s2 = side_samemat(p, 2, p, 1, dx, g_val, k_val[p, 2], k_val[p, 1])

    if voxel_db[p, 18].mat == 0:  # f1

        ua_diag = diag_diffmat(p, 2, p, 16, dx, g_val, k_val[p, 2], htc)

    if voxel_db[p, 18].mat == mat:
        ua_diag = diag_samemat(p, 2, p, 16, dx, g_val, k_val[p, 2], k_val[p, 16])

    i_cor.append(p * 24 + 2)
    j_cor.append((p + nx * ny) * 24 + 6)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 2)
    j_cor.append(p * 24 + 3)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 2)
    j_cor.append(p * 24 + 1)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 2)
    j_cor.append(p * 24 + 16)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 2)
    j_cor.append(p * 24 + 2)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def n4(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("n4")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + nx * ny, 9].mat == 0:  # s4

        ua_perp = perp_diffmat(p, 3, p + nx * ny, 7, dx, g_val, k_val[p, 3], htc)

    if voxel_db[p + nx * ny, 9].mat == mat:
        ua_perp = perp_samemat(p, 3, p + nx * ny, 7, dx, g_val, k_val[p, 3], k_val[p + nx * ny, 7])

    if voxel_db[p, 2].mat == 0:  # n1

        ua_s1 = side_diffmat(p, 3, p, 0, dx, g_val, k_val[p, 3], htc)

    if voxel_db[p, 2].mat == mat:
        ua_s1 = side_samemat(p, 3, p, 0, dx, g_val, k_val[p, 3], k_val[p, 0])

    if voxel_db[p, 4].mat == 0:  # n3

        ua_s2 = side_diffmat(p, 3, p, 2, dx, g_val, k_val[p, 3], htc)

    if voxel_db[p, 4].mat == mat:
        ua_s2 = side_samemat(p, 3, p, 2, dx, g_val, k_val[p, 3], k_val[p, 2])

    if voxel_db[p, 14].mat == 0:  # e1

        ua_diag = diag_diffmat(p, 3, p, 12, dx, g_val, k_val[p, 3], htc)

    if voxel_db[p, 14].mat == mat:
        ua_diag = diag_samemat(p, 3, p, 12, dx, g_val, k_val[p, 3], k_val[p, 12])

    i_cor.append(p * 24 + 3)
    j_cor.append((p + nx * ny) * 24 + 7)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 3)
    j_cor.append(p * 24 + 0)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 3)
    j_cor.append(p * 24 + 2)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 3)
    j_cor.append(p * 24 + 12)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 3)
    j_cor.append(p * 24 + 3)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


# south Tetras #######################################################

def s1(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("s1")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - nx * ny, 2].mat == 0:  # n1

        ua_perp = perp_diffmat(p, 4, p - nx * ny, 0, dx, g_val, k_val[p, 4], htc)

    if voxel_db[p - nx * ny, 2].mat == mat:  # n1

        ua_perp = perp_samemat(p, 4, p - nx * ny, 0, dx, g_val, k_val[p, 4], k_val[p - nx * ny, 0])

    if voxel_db[p, 9].mat == 0:  # s4

        ua_s1 = side_diffmat(p, 4, p, 7, dx, g_val, k_val[p, 4], htc)

    if voxel_db[p, 9].mat == mat:  # s4

        ua_s1 = side_samemat(p, 4, p, 7, dx, g_val, k_val[p, 4], k_val[p, 7])

    if voxel_db[p, 7].mat == 0:  # s2

        ua_s2 = side_diffmat(p, 4, p, 5, dx, g_val, k_val[p, 4], htc)

    if voxel_db[p, 7].mat == mat:  # s2

        ua_s2 = side_samemat(p, 4, p, 5, dx, g_val, k_val[p, 4], k_val[p, 5])

    if voxel_db[p, 24].mat == 0:  # b3

        ua_diag = diag_diffmat(p, 4, p, 22, dx, g_val, k_val[p, 4], htc)

    if voxel_db[p, 24].mat == mat:  # b3

        ua_diag = diag_samemat(p, 4, p, 22, dx, g_val, k_val[p, 4], k_val[p, 22])

    i_cor.append(p * 24 + 4)
    j_cor.append((p - nx * ny) * 24 + 0)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 4)
    j_cor.append(p * 24 + 7)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 4)
    j_cor.append(p * 24 + 5)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 4)
    j_cor.append(p * 24 + 22)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 4)
    j_cor.append(p * 24 + 4)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def s2(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("s2")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - nx * ny, 3].mat != mat:  # n2

        ua_perp = perp_diffmat(p, 5, p - nx * ny, 1, dx, g_val, k_val[p, 5], htc)

    if voxel_db[p - nx * ny, 3].mat == mat:  # n2

        ua_perp = perp_samemat(p, 5, p - nx * ny, 1, dx, g_val, k_val[p, 5], k_val[p - nx * ny, 1])

    if voxel_db[p, 6].mat != mat:  # s1

        ua_s1 = side_diffmat(p, 5, p, 4, dx, g_val, k_val[p, 5], htc)

    if voxel_db[p, 6].mat == mat:  # s1

        ua_s1 = side_samemat(p, 5, p, 4, dx, g_val, k_val[p, 5], k_val[p, 4])

    if voxel_db[p, 8].mat != mat:  # s3

        ua_s2 = side_diffmat(p, 5, p, 6, dx, g_val, k_val[p, 5], htc)

    if voxel_db[p, 8].mat == mat:  # s3

        ua_s2 = side_samemat(p, 5, p, 6, dx, g_val, k_val[p, 5], k_val[p, 6])

    if voxel_db[p, 12].mat != mat:  # w3

        ua_diag = diag_diffmat(p, 5, p, 10, dx, g_val, k_val[p, 5], htc)

    if voxel_db[p, 12].mat == mat:  # w3

        ua_diag = diag_samemat(p, 5, p, 10, dx, g_val, k_val[p, 5], k_val[p, 10])

    i_cor.append(p * 24 + 5)
    j_cor.append((p - nx * ny) * 24 + 1)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 5)
    j_cor.append(p * 24 + 4)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 5)
    j_cor.append(p * 24 + 6)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 5)
    j_cor.append(p * 24 + 10)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 5)
    j_cor.append(p * 24 + 5)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def s3(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("s3")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - nx * ny, 4].mat != mat:  # n3

        ua_perp = perp_diffmat(p, 6, p - nx * ny, 2, dx, g_val, k_val[p, 6], htc)

    if voxel_db[p - nx * ny, 4].mat == mat:  # n3

        ua_perp = perp_samemat(p, 6, p - nx * ny, 2, dx, g_val, k_val[p, 6], k_val[p - nx * ny, 2])

    if voxel_db[p, 9].mat != mat:  # s4

        ua_s1 = side_diffmat(p, 6, p, 7, dx, g_val, k_val[p, 6], htc)

    if voxel_db[p, 9].mat == mat:  # s4

        ua_s1 = side_samemat(p, 6, p, 7, dx, g_val, k_val[p, 6], k_val[p, 7])

    if voxel_db[p, 7].mat != mat:  # s2

        ua_s2 = side_diffmat(p, 6, p, 5, dx, g_val, k_val[p, 6], htc)

    if voxel_db[p, 7].mat == mat:  # s2

        ua_s2 = side_samemat(p, 6, p, 5, dx, g_val, k_val[p, 6], k_val[p, 5])

    if voxel_db[p, 20].mat != mat:  # f3

        ua_diag = diag_diffmat(p, 6, p, 18, dx, g_val, k_val[p, 6], htc)

    if voxel_db[p, 20].mat == mat:  # f3

        ua_diag = diag_samemat(p, 6, p, 18, dx, g_val, k_val[p, 6], k_val[p, 18])

    i_cor.append(p * 24 + 6)
    j_cor.append((p - nx * ny) * 24 + 2)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 6)
    j_cor.append(p * 24 + 7)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 6)
    j_cor.append(p * 24 + 5)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 6)
    j_cor.append(p * 24 + 18)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 6)
    j_cor.append(p * 24 + 6)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def s4(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("s4")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - nx * ny, 5].mat != mat:  # n4

        ua_perp = perp_diffmat(p, 7, p - nx * ny, 3, dx, g_val, k_val[p, 7], htc)

    if voxel_db[p - nx * ny, 5].mat == mat:  # n4

        ua_perp = perp_samemat(p, 7, p - nx * ny, 3, dx, g_val, k_val[p, 7], k_val[p - nx * ny, 3])

    if voxel_db[p, 6].mat != mat:  # s1

        ua_s1 = side_diffmat(p, 7, p, 4, dx, g_val, k_val[p, 7], htc)

    if voxel_db[p, 6].mat == mat:  # s1

        ua_s1 = side_samemat(p, 7, p, 4, dx, g_val, k_val[p, 7], k_val[p, 4])

    if voxel_db[p, 8].mat != mat:  # s3

        ua_s2 = side_diffmat(p, 7, p, 6, dx, g_val, k_val[p, 7], htc)

    if voxel_db[p, 8].mat == mat:  # s3

        ua_s2 = side_samemat(p, 7, p, 6, dx, g_val, k_val[p, 7], k_val[p, 6])

    if voxel_db[p, 16].mat != mat:  # e3

        ua_diag = diag_diffmat(p, 7, p, 14, dx, g_val, k_val[p, 7], htc)

    if voxel_db[p, 16].mat == mat:  # e3

        ua_diag = diag_samemat(p, 7, p, 14, dx, g_val, k_val[p, 7], k_val[p, 14])

    i_cor.append(p * 24 + 7)
    j_cor.append((p - nx * ny) * 24 + 3)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 7)
    j_cor.append(p * 24 + 4)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 7)
    j_cor.append(p * 24 + 6)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 7)
    j_cor.append(p * 24 + 14)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 7)
    j_cor.append(p * 24 + 7)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


# west Tetras #############################################################

def w1(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("w1")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - ny, 14].mat != mat:  # e1

        ua_perp = perp_diffmat(p, 8, p - ny, 12, dx, g_val, k_val[p, 8], htc)

    if voxel_db[p - ny, 14].mat == mat:  # e1

        ua_perp = perp_samemat(p, 8, p - ny, 12, dx, g_val, k_val[p, 8], k_val[p - ny, 12])

    if voxel_db[p, 11].mat != mat:  # w2

        ua_s1 = side_diffmat(p, 8, p, 9, dx, g_val, k_val[p, 8], htc)

    if voxel_db[p, 11].mat == mat:  # w2

        ua_s1 = side_samemat(p, 8, p, 9, dx, g_val, k_val[p, 8], k_val[p, 9])

    if voxel_db[p, 13].mat != mat:  # w4

        ua_s2 = side_diffmat(p, 8, p, 11, dx, g_val, k_val[p, 8], htc)

    if voxel_db[p, 13].mat == mat:  # w4

        ua_s2 = side_samemat(p, 8, p, 11, dx, g_val, k_val[p, 8], k_val[p, 11])

    if voxel_db[p, 3].mat != mat:  # n2

        ua_diag = diag_diffmat(p, 8, p, 1, dx, g_val, k_val[p, 8], htc)

    if voxel_db[p, 3].mat == mat:  # n2

        ua_diag = diag_samemat(p, 8, p, 1, dx, g_val, k_val[p, 8], k_val[p, 1])

    i_cor.append(p * 24 + 8)
    j_cor.append((p - ny) * 24 + 12)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 8)
    j_cor.append(p * 24 + 9)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 8)
    j_cor.append(p * 24 + 11)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 8)
    j_cor.append(p * 24 + 1)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 8)
    j_cor.append(p * 24 + 8)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def w2(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("w2")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - ny, 17].mat != mat:  # e4

        ua_perp = perp_diffmat(p, 9, p - ny, 15, dx, g_val, k_val[p, 9], htc)

    if voxel_db[p - ny, 17].mat == mat:  # e4

        ua_perp = perp_samemat(p, 9, p - ny, 15, dx, g_val, k_val[p, 9], k_val[p - ny, 15])

    if voxel_db[p, 12].mat != mat:  # w3

        ua_s1 = side_diffmat(p, 9, p, 10, dx, g_val, k_val[p, 9], htc)

    if voxel_db[p, 12].mat == mat:  # w3

        ua_s1 = side_samemat(p, 9, p, 10, dx, g_val, k_val[p, 9], k_val[p, 10])

    if voxel_db[p, 10].mat != mat:  # w1

        ua_s2 = side_diffmat(p, 9, p, 8, dx, g_val, k_val[p, 9], htc)

    if voxel_db[p, 10].mat == mat:  # w1

        ua_s2 = side_samemat(p, 9, p, 8, dx, g_val, k_val[p, 9], k_val[p, 8])

    if voxel_db[p, 25].mat != mat:  # b4

        ua_diag = diag_diffmat(p, 9, p, 23, dx, g_val, k_val[p, 9], htc)

    if voxel_db[p, 25].mat == mat:  # b4

        ua_diag = diag_samemat(p, 9, p, 23, dx, g_val, k_val[p, 9], k_val[p, 23])

    i_cor.append(p * 24 + 9)
    j_cor.append((p - ny) * 24 + 15)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 9)
    j_cor.append(p * 24 + 10)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 9)
    j_cor.append(p * 24 + 8)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 9)
    j_cor.append(p * 24 + 23)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 9)
    j_cor.append(p * 24 + 9)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def w3(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("w3")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - ny, 16].mat != mat:  # e3

        ua_perp = perp_diffmat(p, 10, p - ny, 14, dx, g_val, k_val[p, 10], htc)

    if voxel_db[p - ny, 16].mat == mat:  # e3

        ua_perp = perp_samemat(p, 10, p - ny, 14, dx, g_val, k_val[p, 10], k_val[p - ny, 14])

    if voxel_db[p, 13].mat != mat:  # w4

        ua_s1 = side_diffmat(p, 10, p, 11, dx, g_val, k_val[p, 10], htc)

    if voxel_db[p, 13].mat == mat:  # w4

        ua_s1 = side_samemat(p, 10, p, 11, dx, g_val, k_val[p, 10], k_val[p, 11])

    if voxel_db[p, 11].mat != mat:  # w2

        ua_s2 = side_diffmat(p, 10, p, 9, dx, g_val, k_val[p, 10], htc)

    if voxel_db[p, 11].mat == mat:  # w2

        ua_s2 = side_samemat(p, 10, p, 9, dx, g_val, k_val[p, 10], k_val[p, 9])

    if voxel_db[p, 7].mat != mat:  # s2

        ua_diag = diag_diffmat(p, 10, p, 5, dx, g_val, k_val[p, 10], htc)

    if voxel_db[p, 7].mat == mat:  # s2

        ua_diag = diag_samemat(p, 10, p, 5, dx, g_val, k_val[p, 10], k_val[p, 5])

    i_cor.append(p * 24 + 10)
    j_cor.append((p - ny) * 24 + 14)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 10)
    j_cor.append(p * 24 + 11)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 10)
    j_cor.append(p * 24 + 9)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 10)
    j_cor.append(p * 24 + 5)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 10)
    j_cor.append(p * 24 + 10)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def w4(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("w4")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - ny, 15].mat != mat:  # e2

        ua_perp = perp_diffmat(p, 11, p - ny, 13, dx, g_val, k_val[p, 11], htc)

    if voxel_db[p - ny, 15].mat == mat:  # e2

        ua_perp = perp_samemat(p, 11, p - ny, 13, dx, g_val, k_val[p, 11], k_val[p - ny, 13])

    if voxel_db[p, 10].mat != mat:  # w1

        ua_s1 = side_diffmat(p, 11, p, 8, dx, g_val, k_val[p, 11], htc)

    if voxel_db[p, 10].mat == mat:  # w1

        ua_s1 = side_samemat(p, 11, p, 8, dx, g_val, k_val[p, 11], k_val[p, 8])

    if voxel_db[p, 12].mat != mat:  # w3

        ua_s2 = side_diffmat(p, 11, p, 10, dx, g_val, k_val[p, 11], htc)

    if voxel_db[p, 12].mat == mat:  # w3

        ua_s2 = side_samemat(p, 11, p, 10, dx, g_val, k_val[p, 11], k_val[p, 10])

    if voxel_db[p, 19].mat != mat:  # f2

        ua_diag = diag_diffmat(p, 11, p, 17, dx, g_val, k_val[p, 11], htc)

    if voxel_db[p, 19].mat == mat:  # f2

        ua_diag = diag_samemat(p, 11, p, 17, dx, g_val, k_val[p, 11], k_val[p, 17])

    i_cor.append(p * 24 + 11)
    j_cor.append((p - ny) * 24 + 13)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 11)
    j_cor.append(p * 24 + 8)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 11)
    j_cor.append(p * 24 + 10)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 11)
    j_cor.append(p * 24 + 17)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 11)
    j_cor.append(p * 24 + 11)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


# east Tetras ###########################################################

def e1(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("e1")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + ny, 10].mat != mat:  # w1

        ua_perp = perp_diffmat(p, 12, p + ny, 8, dx, g_val, k_val[p, 12], htc)

    if voxel_db[p + ny, 10].mat == mat:  # w1

        ua_perp = perp_samemat(p, 12, p + ny, 8, dx, g_val, k_val[p, 12], k_val[p + ny, 8])

    if voxel_db[p, 15].mat != mat:  # e2

        ua_s1 = side_diffmat(p, 12, p, 13, dx, g_val, k_val[p, 12], htc)

    if voxel_db[p, 15].mat == mat:  # e2

        ua_s1 = side_samemat(p, 12, p, 13, dx, g_val, k_val[p, 12], k_val[p, 13])

    if voxel_db[p, 17].mat != mat:  # e4

        ua_s2 = side_diffmat(p, 12, p, 15, dx, g_val, k_val[p, 12], htc)

    if voxel_db[p, 17].mat == mat:  # e4

        ua_s2 = side_samemat(p, 12, p, 15, dx, g_val, k_val[p, 12], k_val[p, 15])

    if voxel_db[p, 5].mat != mat:  # n4

        ua_diag = diag_diffmat(p, 12, p, 3, dx, g_val, k_val[p, 12], htc)

    if voxel_db[p, 5].mat == mat:  # n4

        ua_diag = diag_samemat(p, 12, p, 3, dx, g_val, k_val[p, 12], k_val[p, 3])

    i_cor.append(p * 24 + 12)
    j_cor.append((p + ny) * 24 + 8)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 12)
    j_cor.append(p * 24 + 13)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 12)
    j_cor.append(p * 24 + 15)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 12)
    j_cor.append(p * 24 + 3)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 12)
    j_cor.append(p * 24 + 12)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def e2(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("e2")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + ny, 13].mat != mat:  # w4

        ua_perp = perp_diffmat(p, 13, p + ny, 11, dx, g_val, k_val[p, 13], htc)

    if voxel_db[p + ny, 13].mat == mat:  # w4

        ua_perp = perp_samemat(p, 13, p + ny, 11, dx, g_val, k_val[p, 13], k_val[p + ny, 11])

    if voxel_db[p, 16].mat != mat:  # e3

        ua_s1 = side_diffmat(p, 13, p, 14, dx, g_val, k_val[p, 13], htc)

    if voxel_db[p, 16].mat == mat:  # e3

        ua_s1 = side_samemat(p, 13, p, 14, dx, g_val, k_val[p, 13], k_val[p, 14])

    if voxel_db[p, 14].mat != mat:  # e1

        ua_s2 = side_diffmat(p, 13, p, 12, dx, g_val, k_val[p, 13], htc)

    if voxel_db[p, 14].mat == mat:  # e1

        ua_s2 = side_samemat(p, 13, p, 12, dx, g_val, k_val[p, 13], k_val[p, 12])

    if voxel_db[p, 21].mat != mat:  # f4

        ua_diag = diag_diffmat(p, 13, p, 19, dx, g_val, k_val[p, 13], htc)

    if voxel_db[p, 21].mat == mat:  # f4

        ua_diag = diag_samemat(p, 13, p, 19, dx, g_val, k_val[p, 13], k_val[p, 19])

    i_cor.append(p * 24 + 13)
    j_cor.append((p + ny) * 24 + 11)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 13)
    j_cor.append(p * 24 + 14)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 13)
    j_cor.append(p * 24 + 12)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 13)
    j_cor.append(p * 24 + 19)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 13)
    j_cor.append(p * 24 + 13)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def e3(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("e3")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + ny, 12].mat != mat:  # w3

        ua_perp = perp_diffmat(p, 14, p + ny, 10, dx, g_val, k_val[p, 14], htc)

    if voxel_db[p + ny, 12].mat == mat:  # w3

        ua_perp = perp_samemat(p, 14, p + ny, 10, dx, g_val, k_val[p, 14], k_val[p + ny, 10])

    if voxel_db[p, 17].mat != mat:  # e4

        ua_s1 = side_diffmat(p, 14, p, 15, dx, g_val, k_val[p, 14], htc)

    if voxel_db[p, 17].mat == mat:  # e4

        ua_s1 = side_samemat(p, 14, p, 15, dx, g_val, k_val[p, 14], k_val[p, 15])

    if voxel_db[p, 15].mat != mat:  # e2

        ua_s2 = side_diffmat(p, 14, p, 13, dx, g_val, k_val[p, 14], htc)

    if voxel_db[p, 15].mat == mat:  # e2

        ua_s2 = side_samemat(p, 14, p, 13, dx, g_val, k_val[p, 14], k_val[p, 13])

    if voxel_db[p, 9].mat != mat:  # s4

        ua_diag = diag_diffmat(p, 14, p, 7, dx, g_val, k_val[p, 14], htc)

    if voxel_db[p, 9].mat == mat:  # s4

        ua_diag = diag_samemat(p, 14, p, 7, dx, g_val, k_val[p, 14], k_val[p, 7])

    i_cor.append((p * 24 + 14))
    j_cor.append((p + ny) * 24 + 10)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append((p * 24 + 14))
    j_cor.append(p * 24 + 15)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append((p * 24 + 14))
    j_cor.append(p * 24 + 13)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append((p * 24 + 14))
    j_cor.append(p * 24 + 7)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append((p * 24 + 14))
    j_cor.append((p * 24 + 14))
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def e4(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("e4")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + ny, 11].mat != mat:  # w2

        ua_perp = perp_diffmat(p, 15, p + ny, 9, dx, g_val, k_val[p, 15], htc)

    if voxel_db[p + ny, 11].mat == mat:  # w2

        ua_perp = perp_samemat(p, 15, p + ny, 9, dx, g_val, k_val[p, 15], k_val[p + ny, 9])

    if voxel_db[p, 14].mat != mat:  # e1

        ua_s1 = side_diffmat(p, 15, p, 12, dx, g_val, k_val[p, 15], htc)

    if voxel_db[p, 14].mat == mat:  # e1

        ua_s1 = side_samemat(p, 15, p, 12, dx, g_val, k_val[p, 15], k_val[p, 12])

    if voxel_db[p, 16].mat != mat:  # e3

        ua_s2 = side_diffmat(p, 15, p, 14, dx, g_val, k_val[p, 15], htc)

    if voxel_db[p, 16].mat == mat:  # e3

        ua_s2 = side_samemat(p, 15, p, 14, dx, g_val, k_val[p, 15], k_val[p, 14])

    if voxel_db[p, 23].mat != mat:  # b2

        ua_diag = diag_diffmat(p, 15, p, 21, dx, g_val, k_val[p, 15], htc)

    if voxel_db[p, 23].mat == mat:  # b2

        ua_diag = diag_samemat(p, 15, p, 21, dx, g_val, k_val[p, 15], k_val[p, 21])

    i_cor.append(p * 24 + 15)
    j_cor.append((p + ny) * 24 + 9)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 15)
    j_cor.append(p * 24 + 12)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 15)
    j_cor.append(p * 24 + 14)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 15)
    j_cor.append(p * 24 + 21)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 15)
    j_cor.append(p * 24 + 15)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


# front Tetras ############################################################

def f1(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("f1")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + 1, 22].mat != mat:  # b1

        ua_perp = perp_diffmat(p, 16, p + 1, 20, dx, g_val, k_val[p, 16], htc)

    if voxel_db[p + 1, 22].mat == mat:  # b1

        ua_perp = perp_samemat(p, 16, p + 1, 20, dx, g_val, k_val[p, 16], k_val[p + 1, 20])

    if voxel_db[p, 19].mat != mat:  # f2

        ua_s1 = side_diffmat(p, 16, p, 17, dx, g_val, k_val[p, 16], htc)

    if voxel_db[p, 19].mat == mat:  # f2

        ua_s1 = side_samemat(p, 16, p, 17, dx, g_val, k_val[p, 16], k_val[p, 17])

    if voxel_db[p, 21].mat != mat:  # f4

        ua_s2 = side_diffmat(p, 16, p, 19, dx, g_val, k_val[p, 16], htc)

    if voxel_db[p, 21].mat == mat:  # f4

        ua_s2 = side_samemat(p, 16, p, 19, dx, g_val, k_val[p, 16], k_val[p, 19])

    if voxel_db[p, 4].mat != mat:  # n3

        ua_diag = diag_diffmat(p, 16, p, 2, dx, g_val, k_val[p, 16], htc)

    if voxel_db[p, 4].mat == mat:  # n3

        ua_diag = diag_samemat(p, 16, p, 2, dx, g_val, k_val[p, 16], k_val[p, 2])

    i_cor.append(p * 24 + 16)
    j_cor.append((p + 1) * 24 + 20)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 16)
    j_cor.append(p * 24 + 17)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 16)
    j_cor.append(p * 24 + 19)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 16)
    j_cor.append(p * 24 + 2)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 16)
    j_cor.append(p * 24 + 16)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def f2(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("f2")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + 1, 25].mat != mat:  # b4

        ua_perp = perp_diffmat(p, 17, p + 1, 23, dx, g_val, k_val[p, 17], htc)

    if voxel_db[p + 1, 25].mat == mat:  # b4

        ua_perp = perp_samemat(p, 17, p + 1, 23, dx, g_val, k_val[p, 17], k_val[p + 1, 23])

    if voxel_db[p, 20].mat != mat:  # f3

        ua_s1 = side_diffmat(p, 17, p, 18, dx, g_val, k_val[p, 17], htc)

    if voxel_db[p, 20].mat == mat:  # f3

        ua_s1 = side_samemat(p, 17, p, 18, dx, g_val, k_val[p, 17], k_val[p, 18])

    if voxel_db[p, 18].mat != mat:  # f1

        ua_s2 = side_diffmat(p, 17, p, 16, dx, g_val, k_val[p, 17], htc)

    if voxel_db[p, 18].mat == mat:  # f1

        ua_s2 = side_samemat(p, 17, p, 16, dx, g_val, k_val[p, 17], k_val[p, 16])

    if voxel_db[p, 13].mat != mat:  # w4

        ua_diag = diag_diffmat(p, 17, p, 11, dx, g_val, k_val[p, 17], htc)

    if voxel_db[p, 13].mat == mat:  # w4

        ua_diag = diag_samemat(p, 17, p, 11, dx, g_val, k_val[p, 17], k_val[p, 11])

    i_cor.append(p * 24 + 17)
    j_cor.append((p + 1) * 24 + 23)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 17)
    j_cor.append(p * 24 + 18)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 17)
    j_cor.append(p * 24 + 16)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 17)
    j_cor.append(p * 24 + 11)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 17)
    j_cor.append(p * 24 + 17)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def f3(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("f3")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + 1, 24].mat != mat:  # b3

        ua_perp = perp_diffmat(p, 18, p + 1, 22, dx, g_val, k_val[p, 18], htc)

    if voxel_db[p + 1, 24].mat == mat:  # b3

        ua_perp = perp_samemat(p, 18, p + 1, 22, dx, g_val, k_val[p, 18], k_val[p + 1, 22])

    if voxel_db[p, 21].mat != mat:  # f4

        ua_s1 = side_diffmat(p, 18, p, 19, dx, g_val, k_val[p, 18], htc)

    if voxel_db[p, 21].mat == mat:  # f4

        ua_s1 = side_samemat(p, 18, p, 19, dx, g_val, k_val[p, 18], k_val[p, 19])

    if voxel_db[p, 19].mat != mat:  # f2

        ua_s2 = side_diffmat(p, 18, p, 17, dx, g_val, k_val[p, 18], htc)

    if voxel_db[p, 19].mat == mat:  # f2

        ua_s2 = side_samemat(p, 18, p, 17, dx, g_val, k_val[p, 18], k_val[p, 17])

    if voxel_db[p, 8].mat != mat:  # s3

        ua_diag = diag_diffmat(p, 18, p, 6, dx, g_val, k_val[p, 18], htc)

    if voxel_db[p, 8].mat == mat:  # s3

        ua_diag = diag_samemat(p, 18, p, 6, dx, g_val, k_val[p, 18], k_val[p, 6])

    i_cor.append(p * 24 + 18)
    j_cor.append((p + 1) * 24 + 22)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 18)
    j_cor.append(p * 24 + 19)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 18)
    j_cor.append(p * 24 + 17)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 18)
    j_cor.append(p * 24 + 6)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 18)
    j_cor.append(p * 24 + 18)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def f4(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("f4")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p + 1, 23].mat != mat:  # b2

        ua_perp = perp_diffmat(p, 19, p + 1, 21, dx, g_val, k_val[p, 19], htc)

    if voxel_db[p + 1, 23].mat == mat:  # b2

        ua_perp = perp_samemat(p, 19, p + 1, 21, dx, g_val, k_val[p, 19], k_val[p + 1, 21])

    if voxel_db[p, 18].mat != mat:  # f1

        ua_s1 = side_diffmat(p, 19, p, 16, dx, g_val, k_val[p, 19], htc)

    if voxel_db[p, 18].mat == mat:  # f1

        ua_s1 = side_samemat(p, 19, p, 16, dx, g_val, k_val[p, 19], k_val[p, 16])

    if voxel_db[p, 20].mat != mat:  # f3

        ua_s2 = side_diffmat(p, 19, p, 18, dx, g_val, k_val[p, 19], htc)

    if voxel_db[p, 20].mat == mat:  # f3

        ua_s2 = side_samemat(p, 19, p, 18, dx, g_val, k_val[p, 19], k_val[p, 18])

    if voxel_db[p, 15].mat != mat:
        ua_diag = diag_diffmat(p, 19, p, 13, dx, g_val, k_val[p, 19], htc)

    if voxel_db[p, 15].mat == mat:
        ua_diag = diag_samemat(p, 19, p, 13, dx, g_val, k_val[p, 19], k_val[p, 13])

    i_cor.append(p * 24 + 19)
    j_cor.append((p + 1) * 24 + 21)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 19)
    j_cor.append(p * 24 + 16)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 19)
    j_cor.append(p * 24 + 18)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 19)
    j_cor.append(p * 24 + 13)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 19)
    j_cor.append(p * 24 + 19)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


# back Tetras #########################################################

def b1(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #   print("b1")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - 1, 18].mat != mat:  # f1

        ua_perp = perp_diffmat(p, 20, p - 1, 16, dx, g_val, k_val[p, 20], htc)

    if voxel_db[p - 1, 18].mat == mat:  # f1

        ua_perp = perp_samemat(p, 20, p - 1, 16, dx, g_val, k_val[p, 20], k_val[p - 1, 16])

    if voxel_db[p, 23].mat != mat:  # b2

        ua_s1 = side_diffmat(p, 20, p, 21, dx, g_val, k_val[p, 20], htc)

    if voxel_db[p, 23].mat == mat:  # b2

        ua_s1 = side_samemat(p, 20, p, 21, dx, g_val, k_val[p, 20], k_val[p, 21])

    if voxel_db[p, 25].mat != mat:  # b4

        ua_s2 = side_diffmat(p, 20, p, 23, dx, g_val, k_val[p, 20], htc)

    if voxel_db[p, 25].mat == mat:  # b4

        ua_s2 = side_samemat(p, 20, p, 23, dx, g_val, k_val[p, 20], k_val[p, 23])

    if voxel_db[p, 2].mat != mat:  # n1

        ua_diag = diag_diffmat(p, 20, p, 0, dx, g_val, k_val[p, 20], htc)

    if voxel_db[p, 2].mat == mat:  # n1

        ua_diag = diag_samemat(p, 20, p, 0, dx, g_val, k_val[p, 20], k_val[p, 0])

    i_cor.append(p * 24 + 20)
    j_cor.append((p - 1) * 24 + 16)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 20)
    j_cor.append(p * 24 + 21)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 20)
    j_cor.append(p * 24 + 23)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 20)
    j_cor.append(p * 24 + 0)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 20)
    j_cor.append(p * 24 + 20)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def b2(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("b2")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - 1, 21].mat != mat:  # f4

        ua_perp = perp_diffmat(p, 21, p - 1, 19, dx, g_val, k_val[p, 21], htc)

    if voxel_db[p - 1, 21].mat == mat:  # f4

        ua_perp = perp_samemat(p, 21, p - 1, 19, dx, g_val, k_val[p, 21], k_val[p - 1, 19])

    if voxel_db[p, 24].mat != mat:  # b3

        ua_s1 = side_diffmat(p, 21, p, 22, dx, g_val, k_val[p, 21], htc)

    if voxel_db[p, 24].mat == mat:  # b3

        ua_s1 = side_samemat(p, 21, p, 22, dx, g_val, k_val[p, 21], k_val[p, 22])

    if voxel_db[p, 22].mat != mat:  # b1

        ua_s2 = side_diffmat(p, 21, p, 20, dx, g_val, k_val[p, 21], htc)

    if voxel_db[p, 22].mat == mat:  # b1

        ua_s2 = side_samemat(p, 21, p, 20, dx, g_val, k_val[p, 21], k_val[p, 20])

    if voxel_db[p, 17].mat != mat:  # e4

        ua_diag = diag_diffmat(p, 21, p, 15, dx, g_val, k_val[p, 21], htc)

    if voxel_db[p, 17].mat == mat:  # e4

        ua_diag = diag_samemat(p, 21, p, 15, dx, g_val, k_val[p, 21], k_val[p, 15])

    i_cor.append(p * 24 + 21)
    j_cor.append((p - 1) * 24 + 19)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 21)
    j_cor.append(p * 24 + 22)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 21)
    j_cor.append(p * 24 + 20)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 21)
    j_cor.append(p * 24 + 15)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 21)
    j_cor.append(p * 24 + 21)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def b3(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("b3")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - 1, 20].mat != mat:  # f3

        ua_perp = perp_diffmat(p, 22, p - 1, 18, dx, g_val, k_val[p, 22], htc)

    if voxel_db[p - 1, 20].mat == mat:  # f3

        ua_perp = perp_samemat(p, 22, p - 1, 18, dx, g_val, k_val[p, 22], k_val[p - 1, 18])

    if voxel_db[p, 25].mat != mat:  # b4

        ua_s1 = side_diffmat(p, 22, p, 23, dx, g_val, k_val[p, 22], htc)

    if voxel_db[p, 25].mat == mat:  # b4

        ua_s1 = side_samemat(p, 22, p, 23, dx, g_val, k_val[p, 22], k_val[p, 23])

    if voxel_db[p, 23].mat != mat:  # b2

        ua_s2 = side_diffmat(p, 22, p, 21, dx, g_val, k_val[p, 22], htc)

    if voxel_db[p, 23].mat == mat:  # b2

        ua_s2 = side_samemat(p, 22, p, 21, dx, g_val, k_val[p, 22], k_val[p, 21])

    if voxel_db[p, 6].mat != mat:  # s3

        ua_diag = diag_diffmat(p, 22, p, 4, dx, g_val, k_val[p, 22], htc)

    if voxel_db[p, 6].mat == mat:  # s3

        ua_diag = diag_samemat(p, 22, p, 4, dx, g_val, k_val[p, 22], k_val[p, 4])

    i_cor.append(p * 24 + 22)
    j_cor.append((p - 1) * 24 + 18)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 22)
    j_cor.append(p * 24 + 23)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 22)
    j_cor.append(p * 24 + 21)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 22)
    j_cor.append(p * 24 + 4)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 22)
    j_cor.append(p * 24 + 22)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def b4(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data):
    #    print("b4")
    ua_perp = ""
    ua_s1 = ""
    ua_s2 = ""
    ua_diag = ""
    if voxel_db[p - 1, 19].mat != mat:  # f2

        ua_perp = perp_diffmat(p, 23, p - 1, 17, dx, g_val, k_val[p, 23], htc)

    if voxel_db[p - 1, 19].mat == mat:  # f2

        ua_perp = perp_samemat(p, 23, p - 1, 17, dx, g_val, k_val[p, 23], k_val[p - 1, 17])

    if voxel_db[p, 22].mat != mat:  # b1

        ua_s1 = side_diffmat(p, 23, p, 20, dx, g_val, k_val[p, 23], htc)

    if voxel_db[p, 22].mat == mat:  # b1

        ua_s1 = side_samemat(p, 23, p, 20, dx, g_val, k_val[p, 23], k_val[p, 20])

    if voxel_db[p, 24].mat != mat:  # b3

        ua_s2 = side_diffmat(p, 23, p, 22, dx, g_val, k_val[p, 23], htc)

    if voxel_db[p, 24].mat == mat:  # b3

        ua_s2 = side_samemat(p, 23, p, 22, dx, g_val, k_val[p, 23], k_val[p, 22])

    if voxel_db[p, 11].mat != mat:  # w2

        ua_diag = diag_diffmat(p, 23, p, 9, dx, g_val, k_val[p, 23], htc)

    if voxel_db[p, 11].mat == mat:  # w2

        ua_diag = diag_samemat(p, 23, p, 9, dx, g_val, k_val[p, 23], k_val[p, 9])

    i_cor.append(p * 24 + 23)
    j_cor.append((p - 1) * 24 + 17)
    ua_data.append(ua_perp)  # ua_perp
    i_cor.append(p * 24 + 23)
    j_cor.append(p * 24 + 20)
    ua_data.append(ua_s1)  # ua_s1
    i_cor.append(p * 24 + 23)
    j_cor.append(p * 24 + 22)
    ua_data.append(ua_s2)  # ua_s2
    i_cor.append(p * 24 + 23)
    j_cor.append(p * 24 + 9)
    ua_data.append(ua_diag)  # ua_diag
    i_cor.append(p * 24 + 23)
    j_cor.append(p * 24 + 23)
    ua_data.append(-(ua_perp + ua_s1 + ua_s2 + ua_diag))  # sum of ua
    q_val.append(-q)


def func(tag, voxel_db, p, nx, ny, nz, dx, mat, g_val, htc, k_val, q, q_val, poi, i_cor, j_cor, ua_data):
    if tag == 'n1':
        (n1(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'n2':
        (n2(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'n3':
        (n3(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'n4':
        (n4(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))

    elif tag == 's1':
        (s1(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 's2':
        (s2(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 's3':
        (s3(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 's4':
        (s4(voxel_db, p, nx, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))

    elif tag == 'w1':
        (w1(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'w2':
        (w2(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'w3':
        (w3(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'w4':
        (w4(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))

    elif tag == 'e1':
        (e1(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'e2':
        (e2(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'e3':
        (e3(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'e4':
        (e4(voxel_db, p, ny, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))

    elif tag == 'f1':
        (f1(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'f2':
        (f2(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'f3':
        (f3(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'f4':
        (f4(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))

    elif tag == 'b1':
        (b1(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'b2':
        (b2(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'b3':
        (b3(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))
    elif tag == 'b4':
        (b4(voxel_db, p, dx, mat, g_val, htc, k_val, q, q_val, i_cor, j_cor, ua_data))


# ua matrix generator function

def matrixgenerator(voxel_db, k_val, htc, dx, g_val, nx, ny, nz, q, tamb):
    i_cor = []
    j_cor = []
    ua_data = []
    q_val = []
    vc = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                for t in range(24):
                    poi = 24 * i + 24 * ny * j + 24 * ny * nz * k + t
                    if voxel_db[vc, t + 2].mat == 0:
                        i_cor.append(vc * 24 + t)
                        j_cor.append(vc * 24 + t)
                        ua_data.append(1.0)
                        q_val.append(tamb)

                    if voxel_db[vc, t + 2].mat == 1:
                        loc = voxel_db[vc, (t + 2)].pos
                        func(loc, voxel_db,
                             vc, nx, ny, nz, dx,
                             voxel_db[vc, (t + 2)].mat,
                             g_val, htc, k_val, q, q_val, poi, i_cor,
                             j_cor, ua_data)

                vc = vc + 1

    ua_csc = sp.csc_matrix((ua_data, (i_cor, j_cor)))
    print(len(i_cor), len(j_cor))
    print("memory of i and j", sys.getsizeof(i_cor), sys.getsizeof(j_cor))
    return ua_csc, np.asarray(q_val)
