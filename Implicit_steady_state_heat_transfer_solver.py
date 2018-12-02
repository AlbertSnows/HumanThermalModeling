import faulthandler
import sys
import numpy
import scipy.sparse.linalg as spla
import Ellipsoid_optimized_full as Ellipse
import ua_matrix_generator as ua_m

# import scipy.sparse as sp
# import xlsxwriter as xl
# from numpy import linalg as la
# matplotlib.use('Agg')
faulthandler.enable()


def heatsolver(heat_rate, th_cond, htc, tamb):
    voxel_db = Ellipse.voxel_db

    some_val_k = numpy.zeros((Ellipse.voxel_n, 24), dtype=float)

    for voxel_n_element in range(Ellipse.voxel_n):
        for number in range(24):
            if voxel_db[voxel_n_element, (number + 2)].mat != 0:
                some_val_k[voxel_n_element, number] = th_cond

    nx = Ellipse.nx
    ny = Ellipse.ny
    nz = Ellipse.nz

    # n_tetra = nx * ny * nz * 24

    tetra_vol = Ellipse.dx * Ellipse.dy * Ellipse.dz / 24

    q = tetra_vol * heat_rate

    # some_val_q = np.zeros((n_tetra, 1), dtype=float)

    # Creating the UA matrix

    ua, some_val_q = \
        ua_m.matrixgenerator(
            voxel_db, some_val_k, htc,
            Ellipse.dx, Ellipse.g_coord,
            nx, ny, nz, q, tamb)

    return ua, some_val_q


Heat_rate = 1000
HTC = 2.0
k_val = 0.3
Tamb = 20.0

UA, Q = heatsolver(Heat_rate, k_val, HTC, Tamb)

print("Matrix generated, now solving it")

print("memory = ", sys.getsizeof(UA))

T = spla.spsolve(UA, Q)

Tdomtetra = numpy.zeros((Ellipse.voxel_n, 24), dtype=float)
count = 0
for vc in range(Ellipse.voxel_n):
    for i in range(24):
        Tdomtetra[vc, i] = T[count]
        count = count + 1
