import math
import numpy
# import dispy
# import cProfile
import visualization
import dictionary_surfacearea as dicsurfarea


# import timeit as timeit
# start = timeit.default_timer()


class Coordinate:
    def __init__(self, x_coord, y_coord, z_coord):
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord


class Tetra:
    def __init__(self, point1, point2, point3, point4, material, position):
        self.p1 = point1
        self.p2 = point2
        self.p3 = point3
        self.p4 = point4
        self.mat = material
        self.pos = position


decimal = 5
a = b = c = 1  # radius
smoothness = dx = dy = dz = 0.15  # smoothness
dr = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
ne = 3
nx = int(2 * a / dx + ne)
ny = int(2 * b / dy + ne)
nz = int(2 * c / dz + ne)

cx = int(nx / 2)
cy = int(ny / 2)
cz = int(nz / 2)
voxel = numpy.zeros((nx, ny, nz), dtype=int)
voxel_n = nx * ny * nz
voxel_db = numpy.empty((voxel_n, 26), dtype=object)
g_coord = -1
# cProfile.run('main()')


def main():
    global voxel
    global voxel_n
    global voxel_db
    global g_coord
    calculate_surface_area = True
    use_visualization = True
    calculate_volume = True
    tetra_mode = True

    voxel_count = numpy.copy(voxel)
    count = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                voxel_count[i, j, k] = count
                count = count + 1
                if (((i - cx) / (a / dx)) ** 2) + (((j - cy) / (b / dy)) ** 2) + (((k - cz) / (c / dz)) ** 2) <= 1:
                    # if(i>(0.5*a/dx) and i<2*a/(dx) and j>(0.5*b/dx) and j<2*b/dx and k>0.5*c/dz and k<2*c/dz):
                    voxel[i, j, k] = 1

    # Smoothen
    # tagging the voxel on sides exposed
    side_exposed = numpy.copy(voxel)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if voxel[i, j, k] != 0:
                    count = 0
                    if voxel[i - 1, j, k] == 0:
                        count = count + 1
                    if voxel[i + 1, j, k] == 0:
                        count = count + 1
                    if voxel[i, j - 1, k] == 0:
                        count = count + 1
                    if voxel[i, j + 1, k] == 0:
                        count = count + 1
                    if voxel[i, j, k - 1] == 0:
                        count = count + 1
                    if voxel[i, j, k + 1] == 0:
                        count = count + 1
                    side_exposed[i, j, k] = count

    # Removing 5 sides exposed cubes
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if side_exposed[i, j, k] == 5:
                    voxel[i, j, k] = 0
                    side_exposed[i, j, k] = 0

    # Recounting
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if voxel[i, j, k] != 0:
                    count = 0
                    if voxel[i - 1, j, k] == 0:
                        count = count + 1
                    if voxel[i + 1, j, k] == 0:
                        count = count + 1
                    if voxel[i, j - 1, k] == 0:
                        count = count + 1
                    if voxel[i, j + 1, k] == 0:
                        count = count + 1
                    if voxel[i, j, k - 1] == 0:
                        count = count + 1
                    if voxel[i, j, k + 1] == 0:
                        count = count + 1
                    side_exposed[i, j, k] = count

    # Creating the coordinate matrix
    # 0  : voxel(x,y,z)
    # 1  : Material tag

    # 2  : N1
    # 3  : N2
    # 4  : N3
    # 5  : N4

    # 6  : S1
    # 7  : S2
    # 8  : S3
    # 9  : S4

    # 10 : W1
    # 11 : W2
    # 12 : W3
    # 13 : W4

    # 14 : E1
    # 15 : E2
    # 16 : E3
    # 17 : E4

    # 18 : F1
    # 19 : F2
    # 20 : F3
    # 21 : F4

    # 22 : B1
    # 23 : B2
    # 24 : B3
    # 25 : B4

    v_count = 0
    print("nz", nz)
    print("ny", ny)
    print("nx", nz)
    print("v_count", v_count)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                voxel_db[v_count, 0] = Coordinate(i * dx, j * dy, k * dz)
                voxel_db[v_count, 1] = voxel[i, j, k]

                # Breaking down the i,j,k in vertex coordinates
                a_coord = Coordinate(round((i * dx - dx / 2.0), decimal), round((j * dy + dy / 2.0), decimal),
                                     round((k * dz + dz / 2.0), decimal))
                b_coord = Coordinate(round((i * dx - dx / 2.0), decimal), round((j * dy - dy / 2.0), decimal),
                                     round((k * dz + dz / 2.0), decimal))
                c_coord = Coordinate(round((i * dx + dx / 2.0), decimal), round((j * dy - dy / 2.0), decimal),
                                     round((k * dz + dz / 2.0), decimal))
                d_coord = Coordinate(round((i * dx + dx / 2.0), decimal), round((j * dy + dy / 2.0), decimal),
                                     round((k * dz + dz / 2.0), decimal))
                e_coord = Coordinate(round((i * dx - dx / 2.0), decimal), round((j * dy + dy / 2.0), decimal),
                                     round((k * dz - dz / 2.0), decimal))
                f_coord = Coordinate(round((i * dx - dx / 2.0), decimal), round((j * dy - dy / 2.0), decimal),
                                     round((k * dz - dz / 2.0), decimal))
                g_coord = Coordinate(round((i * dx + dx / 2.0), decimal), round((j * dy - dy / 2.0), decimal),
                                     round((k * dz - dz / 2.0), decimal))
                h_coord = Coordinate(round((i * dx + dx / 2.0), decimal), round((j * dy + dy / 2.0), decimal),
                                     round((k * dz - dz / 2.0), decimal))

                i_val = Coordinate(round((i * dx), decimal), round(j * dy, decimal), round(k * dz, decimal))
                nc = Coordinate(round(i * dx, decimal), round(j * dy, decimal), round((k * dz + dz / 2.0), decimal))
                sc = Coordinate(round(i * dx, decimal), round(j * dy, decimal), round((k * dz - dz / 2.0), decimal))
                wc = Coordinate(round(i * dx, decimal), round((j * dy - dy / 2.0), decimal), round(k * dz, decimal))
                ec = Coordinate(round(i * dx, decimal), round((j * dy + dy / 2.0), decimal), round(k * dz, decimal))
                fc = Coordinate(round((i * dx + dx / 2.0), decimal), round(j * dy, decimal), round(k * dz, decimal))
                bc = Coordinate(round((i * dx - dx / 2.0), decimal), round(j * dy, decimal), round(k * dz, decimal))
                mat_tag = voxel[i, j, k]

                # North Tetras

                voxel_db[v_count, 2] = Tetra(nc, a_coord, b_coord, i_val, mat_tag, 'N1')  # N1 coordinates 2
                voxel_db[v_count, 3] = Tetra(nc, b_coord, c_coord, i_val, mat_tag, 'N2')  # N2 coordinates 3
                voxel_db[v_count, 4] = Tetra(nc, c_coord, d_coord, i_val, mat_tag, 'N3')  # N3 coordinates 4
                voxel_db[v_count, 5] = Tetra(nc, d_coord, a_coord, i_val, mat_tag, 'N4')  # N4 coordinates 5

                # South Tetras
                voxel_db[v_count, 6] = Tetra(sc, e_coord, f_coord, i_val, mat_tag, 'S1')  # S1 coordinates 6
                voxel_db[v_count, 7] = Tetra(sc, f_coord, g_coord, i_val, mat_tag, 'S2')  # S2 coordinates 7
                voxel_db[v_count, 8] = Tetra(sc, g_coord, h_coord, i_val, mat_tag, 'S3')  # S3 coordinates 8
                voxel_db[v_count, 9] = Tetra(sc, h_coord, e_coord, i_val, mat_tag, 'S4')  # S4 coordinates 9

                # West Tetras
                voxel_db[v_count, 10] = Tetra(wc, c_coord, b_coord, i_val, mat_tag, 'W1')  # W1 coordinates 10
                voxel_db[v_count, 11] = Tetra(wc, b_coord, f_coord, i_val, mat_tag, 'W2')  # W2 coordinates 11
                voxel_db[v_count, 12] = Tetra(wc, f_coord, g_coord, i_val, mat_tag, 'W3')  # W3 coordinates 12
                voxel_db[v_count, 13] = Tetra(wc, g_coord, c_coord, i_val, mat_tag, 'W4')  # W4 coordinates 13

                # East Tetras
                voxel_db[v_count, 14] = Tetra(ec, a_coord, d_coord, i_val, mat_tag, 'E1')  # E1 coordinates 14
                voxel_db[v_count, 15] = Tetra(ec, d_coord, h_coord, i_val, mat_tag, 'E2')  # E2 coordinates 15
                voxel_db[v_count, 16] = Tetra(ec, h_coord, e_coord, i_val, mat_tag, 'E3')  # E3 coordinates 16
                voxel_db[v_count, 17] = Tetra(ec, e_coord, a_coord, i_val, mat_tag, 'E4')  # E4 coordinates 17

                # Front Tetras
                voxel_db[v_count, 18] = Tetra(fc, d_coord, c_coord, i_val, mat_tag, 'F1')  # F1 coordinates 18
                voxel_db[v_count, 19] = Tetra(fc, c_coord, g_coord, i_val, mat_tag, 'F2')  # F2 coordinates 19
                voxel_db[v_count, 20] = Tetra(fc, g_coord, h_coord, i_val, mat_tag, 'F3')  # F3 coordinates 20
                voxel_db[v_count, 21] = Tetra(fc, h_coord, d_coord, i_val, mat_tag, 'F4')  # F4 coordinates 21

                # Back Tetras
                voxel_db[v_count, 22] = Tetra(bc, b_coord, a_coord, i_val, mat_tag, 'B1')  # B1 coordinates 22
                voxel_db[v_count, 23] = Tetra(bc, a_coord, e_coord, i_val, mat_tag, 'B2')  # B2 coordinates 23
                voxel_db[v_count, 24] = Tetra(bc, e_coord, f_coord, i_val, mat_tag, 'B3')  # B3 coordinates 24
                voxel_db[v_count, 25] = Tetra(bc, f_coord, b_coord, i_val, mat_tag, 'B4')  # B4 coordinates 25

                v_count = v_count + 1

    # Add centroid

    g_coord = numpy.zeros((voxel_n, 24), dtype=object)
    for vc in range(voxel_n):
        for i in range(24):
            gx = (voxel_db[vc, i + 2].p1.x + voxel_db[vc, i + 2].p2.x + voxel_db[vc, i + 2].p3.x + voxel_db[
                vc, i + 2].p4.x) / 4
            gy = (voxel_db[vc, i + 2].p1.y + voxel_db[vc, i + 2].p2.y + voxel_db[vc, i + 2].p3.y + voxel_db[
                vc, i + 2].p4.y) / 4
            gz = (voxel_db[vc, i + 2].p1.z + voxel_db[vc, i + 2].p2.z + voxel_db[vc, i + 2].p3.z + voxel_db[
                vc, i + 2].p4.z) / 4
            g_coord[vc, i] = Coordinate(gx, gy, gz)

    # Triangle smoothening
    if tetra_mode is True:
        v_count = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    # Three Sides Exposed
                    if side_exposed[i, j, k] > 1:

                        # North Side
                        if voxel[i, j, k + 1] == 0:
                            voxel_db[v_count, 2].mat = 0
                            voxel_db[v_count, 3].mat = 0
                            voxel_db[v_count, 4].mat = 0
                            voxel_db[v_count, 5].mat = 0

                            if side_exposed[i, j + 1, k] > 1:
                                voxel_db[v_count, 14].mat = 0
                            if side_exposed[i, j - 1, k] > 1:
                                voxel_db[v_count, 10].mat = 0
                            if side_exposed[i + 1, j, k] > 1:
                                voxel_db[v_count, 18].mat = 0
                            if side_exposed[i - 1, j, k] > 1:
                                voxel_db[v_count, 22].mat = 0

                        # South Side
                        if voxel[i, j, k - 1] == 0:
                            voxel_db[v_count, 6].mat = 0
                            voxel_db[v_count, 7].mat = 0
                            voxel_db[v_count, 8].mat = 0
                            voxel_db[v_count, 9].mat = 0

                            if side_exposed[i, j + 1, k] > 1:
                                voxel_db[v_count, 16].mat = 0
                            if side_exposed[i, j - 1, k] > 1:
                                voxel_db[v_count, 12].mat = 0
                            if side_exposed[i + 1, j, k] > 1:
                                voxel_db[v_count, 20].mat = 0
                            if side_exposed[i - 1, j, k] > 1:
                                voxel_db[v_count, 24].mat = 0

                        # West Side
                        if voxel[i, j - 1, k] == 0:
                            voxel_db[v_count, 10].mat = 0
                            voxel_db[v_count, 11].mat = 0
                            voxel_db[v_count, 12].mat = 0
                            voxel_db[v_count, 13].mat = 0

                            if side_exposed[i, j, k + 1] > 1:
                                voxel_db[v_count, 3].mat = 0
                            if side_exposed[i, j, k - 1] > 1:
                                voxel_db[v_count, 7].mat = 0
                            if side_exposed[i + 1, j, k] > 1:
                                voxel_db[v_count, 19].mat = 0
                            if side_exposed[i - 1, j, k] > 1:
                                voxel_db[v_count, 25].mat = 0

                        # East Side
                        if voxel[i, j + 1, k] == 0:
                            voxel_db[v_count, 14].mat = 0
                            voxel_db[v_count, 15].mat = 0
                            voxel_db[v_count, 16].mat = 0
                            voxel_db[v_count, 17].mat = 0

                            if side_exposed[i, j, k + 1] > 1:
                                voxel_db[v_count, 5].mat = 0
                            if side_exposed[i, j, k - 1] > 1:
                                voxel_db[v_count, 9].mat = 0
                            if side_exposed[i + 1, j, k] > 1:
                                voxel_db[v_count, 21].mat = 0
                            if side_exposed[i - 1, j, k] > 1:
                                voxel_db[v_count, 23].mat = 0

                        # Front Side
                        if voxel[i + 1, j, k] == 0:
                            voxel_db[v_count, 18].mat = 0
                            voxel_db[v_count, 19].mat = 0
                            voxel_db[v_count, 20].mat = 0
                            voxel_db[v_count, 21].mat = 0

                            if side_exposed[i, j + 1, k] > 1:
                                voxel_db[v_count, 15].mat = 0
                            if side_exposed[i, j - 1, k] > 1:
                                voxel_db[v_count, 13].mat = 0
                            if side_exposed[i, j, k + 1] > 1:
                                voxel_db[v_count, 4].mat = 0
                            if side_exposed[i, j, k - 1] > 1:
                                voxel_db[v_count, 8].mat = 0

                        # Back Side
                        if voxel[i - 1, j, k] == 0:
                            voxel_db[v_count, 22].mat = 0
                            voxel_db[v_count, 23].mat = 0
                            voxel_db[v_count, 24].mat = 0
                            voxel_db[v_count, 25].mat = 0

                            if side_exposed[i, j + 1, k] > 1:
                                voxel_db[v_count, 17].mat = 0
                            if side_exposed[i, j - 1, k] > 1:
                                voxel_db[v_count, 11].mat = 0
                            if side_exposed[i, j, k + 1] > 1:
                                voxel_db[v_count, 2].mat = 0
                            if side_exposed[i, j, k - 1] > 1:
                                voxel_db[v_count, 6].mat = 0

                    v_count = v_count + 1

    # Visualization
    if use_visualization is True:
        visualization.visualize(voxel_db, side_exposed, nx, ny, nz)

    # Calculate the surface area
    if calculate_surface_area is True:
        d = 0
        first_area_sum = 0.0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if side_exposed[i, j, k] != 0:
                        for m in range(2, 26):
                            if voxel_db[d, m].mat != 0:
                                first_area_sum = first_area_sum + dicsurfarea.func(voxel_db[d, m].pos, voxel_db, d, nx,
                                                                                   ny, nz, dx,
                                                                                   voxel_db[d, m].mat)
                    d = d + 1
        print("Tetra Area = ", first_area_sum)

        voxel_area = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if voxel[i, j, k] == 1:
                        voxel_area = side_exposed[i, j, k] + voxel_area
        main_area_sum = 4 * math.pi * (((a * b) ** 1.6 + (a * c) ** 1.6 + (b * c) ** 1.6) / 3.0) ** (1 / 1.6)
        print("Actual Area = ", main_area_sum)
        voxel_area = voxel_area * dx * dx
        print("voxel area = ", voxel_area)

    # Calculate volume
    if calculate_volume is True:

        volume = 0
        print("voxel_n is", voxel_n)
        # print("vc is", vc)
        for vc in range(voxel_n):
            for i in range(2, 26):
                # print("voxel value", voxel_db[vc, i])
                if voxel_db[vc, i].mat == 1:
                    volume = volume + 1
        if volume == 0:
            raise ValueError('Error, 0 volume!')
        print("Vals")
        print(volume, dx, dy, dz)
        tetra_volume = volume * dx * dy * dz / 24
        print("\ntetra volume = ", tetra_volume)

        voxel_volume = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if voxel[i, j, k] == 1:
                        voxel_volume = voxel_volume + dx * dy * dz

        print("voxel volume = ", voxel_volume)

        main_volume = 4 / 3 * math.pi * a * b * c
        print("actual volume = ", main_volume, "\n")


# stop = timeit.default_timer()
# print("dx = ",dx,"time = ",stop-start,"seconds")
# print("domain size =",voxel_n*24)
# print("matrix size =",voxel_n*24*voxel_n*24)

# Steady State Heat Transfer Solver

def voxel_database():
    return voxel_db, voxel, dx, dy, dz, voxel_n, g_coord


if __name__ == "__main__":
    main()
