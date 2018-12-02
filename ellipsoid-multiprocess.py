import math as mt
import numpy as np
import dictionary_surfacearea as dicsurfarea
import visualization
import timeit as timeit
import multiprocessing as mproc

# initializes timer
start = timeit.default_timer()

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

# Smoothening
# tagging the voxel on sides exposed
def smoothening():
    side_exposed = np.copy(voxel)

    for i in range(nx):
        x1 = i-1
        y1 = i+1
        for j in range(ny):
            x2 = j-1
            y2 = j+1
            for k in range(nz):
                x3 = k-1
                y3 = k+1
                if voxel[i, j, k] != 0:
                    count = 0
                    if voxel[x1, j, k] == 0:
                        count = count + 1
                    if voxel[y1, j, k] == 0:
                        count = count + 1
                    if voxel[i, x2, k] == 0:
                        count = count + 1
                    if voxel[i, y2, k] == 0:
                        count = count + 1
                    if voxel[i, j, x3] == 0:
                        count = count + 1
                    if voxel[i, j, y3] == 0:
                        count = count + 1
                    side_exposed[i, j, k] = count
                    return side_exposed

# Removing 5 sides exposed cubes
def remove_exposed(side_exposed):
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if side_exposed[i, j, k] == 5:
                    voxel[i, j, k] = 0
                    side_exposed[i, j, k] = 0

# Recounting
def recounting():
    for i in range(nx):
        a1 = i-1
        b1 = i+1
        for j in range(ny):
            a2 = j-1
            b2 = j+1
            for k in range(nz):
                a3 = k-1
                b3 = k+1
                if voxel[i, j, k] != 0:
                    count = 0
                    if voxel[a1, j, k] == 0:
                        count = count + 1
                    if voxel[b1, j, k] == 0:
                        count = count + 1
                    if voxel[i, a2, k] == 0:
                        count = count + 1
                    if voxel[i, b2, k] == 0:
                        count = count + 1
                    if voxel[i, j, a3] == 0:
                        count = count + 1
                    if voxel[i, j, b3] == 0:
                        count = count + 1
                    side_exposed[i, j, k] = count

    return nx * ny * nz

# Creating the coordinate matrix
def create_matrix(voxel_n):
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
    global voxel_db
    voxel_db = np.empty((voxel_n, 26), dtype=object)

    v_count = 0
    divx = dx / 2.0
    divy = dy / 2.0
    divz = dz / 2.0

    for i in range(nx):
        multi = i * dx
        for j in range(ny):
            multj = j * dy
            for k in range(nz):
                multk = k * dz
                voxel_db[v_count, 0] = Coordinate(i * dx, j * dy, k * dz)
                voxel_db[v_count, 1] = voxel[i, j, k]

                # Breaking down the i,j,k in vertex coordinates
                A = Coordinate(round((multi - divx), decimal), round((multj + divy), decimal), round((multk + divz), decimal))
                B = Coordinate(round((multi - divx), decimal), round((multj - divy), decimal), round((multk + divz), decimal))
                C = Coordinate(round((multi + divx), decimal), round((multj - divy), decimal), round((multk + divz), decimal))
                D = Coordinate(round((multi + divx), decimal), round((multj + divy), decimal), round((multk + divz), decimal))
                E = Coordinate(round((multi - divx), decimal), round((multj + divy), decimal), round((multk - divz), decimal))
                F = Coordinate(round((multi - divx), decimal), round((multj - divy), decimal), round((multk - divz), decimal))
                G = Coordinate(round((multi + divx), decimal), round((multj - divy), decimal), round((multk - divz), decimal))
                H = Coordinate(round((multi + divx), decimal), round((multj + divy), decimal), round((multk - divz), decimal))
                
                I_val = Coordinate(round(multi, decimal), round(multj, decimal), round(multk, decimal))
                Nc = Coordinate(round(multi, decimal), round(multj, decimal), round((multk + divz), decimal))
                Sc = Coordinate(round(multi, decimal), round(multj, decimal), round((multk - divz), decimal))
                Wc = Coordinate(round(multi, decimal), round((multj - divy), decimal), round(multk, decimal))
                Ec = Coordinate(round(multi, decimal), round((multj + divy), decimal), round(multk, decimal))
                Fc = Coordinate(round((multi + divx), decimal), round(multj, decimal), round(multk, decimal))
                Bc = Coordinate(round((multi - divx), decimal), round(multj, decimal), round(multk, decimal))
                mat_tag = voxel[i, j, k]

                # North Tetras

                voxel_db[v_count, 2] = Tetra(Nc, A, B, I_val, mat_tag, 'N1')  # N1 coordinates 2
                voxel_db[v_count, 3] = Tetra(Nc, B, C, I_val, mat_tag, 'N2')  # N2 coordinates 3
                voxel_db[v_count, 4] = Tetra(Nc, C, D, I_val, mat_tag, 'N3')  # N3 coordinates 4
                voxel_db[v_count, 5] = Tetra(Nc, D, A, I_val, mat_tag, 'N4')  # N4 coordinates 5

                # South Tetras
                voxel_db[v_count, 6] = Tetra(Sc, E, F, I_val, mat_tag, 'S1')  # S1 coordinates 6
                voxel_db[v_count, 7] = Tetra(Sc, F, G, I_val, mat_tag, 'S2')  # S2 coordinates 7
                voxel_db[v_count, 8] = Tetra(Sc, G, H, I_val, mat_tag, 'S3')  # S3 coordinates 8
                voxel_db[v_count, 9] = Tetra(Sc, H, E, I_val, mat_tag, 'S4')  # S4 coordinates 9

                # West Tetras
                voxel_db[v_count, 10] = Tetra(Wc, C, B, I_val, mat_tag, 'W1')  # W1 coordinates 10
                voxel_db[v_count, 11] = Tetra(Wc, B, F, I_val, mat_tag, 'W2')  # W2 coordinates 11
                voxel_db[v_count, 12] = Tetra(Wc, F, G, I_val, mat_tag, 'W3')  # W3 coordinates 12
                voxel_db[v_count, 13] = Tetra(Wc, G, C, I_val, mat_tag, 'W4')  # W4 coordinates 13

                # East Tetras
                voxel_db[v_count, 14] = Tetra(Ec, A, D, I_val, mat_tag, 'E1')  # E1 coordinates 14
                voxel_db[v_count, 15] = Tetra(Ec, D, H, I_val, mat_tag, 'E2')  # E2 coordinates 15
                voxel_db[v_count, 16] = Tetra(Ec, H, E, I_val, mat_tag, 'E3')  # E3 coordinates 16
                voxel_db[v_count, 17] = Tetra(Ec, E, A, I_val, mat_tag, 'E4')  # E4 coordinates 17

                # Front Tetras
                voxel_db[v_count, 18] = Tetra(Fc, D, C, I_val, mat_tag, 'F1')  # F1 coordinates 18
                voxel_db[v_count, 19] = Tetra(Fc, C, G, I_val, mat_tag, 'F2')  # F2 coordinates 19
                voxel_db[v_count, 20] = Tetra(Fc, G, H, I_val, mat_tag, 'F3')  # F3 coordinates 20
                voxel_db[v_count, 21] = Tetra(Fc, H, D, I_val, mat_tag, 'F4')  # F4 coordinates 21

                # Back Tetras
                voxel_db[v_count, 22] = Tetra(Bc, B, A, I_val, mat_tag, 'B1')  # B1 coordinates 22
                voxel_db[v_count, 23] = Tetra(Bc, A, E, I_val, mat_tag, 'B2')  # B2 coordinates 23
                voxel_db[v_count, 24] = Tetra(Bc, E, F, I_val, mat_tag, 'B3')  # B3 coordinates 24
                voxel_db[v_count, 25] = Tetra(Bc, F, B, I_val, mat_tag, 'B4')  # B4 coordinates 25

                v_count = v_count + 1

# Add centroid
def add_centroid():
    G = np.zeros((voxel_n, 24), dtype=object)
    for vc in range(voxel_n):
        for i in range(24):
            gx = (voxel_db[vc, i + 2].p1.x + voxel_db[vc, i + 2].p2.x + voxel_db[vc, i + 2].p3.x + voxel_db[
                vc, i + 2].p4.x) / 4
            gy = (voxel_db[vc, i + 2].p1.y + voxel_db[vc, i + 2].p2.y + voxel_db[vc, i + 2].p3.y + voxel_db[
                vc, i + 2].p4.y) / 4
            gz = (voxel_db[vc, i + 2].p1.z + voxel_db[vc, i + 2].p2.z + voxel_db[vc, i + 2].p3.z + voxel_db[
                vc, i + 2].p4.z) / 4
            G[vc, i] = Coordinate(gx, gy, gz)

# Triangle smoothening
def triangle_smoothening():
    if tetramode == 1:
        v_count = 0
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
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

# Calculate the surface area
def calc_surface_area():
    if surfaceareacalc == 1:
        d = 0
        areasum1 = 0.0
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if side_exposed[i, j, k] != 0:
                        for m in range(2, 26):
                            if voxel_db[d, m].mat != 0:
                                areasum1 = areasum1 + dicsurfarea.func(voxel_db[d, m].pos, voxel_db, d, nx, ny, nz, dx,
                                                                    voxel_db[d, m].mat)
                    d = d + 1
        

        voxelarea = float(0)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if voxel[i, j, k] == 1:
                        voxelarea = side_exposed[i, j, k] + voxelarea
        actualarea = float(4 * mt.pi * (((a * b) ** 1.6 + (a * c) ** 1.6 + (b * c) ** 1.6) / 3.0) ** (1 / 1.6))
        
        voxelarea = voxelarea * dx * dx
        print("Tetra Area = ", areasum1)
        print("Actual Area = ", actualarea)
        print("voxel area = ", voxelarea)

# Calculate volume
def calc_volume():
    if volumecalc == 1:

        volume = 0
        # print("vc is", vc)
        for vc in range(voxel_n):
            for i in range(2, 26):
                # print("voxel value", voxel_db[vc, i])
                if voxel_db[vc, i].mat == 1:
                    volume = volume + 1

        if volume == 0:
            raise ValueError('Error, 0 volume!')

        tetra_volume = volume * dx * dy * dz / 24        

        voxel_volume = float(0)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if voxel[i, j, k] == 1:
                        voxel_volume = voxel_volume + dx * dy * dz

        actualvolume = 4 / 3 * mt.pi * a * b * c
        print("\nTetra volume = ", tetra_volume)
        print("Voxel volume = ", voxel_volume)
        print("Actual volume = ", actualvolume, "\n")

# stops timer and displays the calculated values
def print_val():
    stop = timeit.default_timer()
    print("dx = ",dx,"time = ",stop-start,"seconds")
    print("domain size =",voxel_n*24)
    print("matrix size =",voxel_n*24*voxel_n*24)

# toggles visualization 
def display_visual(use_visualization):
    if(use_visualization == 1):
        visualization.visualize(voxel_db, side_exposed, nx, ny, nz)

# Steady State Heat Transfer Solver
def voxeldatabase():
    return voxel_db, voxel, dx, dy, dz, voxel_n, G

#---------------------------------------------
# declares variables and starts the program
#----------------------------------------------
global side_exposed

decimal = 5
a = 2.1
b = 2.1
c = 2.1
# R = 0.1
dx = dy = dz = 0.06
dr = mt.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

# 1 for Yes 0 for No
tetramode = 1
surfaceareacalc = 1
volumecalc = 1
use_visualization = 1

ne = 3
nx = int(2 * a / dx + ne)
ny = int(2 * b / dy + ne)
nz = int(2 * c / dz + ne)

cx = int(nx / 2)
cy = int(ny / 2)
cz = int(nz / 2)

voxel = np.zeros((nx, ny, nz), dtype=int)
voxel_count = np.copy(voxel)
count = 0

# builds sphere
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            voxel_count[i, j, k] = count
            count = count + 1
            if (((i - cx) / (a / dx)) ** 2) + (((j - cy) / (b / dy)) ** 2) + (((k - cz) / (c / dz)) ** 2) <= 1:
                # if(i>(0.5*a/dx) and i<2*a/(dx) and j>(0.5*b/dx) and j<2*b/dx and k>0.5*c/dz and k<2*c/dz):
                voxel[i, j, k] = 1

# runs the methods declared above
side_exposed = smoothening() #phase 1 
remove_exposed(side_exposed) #phase 2
voxel_n = recounting() #phase 3
create_matrix(voxel_n) #phase 4
triangle_smoothening() #phase 5

# multiprocessing
if __name__ == "__main__":
    phase6 = mproc.Process(target=add_centroid)
    phase7 = mproc.Process(target=calc_surface_area)
    phase8 = mproc.Process(target=calc_volume)

    phase6.start()
    phase7.start()
    phase8.start()
  
    phase6.join()
    phase7.join()
    phase8.join()

print_val()
display_visual(use_visualization)