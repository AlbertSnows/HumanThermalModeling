import numpy as np
import math as mt
#import visualization as v
import dictionary_surfacearea as dicsurfarea
import timeit as timeit

start = timeit.default_timer()

class coordinate:
    def __init__(self,x_coord,y_coord,z_coord):
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord

class tetra:
    def __init__(self,point1,point2,point3,point4,material,position):
        self.p1 = point1
        self.p2 = point2
        self.p3 = point3
        self.p4 = point4
        self.mat = material
        self.pos = position

decimal = 5
a = 2.0
b = 2.0
c = 2.0
#R = 0.1
dx = dy = dz = 0.05
dr = mt.sqrt(dx**2 + dy**2 + dz**2)

# 1 for Yes 0 for No
tetramode = 1
surfaceareacalc = 1
volumecalc = 1
visualization = 1

ne = 3
nx = int(2*a/dx + ne)
ny = int(2*b/dy + ne)
nz = int(2*c/dz + ne)

cx = int(nx/2)
cy = int(ny/2)
cz = int(nz/2)

voxel = np.zeros((nx,ny,nz),dtype=int)
voxel_count = np.copy(voxel)
count = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            voxel_count[i,j,k] = count 
            count = count + 1
            if((i-cx)**2/(a/dx)**2 + (j-cy)**2/(b/dy)**2 +(k-cz)**2/(c/dz)**2 <= 1):
#            if(i>(0.5*a/dx) and i<2*a/(dx) and j>(0.5*b/dx) and j<2*b/dx and k>0.5*c/dz and k<2*c/dz):
                voxel[i,j,k] = 1

# Smoothening
# tagging the voxel on sides exposed
side_exposed = np.copy(voxel)

for k in range (nz):
    for j in range(ny):
        for i in range(nx):
            if(voxel[i,j,k] != 0):
                count = 0
                if(voxel[i-1,j,k] == 0):
                    count = count + 1
                if(voxel[i+1,j,k] == 0):
                    count = count + 1
                if(voxel[i,j-1,k] == 0):
                    count = count + 1
                if(voxel[i,j+1,k] == 0):
                    count = count + 1
                if(voxel[i,j,k-1] == 0):
                    count = count + 1
                if(voxel[i,j,k+1] == 0):
                    count = count + 1
                side_exposed[i,j,k] = count

# Removing 5 sides exposed cubes
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            if(side_exposed[i,j,k] == 5):
                voxel[i,j,k] = 0
                side_exposed[i,j,k] = 0

# Recounting
for k in range (nz):
    for j in range(ny):
        for i in range(nx):
            if(voxel[i,j,k] != 0):
                count = 0
                if(voxel[i-1,j,k] == 0):
                    count = count + 1
                if(voxel[i+1,j,k] == 0):
                    count = count + 1
                if(voxel[i,j-1,k] == 0):
                    count = count + 1
                if(voxel[i,j+1,k] == 0):
                    count = count + 1
                if(voxel[i,j,k-1] == 0):
                    count = count + 1
                if(voxel[i,j,k+1] == 0):
                    count = count + 1
                side_exposed[i,j,k] = count
                
voxel_n = nx*ny*nz

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

voxel_db = np.empty((voxel_n,26),dtype=object)

v_count = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            voxel_db[v_count,0] = coordinate(i*dx,j*dy,k*dz)
            voxel_db[v_count,1] = voxel[i,j,k]
            
            # Breaking down the i,j,k in vertex coordinates
            A =  coordinate(round((i*dx - dx/2.0),decimal), round((j*dy + dy/2.0),decimal), round((k*dz + dz/2.0),decimal))
            B =  coordinate(round((i*dx - dx/2.0),decimal), round((j*dy - dy/2.0),decimal), round((k*dz + dz/2.0),decimal))
            C =  coordinate(round((i*dx + dx/2.0),decimal), round((j*dy - dy/2.0),decimal), round((k*dz + dz/2.0),decimal))
            D =  coordinate(round((i*dx + dx/2.0),decimal), round((j*dy + dy/2.0),decimal), round((k*dz + dz/2.0),decimal))
            E =  coordinate(round((i*dx - dx/2.0),decimal), round((j*dy + dy/2.0),decimal), round((k*dz - dz/2.0),decimal))
            F =  coordinate(round((i*dx - dx/2.0),decimal) ,round((j*dy - dy/2.0),decimal), round((k*dz - dz/2.0),decimal))
            G =  coordinate(round((i*dx + dx/2.0),decimal), round((j*dy - dy/2.0),decimal), round((k*dz - dz/2.0),decimal))
            H =  coordinate(round((i*dx + dx/2.0),decimal), round((j*dy + dy/2.0),decimal), round((k*dz - dz/2.0),decimal))
            I =  coordinate(round((i*dx),decimal), round(j*dy,decimal), round(k*dz,decimal))
            Nc = coordinate(round(i*dx,decimal), round(j*dy,decimal), round((k*dz + dz/2.0),decimal))
            Sc = coordinate(round(i*dx,decimal), round(j*dy,decimal), round((k*dz - dz/2.0),decimal))
            Wc = coordinate(round(i*dx,decimal), round((j*dy - dy/2.0),decimal), round(k*dz,decimal))
            Ec = coordinate(round(i*dx,decimal),round((j*dy + dy/2.0),decimal), round(k*dz,decimal))
            Fc = coordinate(round((i*dx + dx/2.0),decimal), round(j*dy,decimal), round(k*dz,decimal))
            Bc = coordinate(round((i*dx - dx/2.0),decimal), round(j*dy,decimal), round(k*dz,decimal))
            mat_tag = voxel[i,j,k]
            
            # North Tetras
            
            voxel_db[v_count,2] = tetra(Nc,A,B,I,mat_tag,'N1') # N1 coordinates 2        
            voxel_db[v_count,3] = tetra(Nc,B,C,I,mat_tag,'N2') # N2 coordinates 3
            voxel_db[v_count,4] = tetra(Nc,C,D,I,mat_tag,'N3') # N3 coordinates 4
            voxel_db[v_count,5] = tetra(Nc,D,A,I,mat_tag,'N4') # N4 coordinates 5
            
            # South Tetras
            voxel_db[v_count,6] = tetra(Sc,E,F,I,mat_tag,'S1') # S1 coordinates 6
            voxel_db[v_count,7] = tetra(Sc,F,G,I,mat_tag,'S2') # S2 coordinates 7
            voxel_db[v_count,8] = tetra(Sc,G,H,I,mat_tag,'S3') # S3 coordinates 8
            voxel_db[v_count,9] = tetra(Sc,H,E,I,mat_tag,'S4') # S4 coordinates 9
            
            # West Tetras
            voxel_db[v_count,10] = tetra(Wc,C,B,I,mat_tag,'W1') # W1 coordinates 10
            voxel_db[v_count,11] = tetra(Wc,B,F,I,mat_tag,'W2') # W2 coordinates 11
            voxel_db[v_count,12] = tetra(Wc,F,G,I,mat_tag,'W3') # W3 coordinates 12
            voxel_db[v_count,13] = tetra(Wc,G,C,I,mat_tag,'W4') # W4 coordinates 13
            
            # East Tetras
            voxel_db[v_count,14] = tetra(Ec,A,D,I,mat_tag,'E1') # E1 coordinates 14
            voxel_db[v_count,15] = tetra(Ec,D,H,I,mat_tag,'E2') # E2 coordinates 15
            voxel_db[v_count,16] = tetra(Ec,H,E,I,mat_tag,'E3') # E3 coordinates 16
            voxel_db[v_count,17] = tetra(Ec,E,A,I,mat_tag,'E4') # E4 coordinates 17
            
            # Front Tetras
            voxel_db[v_count,18] = tetra(Fc,D,C,I,mat_tag,'F1') # F1 coordinates 18
            voxel_db[v_count,19] = tetra(Fc,C,G,I,mat_tag,'F2') # F2 coordinates 19
            voxel_db[v_count,20] = tetra(Fc,G,H,I,mat_tag,'F3') # F3 coordinates 20
            voxel_db[v_count,21] = tetra(Fc,H,D,I,mat_tag,'F4') # F4 coordinates 21
            
            # Back Tetras
            voxel_db[v_count,22] = tetra(Bc,B,A,I,mat_tag,'B1') # B1 coordinates 22
            voxel_db[v_count,23] = tetra(Bc,A,E,I,mat_tag,'B2') # B2 coordinates 23
            voxel_db[v_count,24] = tetra(Bc,E,F,I,mat_tag,'B3') # B3 coordinates 24
            voxel_db[v_count,25] = tetra(Bc,F,B,I,mat_tag,'B4') # B4 coordinates 25
            
            v_count = v_count + 1

# Add centroid

G = np.zeros((voxel_n,24),dtype = object)
for vc in range(voxel_n):
    for i in range(24):
        gx = (voxel_db[vc,i+2].p1.x + voxel_db[vc,i+2].p2.x + voxel_db[vc,i+2].p3.x + voxel_db[vc,i+2].p4.x)/4
        gy = (voxel_db[vc,i+2].p1.y + voxel_db[vc,i+2].p2.y + voxel_db[vc,i+2].p3.y + voxel_db[vc,i+2].p4.y)/4
        gz = (voxel_db[vc,i+2].p1.z + voxel_db[vc,i+2].p2.z + voxel_db[vc,i+2].p3.z + voxel_db[vc,i+2].p4.z)/4
        G[vc,i] = coordinate(gx,gy,gz)

        
# Triangle smoothening
if (tetramode == 1):
    v_count = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # Three Sides Exposed
                if (side_exposed[i,j,k] > 1):
                    
                    # North Side
                    if(voxel[i,j,k+1] == 0):
                        voxel_db[v_count,2].mat = 0 
                        voxel_db[v_count,3].mat = 0 
                        voxel_db[v_count,4].mat = 0 
                        voxel_db[v_count,5].mat = 0 
                        
                        if(side_exposed[i,j+1,k] > 1):
                            voxel_db[v_count,14].mat = 0
                        if(side_exposed[i,j-1,k] > 1):
                            voxel_db[v_count,10].mat = 0
                        if(side_exposed[i+1,j,k] > 1):
                            voxel_db[v_count,18].mat = 0
                        if(side_exposed[i-1,j,k] > 1):
                            voxel_db[v_count,22].mat = 0
                    
                    # South Side
                    if(voxel[i,j,k-1] == 0):
                        voxel_db[v_count,6].mat = 0 
                        voxel_db[v_count,7].mat = 0 
                        voxel_db[v_count,8].mat = 0 
                        voxel_db[v_count,9].mat = 0 
                        
                        if(side_exposed[i,j+1,k] > 1):
                            voxel_db[v_count,16].mat = 0
                        if(side_exposed[i,j-1,k] > 1):
                            voxel_db[v_count,12].mat = 0
                        if(side_exposed[i+1,j,k] > 1):
                            voxel_db[v_count,20].mat = 0
                        if(side_exposed[i-1,j,k] > 1):
                            voxel_db[v_count,24].mat = 0
                    
                    # West Side
                    if(voxel[i,j-1,k] == 0):
                        voxel_db[v_count,10].mat = 0 
                        voxel_db[v_count,11].mat = 0 
                        voxel_db[v_count,12].mat = 0 
                        voxel_db[v_count,13].mat = 0 
                        
                        if(side_exposed[i,j,k+1] > 1):
                            voxel_db[v_count,3].mat = 0
                        if(side_exposed[i,j,k-1] > 1):
                            voxel_db[v_count,7].mat = 0
                        if(side_exposed[i+1,j,k] > 1):
                            voxel_db[v_count,19].mat = 0
                        if(side_exposed[i-1,j,k] > 1):
                            voxel_db[v_count,25].mat = 0
                    
                    # East Side
                    if(voxel[i,j+1,k] == 0):
                        voxel_db[v_count,14].mat = 0 
                        voxel_db[v_count,15].mat = 0 
                        voxel_db[v_count,16].mat = 0 
                        voxel_db[v_count,17].mat = 0 
                        
                        if(side_exposed[i,j,k+1] > 1):
                            voxel_db[v_count,5].mat = 0
                        if(side_exposed[i,j,k-1] > 1):
                            voxel_db[v_count,9].mat = 0
                        if(side_exposed[i+1,j,k] > 1):
                            voxel_db[v_count,21].mat = 0
                        if(side_exposed[i-1,j,k] > 1):
                            voxel_db[v_count,23].mat = 0
                    
                    # Front Side
                    if(voxel[i+1,j,k] == 0):
                        voxel_db[v_count,18].mat = 0 
                        voxel_db[v_count,19].mat = 0 
                        voxel_db[v_count,20].mat = 0 
                        voxel_db[v_count,21].mat = 0  
                        
                        if(side_exposed[i,j+1,k] > 1):
                            voxel_db[v_count,15].mat = 0
                        if(side_exposed[i,j-1,k] > 1):
                            voxel_db[v_count,13].mat = 0
                        if(side_exposed[i,j,k+1] > 1):
                            voxel_db[v_count,4].mat = 0
                        if(side_exposed[i,j,k-1] > 1):
                            voxel_db[v_count,8].mat = 0
                        
                    
                    # Back Side
                    if(voxel[i-1,j,k] == 0):
                        voxel_db[v_count,22].mat = 0 
                        voxel_db[v_count,23].mat = 0 
                        voxel_db[v_count,24].mat = 0 
                        voxel_db[v_count,25].mat = 0
                        
                        if(side_exposed[i,j+1,k] > 1):
                            voxel_db[v_count,17].mat = 0
                        if(side_exposed[i,j-1,k] > 1):
                            voxel_db[v_count,11].mat = 0
                        if(side_exposed[i,j,k+1] > 1):
                            voxel_db[v_count,2].mat = 0
                        if(side_exposed[i,j,k-1] > 1):
                            voxel_db[v_count,6].mat = 0
                                      
    
                v_count = v_count + 1

# Visualization
#if(visualization == 1):
#   v.visualize(voxel_db, side_exposed, nx, ny, nz)

# Calculate the surface area
if(surfaceareacalc == 1):
             
    d = 0
    areasum1 = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if(side_exposed[i,j,k] != 0):
                    for m in range(2,26):
                        if(voxel_db[d,m].mat != 0):
                            areasum1 = areasum1 + dicsurfarea.func(voxel_db[d,m].pos,voxel_db,d,nx,ny,dx,voxel_db[d,m].mat)
                d = d + 1
    print("Tetra Area = ", areasum1)
    
    voxelarea = float(0)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if(voxel[i,j,k] == 1):
                    voxelarea = side_exposed[i,j,k] + voxelarea
    actualarea = float(4*mt.pi*(((a*b)**1.6+(a*c)**1.6+(b*c)**1.6)/3.0)**(1/1.6))
    print("Actual Area = ", actualarea)
    voxelarea = voxelarea*dx*dx
    print("voxel area = ", voxelarea)
    
# Calculate volume
if(volumecalc == 1):
    
    volume = 0
    for vc in range(voxel_n):
        for i in range(2,26):
            if(voxel_db[vc,i].mat == 1):
                volume = volume + 1
    
    tetra_volume = volume * dx*dy*dz / 24
    print("\ntetra volume = ",tetra_volume)
    
    voxel_volume = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if(voxel[i,j,k] == 1):
                    voxel_volume = voxel_volume + dx*dy*dz
    
    print("voxel volume = ", voxel_volume)
    
    actualvolume = 4/3*mt.pi*a*b*c
    print("actual volume = ",actualvolume,"\n")

stop = timeit.default_timer()
print("dx = ",dx,"time = ",stop-start,"seconds")
print("domain size =",voxel_n*24)
print("matrix size =",voxel_n*24*voxel_n*24)
print("DATA: ", stop - start, "seconds")

# Steady State Heat Transfer Solver

def voxeldatabase():
    return(voxel_db,voxel,dx,dy,dz,voxel_n,G)