import numpy as np
import math as mt
#from numpy import linalg as la
import Ellipsoid_optimized_full as ellipse
import scipy.sparse as sp
import sys
import faulthandler; faulthandler.enable()

#######################################

deci = 5

def perp_samemat(vc1,t1,vc2,t2,dx,G,k1,k2):

    L = round(mt.sqrt((G[vc1,t1].x - G[vc2,t2].x)**2 + (G[vc1,t1].y - G[vc2,t2].y)**2 + (G[vc1,t1].z - G[vc2,t2].z)**2)/2,deci)
    UA_perp = 1/(L/k1 + L/k2)*(0.25*dx**2)
    #print("perp distance", L)
#    print("U same perp",UA_perp/(0.25*dx**2))
    return(UA_perp)

def perp_diffmat(vc1,t1,vc2,t2,dx,G,k1,h):
   
    L = round(mt.sqrt((G[vc1,t1].x - G[vc2,t2].x)**2 + (G[vc1,t1].y - G[vc2,t2].y)**2 + (G[vc1,t1].z - G[vc2,t2].z)**2)/2,deci)
    UA_perp = 1/(L/k1 +1/h)*(0.25*dx**2)
    #print("perp distance", L)
#    print("U diff perp",UA_perp/(0.25*dx**2))
    return(UA_perp)

def side_samemat(vc1,t1,vc2,t2,dx,G,k1,k2):    

    L = round(mt.sqrt((G[vc1,t1].x - G[vc2,t2].x)**2 + (G[vc1,t1].y - G[vc2,t2].y)**2 + (G[vc1,t1].z - G[vc2,t2].z)**2)/2,deci)
    UA_side = 1/(L/k1 + L/k2)*(0.17677*dx**2)
    #print("side distance", L)
#    print("U same side",UA_side/(0.17677*dx**2))
    return(UA_side)

def side_diffmat(vc1,t1,vc2,t2,dx,G,k1,h):    

    L = round(mt.sqrt((G[vc1,t1].x - G[vc2,t2].x)**2 + (G[vc1,t1].y - G[vc2,t2].y)**2 + (G[vc1,t1].z - G[vc2,t2].z)**2)/2,deci)
    UA_side = 1/(L/k1 + 1/h)*(0.17677*dx**2)
   # print("side distance", L)
#    print("U diff side",UA_side/(0.17677*dx**2))
    return(UA_side)
        
def diag_samemat(vc1,t1,vc2,t2,dx,G,k1,k2):
    
    L = round(mt.sqrt((G[vc1,t1].x - G[vc2,t2].x)**2 + (G[vc1,t1].y - G[vc2,t2].y)**2 + (G[vc1,t1].z - G[vc2,t2].z)**2)/2,deci)
    UA_diag = 1/(L/k1 + L/k2)*(0.3535*dx**2)
    #print("diag distance", L)
#    print("U same diag",UA_diag/(0.3535*dx**2))
    return(UA_diag)
    
def diag_diffmat(vc1,t1,vc2,t2,dx,G,k1,h):
    
    L = round(mt.sqrt((G[vc1,t1].x - G[vc2,t2].x)**2 + (G[vc1,t1].y - G[vc2,t2].y)**2 + (G[vc1,t1].z - G[vc2,t2].z)**2)/2,deci)
    UA_diag = 1/(L/k1 + 1/h)*(0.3535*dx**2)
    #print("diag distance", L)
#    print("U diff diag",UA_diag/(0.3535*dx**2))
    return(UA_diag)
    
        
 # North Tetras ####################################################

def N1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("N1")

    if(voxel_db[p + nx*ny,6].mat == 0):       #S1    
    
        UA_perp = perp_diffmat(p,0,p+nx*ny,4,dx,G,K[p,0],HTC)
        
    if(voxel_db[p + nx*ny,6].mat == mat):
        
        UA_perp = perp_samemat(p,0,p+nx*ny,4,dx,G,K[p,0],K[p+nx*ny,4])
               
    if(voxel_db[p,3].mat == 0):               #N2
        
        UA_s1 = side_diffmat(p,0,p,1,dx,G,K[p,0],HTC)
        
    if(voxel_db[p,3].mat == mat):
        
        UA_s1 = side_samemat(p,0,p,1,dx,G,K[p,0],K[p,1])
        
    if(voxel_db[p,5].mat == 0):               #N4        
        
        UA_s2 = side_diffmat(p,0,p,3,dx,G,K[p,0],HTC)
        
    if(voxel_db[p,5].mat == mat):
        
        UA_s2 = side_samemat(p,0,p,3,dx,G,K[p,0],K[p,3])
        
    if(voxel_db[p,22].mat == 0):              #B1
        
        UA_diag = diag_diffmat(p,0,p,20,dx,G,K[p,0],HTC)
    
    if(voxel_db[p,22].mat == mat):
        
        UA_diag = diag_samemat(p,0,p,20,dx,G,K[p,0],K[p,20])
    
    i_cor.append(p*24 + 0) ; j_cor.append(((p+nx*ny)*24+4)) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24 + 0) ; j_cor.append(p*24 + 1) ; ua_data.append(UA_s1) # UA_s1
    i_cor.append(p*24 + 0) ; j_cor.append(p*24+3) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24 + 0) ; j_cor.append(p*24+20) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24 + 0) ; j_cor.append(p*24+0) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
       

def N2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("N2")

    if(voxel_db[p + nx*ny,7].mat == 0):       #S2
        
        UA_perp = perp_diffmat(p,1,p+nx*ny,5,dx,G,K[p,1],HTC)
        
    if(voxel_db[p + nx*ny,7].mat == mat): 
        
        UA_perp = perp_samemat(p,1,p+nx*ny,5,dx,G,K[p,1],K[p+nx*ny,5])
        
    if(voxel_db[p,4].mat == 0):               #N3
        
        UA_s1 = side_diffmat(p,1,p,2,dx,G,K[p,1],HTC)
    
    if(voxel_db[p,4].mat == mat):
        
        UA_s1 = side_samemat(p,1,p,2,dx,G,K[p,1],K[p,2])
        
    if(voxel_db[p,2].mat == 0):               #N1
        
        UA_s2 = side_diffmat(p,1,p,0,dx,G,K[p,1],HTC)
        
    if(voxel_db[p,2].mat == mat): 
        
        UA_s2 = side_samemat(p,1,p,0,dx,G,K[p,1],K[p,0])
        
    if(voxel_db[p,10].mat == 0):              #W1
        
        UA_diag = diag_diffmat(p,1,p,8,dx,G,K[p,1],HTC)
    
    if(voxel_db[p,10].mat == mat):
        
        UA_diag = diag_samemat(p,1,p,8,dx,G,K[p,1],K[p,8])
    
    i_cor.append(p*24+1) ; j_cor.append((p+nx*ny)*24+5) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+1) ; j_cor.append(p*24+2) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+1) ; j_cor.append(p*24+0) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+1) ; j_cor.append(p*24+8) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+1) ; j_cor.append(p*24+1) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    

def N3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("N3")

    if(voxel_db[p + nx*ny,8].mat == 0):       #S3
        
        UA_perp = perp_diffmat(p,2,(p+nx*ny),6,dx,G,K[p,2],HTC)
        
    if(voxel_db[p + nx*ny,8].mat == mat):
        
        UA_perp = perp_samemat(p,2,(p+nx*ny),6,dx,G,K[p,2],K[p+nx*ny,6])
        
        
    if(voxel_db[p,5].mat == 0):               #N4
        
        UA_s1 = side_diffmat(p,2,p,3,dx,G,K[p,2],HTC)
        
    if(voxel_db[p,5].mat == mat):
        
        UA_s1 = side_samemat(p,2,p,3,dx,G,K[p,2],K[p,3])
        
        
    if(voxel_db[p,3].mat == 0):               #N2
        
        UA_s2 = side_diffmat(p,2,p,1,dx,G,K[p,2],HTC)
        
    if(voxel_db[p,3].mat == mat):
        
        UA_s2 = side_samemat(p,2,p,1,dx,G,K[p,2],K[p,1])
        
        
    if(voxel_db[p,18].mat == 0):              #F1
        
        UA_diag = diag_diffmat(p,2,p,16,dx,G,K[p,2],HTC)
        
    if(voxel_db[p,18].mat == mat):
        
        UA_diag = diag_samemat(p,2,p,16,dx,G,K[p,2],K[p,16])
    
    i_cor.append(p*24+2) ; j_cor.append((p+nx*ny)*24+6) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+2) ; j_cor.append(p*24+3) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+2) ; j_cor.append(p*24+1) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+2) ; j_cor.append(p*24+16) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+2) ; j_cor.append(p*24+2) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    

def N4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("N4")

    if(voxel_db[p + nx*ny,9].mat == 0):       #S4
        
        UA_perp = perp_diffmat(p,3,p+nx*ny,7,dx,G,K[p,3],HTC)
        
    if(voxel_db[p + nx*ny,9].mat == mat):
        
        UA_perp = perp_samemat(p,3,p+nx*ny,7,dx,G,K[p,3],K[p+nx*ny,7])
        
        
    if(voxel_db[p,2].mat == 0):              #N1
        
        UA_s1 = side_diffmat(p,3,p,0,dx,G,K[p,3],HTC)
        
    if(voxel_db[p,2].mat == mat): 
        
        UA_s1 = side_samemat(p,3,p,0,dx,G,K[p,3],K[p,0])
        
        
    if(voxel_db[p,4].mat == 0):              #N3
        
        UA_s2 = side_diffmat(p,3,p,2,dx,G,K[p,3],HTC)
        
    if(voxel_db[p,4].mat == mat):
        
        UA_s2 = side_samemat(p,3,p,2,dx,G,K[p,3],K[p,2])
        
    if(voxel_db[p,14].mat == 0):             #E1
        
        UA_diag = diag_diffmat(p,3,p,12,dx,G,K[p,3],HTC)
        
    if(voxel_db[p,14].mat == mat):
        
        UA_diag = diag_samemat(p,3,p,12,dx,G,K[p,3],K[p,12])
    
    i_cor.append(p*24+3) ; j_cor.append((p+nx*ny)*24+7) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+3) ; j_cor.append(p*24+0) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+3) ; j_cor.append(p*24+2) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+3) ; j_cor.append(p*24+12) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+3) ; j_cor.append(p*24+3) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
# South Tetras #######################################################
    
def S1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("S1")
    if(voxel_db[p - nx*ny,2].mat == 0):       #N1
        
        UA_perp = perp_diffmat(p,4,p - nx*ny,0,dx,G,K[p,4],HTC)
        
    if(voxel_db[p - nx*ny,2].mat == mat):       #N1
        
        UA_perp = perp_samemat(p,4,p-nx*ny,0,dx,G,K[p,4],K[p-nx*ny,0])
        
    if(voxel_db[p,9].mat == 0):              #S4
        
        UA_s1 = side_diffmat(p,4,p,7,dx,G,K[p,4],HTC)
        
    if(voxel_db[p,9].mat == mat):              #S4
        
        UA_s1 = side_samemat(p,4,p,7,dx,G,K[p,4],K[p,7])       
        
    if(voxel_db[p,7].mat == 0):               #S2
        
        UA_s2 = side_diffmat(p,4,p,5,dx,G,K[p,4],HTC)
        
    if(voxel_db[p,7].mat == mat):               #S2
        
        UA_s2 = side_samemat(p,4,p,5,dx,G,K[p,4],K[p,5])
        
    if(voxel_db[p,24].mat == 0):             #B3
        
        UA_diag = diag_diffmat(p,4,p,22,dx,G,K[p,4],HTC)
        
    if(voxel_db[p,24].mat == mat):             #B3
        
        UA_diag = diag_samemat(p,4,p,22,dx,G,K[p,4],K[p,22])
        
    i_cor.append(p*24+4) ; j_cor.append((p-nx*ny)*24+0) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+4) ; j_cor.append(p*24+7) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+4) ; j_cor.append(p*24+5) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+4) ; j_cor.append(p*24+22) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+4) ; j_cor.append(p*24+4) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)

def S2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("S2")
    if(voxel_db[p - nx*ny,3].mat != mat):       #N2
        
        UA_perp = perp_diffmat(p,5,p-nx*ny,1,dx,G,K[p,5],HTC)
        
    if(voxel_db[p - nx*ny,3].mat == mat):       #N2
        
        UA_perp = perp_samemat(p,5,p-nx*ny,1,dx,G,K[p,5],K[p-nx*ny,1])
        
    if(voxel_db[p,6].mat != mat):                #S1
        
        UA_s1 = side_diffmat(p,5,p,4,dx,G,K[p,5],HTC)
        
    if(voxel_db[p,6].mat == mat):                #S1
        
        UA_s1 = side_samemat(p,5,p,4,dx,G,K[p,5],K[p,4])
        
    if(voxel_db[p,8].mat != mat):                #S3
        
        UA_s2 = side_diffmat(p,5,p,6,dx,G,K[p,5],HTC)
        
    if(voxel_db[p,8].mat == mat):                #S3
        
        UA_s2 = side_samemat(p,5,p,6,dx,G,K[p,5],K[p,6])
        
    if(voxel_db[p,12].mat != mat):              #W3
        
        UA_diag = diag_diffmat(p,5,p,10,dx,G,K[p,5],HTC)
        
    if(voxel_db[p,12].mat == mat):              #W3
        
        UA_diag = diag_samemat(p,5,p,10,dx,G,K[p,5],K[p,10])
   
    i_cor.append(p*24+5) ; j_cor.append((p-nx*ny)*24+1) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+5) ; j_cor.append(p*24+4) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+5) ; j_cor.append(p*24+6) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+5) ; j_cor.append(p*24+10) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+5) ; j_cor.append(p*24+5) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
    
def S3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("S3")
    if(voxel_db[p - nx*ny,4].mat != mat):        #N3
        
        UA_perp = perp_diffmat(p,6,p-nx*ny,2,dx,G,K[p,6],HTC)
        
    if(voxel_db[p - nx*ny,4].mat == mat):        #N3
        
        UA_perp = perp_samemat(p,6,p-nx*ny,2,dx,G,K[p,6],K[p-nx*ny,2])
        
    if(voxel_db[p,9].mat != mat):               #S4
        
        UA_s1 = side_diffmat(p,6,p,7,dx,G,K[p,6],HTC)
        
    if(voxel_db[p,9].mat == mat):               #S4
        
        UA_s1 = side_samemat(p,6,p,7,dx,G,K[p,6],K[p,7])
        
    if(voxel_db[p,7].mat != mat):               #S2
        
        UA_s2 = side_diffmat(p,6,p,5,dx,G,K[p,6],HTC)
        
    if(voxel_db[p,7].mat == mat):               #S2
        
        UA_s2 = side_samemat(p,6,p,5,dx,G,K[p,6],K[p,5])
        
    if(voxel_db[p,20].mat != mat):               #F3
        
        UA_diag = diag_diffmat(p,6,p,18,dx,G,K[p,6],HTC)
        
    if(voxel_db[p,20].mat == mat):               #F3
        
        UA_diag = diag_samemat(p,6,p,18,dx,G,K[p,6],K[p,18])
    
    i_cor.append(p*24+6) ; j_cor.append((p-nx*ny)*24+2) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+6) ; j_cor.append(p*24+7) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+6) ; j_cor.append(p*24+5) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+6) ; j_cor.append(p*24+18) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+6) ; j_cor.append(p*24+6) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
    

def S4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("S4")
    if(voxel_db[p - nx*ny,5].mat != mat):        #N4
        
        UA_perp = perp_diffmat(p,7,p-nx*ny,3,dx,G,K[p,7],HTC)
        
    if(voxel_db[p - nx*ny,5].mat == mat):        #N4
        
        UA_perp = perp_samemat(p,7,p-nx*ny,3,dx,G,K[p,7],K[p-nx*ny,3])
        
    if(voxel_db[p,6].mat != mat):               #S1
        
        UA_s1 = side_diffmat(p,7,p,4,dx,G,K[p,7],HTC)
        
    if(voxel_db[p,6].mat == mat):               #S1
        
        UA_s1 = side_samemat(p,7,p,4,dx,G,K[p,7],K[p,4])
        
    if(voxel_db[p,8].mat != mat):                #S3
        
        UA_s2 = side_diffmat(p,7,p,6,dx,G,K[p,7],HTC)
        
    if(voxel_db[p,8].mat == mat):                #S3
        
        UA_s2 = side_samemat(p,7,p,6,dx,G,K[p,7],K[p,6])
        
    if(voxel_db[p,16].mat != mat):               #E3
        
        UA_diag = diag_diffmat(p,7,p,14,dx,G,K[p,7],HTC)
        
    if(voxel_db[p,16].mat == mat):               #E3
        
        UA_diag = diag_samemat(p,7,p,14,dx,G,K[p,7],K[p,14])
    
    i_cor.append(p*24+7) ; j_cor.append((p-nx*ny)*24+3) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+7) ; j_cor.append(p*24+4) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+7) ; j_cor.append(p*24+6) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+7) ; j_cor.append(p*24+14) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+7) ; j_cor.append(p*24+7) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
   

# West Tetras #############################################################    

def W1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("W1")
    if(voxel_db[p - ny,14].mat != mat):         #E1
        
        UA_perp = perp_diffmat(p,8,p-ny,12,dx,G,K[p,8],HTC)
        
    if(voxel_db[p - ny,14].mat == mat):         #E1
        
        UA_perp = perp_samemat(p,8,p-ny,12,dx,G,K[p,8],K[p-ny,12])
        
    if(voxel_db[p,11].mat != mat):              #W2
        
        UA_s1 = side_diffmat(p,8,p,9,dx,G,K[p,8],HTC)
        
    if(voxel_db[p,11].mat == mat):              #W2
        
        UA_s1 = side_samemat(p,8,p,9,dx,G,K[p,8],K[p,9])
        
    if(voxel_db[p,13].mat != mat):              #W4
        
        UA_s2 = side_diffmat(p,8,p,11,dx,G,K[p,8],HTC)
        
    if(voxel_db[p,13].mat == mat):              #W4
        
        UA_s2 = side_samemat(p,8,p,11,dx,G,K[p,8],K[p,11])
        
    if(voxel_db[p,3].mat != mat):             #N2
        
        UA_diag = diag_diffmat(p,8,p,1,dx,G,K[p,8],HTC)
        
    if(voxel_db[p,3].mat == mat):             #N2
        
        UA_diag = diag_samemat(p,8,p,1,dx,G,K[p,8],K[p,1])
    
    i_cor.append(p*24+8) ; j_cor.append((p-ny)*24+12) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+8) ; j_cor.append(p*24+9) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+8) ; j_cor.append(p*24+11) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+8) ; j_cor.append(p*24+1) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+8) ; j_cor.append(p*24+8) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
    
def W2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("W2")
    if(voxel_db[p - ny,17].mat != mat):          #E4
        
        UA_perp = perp_diffmat(p,9,p-ny,15,dx,G,K[p,9],HTC)
        
    if(voxel_db[p - ny,17].mat == mat):          #E4
        
        UA_perp = perp_samemat(p,9,p-ny,15,dx,G,K[p,9],K[p-ny,15])
        
    if(voxel_db[p,12].mat != mat):              #W3
        
        UA_s1 = side_diffmat(p,9,p,10,dx,G,K[p,9],HTC)
        
    if(voxel_db[p,12].mat == mat):              #W3
        
        UA_s1 = side_samemat(p,9,p,10,dx,G,K[p,9],K[p,10])
        
    if(voxel_db[p,10].mat != mat):               #W1
        
        UA_s2 = side_diffmat(p,9,p,8,dx,G,K[p,9],HTC)
        
    if(voxel_db[p,10].mat == mat):               #W1
        
        UA_s2 = side_samemat(p,9,p,8,dx,G,K[p,9],K[p,8])
        
    if(voxel_db[p,25].mat != mat):               #B4
        
        UA_diag = diag_diffmat(p,9,p,23,dx,G,K[p,9],HTC)
        
    if(voxel_db[p,25].mat == mat):               #B4
        
        UA_diag = diag_samemat(p,9,p,23,dx,G,K[p,9],K[p,23])
    
    i_cor.append(p*24+9) ; j_cor.append((p-ny)*24+15) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+9) ; j_cor.append(p*24+10) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+9) ; j_cor.append(p*24+8) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+9) ; j_cor.append(p*24+23) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+9) ; j_cor.append(p*24+9) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
   

def W3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("W3")
    if(voxel_db[p - ny,16].mat != mat):          #E3
        
        UA_perp = perp_diffmat(p,10,p-ny,14,dx,G,K[p,10],HTC)
    
    if(voxel_db[p - ny,16].mat == mat):          #E3
        
        UA_perp = perp_samemat(p,10,p-ny,14,dx,G,K[p,10],K[p-ny,14])
        
    if(voxel_db[p,13].mat != mat):              #W4
        
        UA_s1 = side_diffmat(p,10,p,11,dx,G,K[p,10],HTC)
        
    if(voxel_db[p,13].mat == mat):              #W4
        
        UA_s1 = side_samemat(p,10,p,11,dx,G,K[p,10],K[p,11])
        
    if(voxel_db[p,11].mat != mat):               #W2
        
        UA_s2 = side_diffmat(p,10,p,9,dx,G,K[p,10],HTC)
        
    if(voxel_db[p,11].mat == mat):               #W2
        
        UA_s2 = side_samemat(p,10,p,9,dx,G,K[p,10],K[p,9])
        
    if(voxel_db[p,7].mat != mat):                #S2
        
        UA_diag = diag_diffmat(p,10,p,5,dx,G,K[p,10],HTC)
        
    if(voxel_db[p,7].mat == mat):                #S2
        
        UA_diag = diag_samemat(p,10,p,5,dx,G,K[p,10],K[p,5])
    
    i_cor.append(p*24+10) ; j_cor.append((p-ny)*24+14) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+10) ; j_cor.append(p*24+11) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+10) ; j_cor.append(p*24+9) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+10) ; j_cor.append(p*24+5) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+10) ; j_cor.append(p*24+10) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
def W4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("W4")
    if(voxel_db[p - ny,15].mat != mat):          #E2
        
        UA_perp = perp_diffmat(p,11,p-ny,13,dx,G,K[p,11],HTC)
        
    if(voxel_db[p - ny,15].mat == mat):          #E2
        
        UA_perp = perp_samemat(p,11,p-ny,13,dx,G,K[p,11],K[p-ny,13])
        
    if(voxel_db[p,10].mat != mat):               #W1
        
        UA_s1 = side_diffmat(p,11,p,8,dx,G,K[p,11],HTC)
        
    if(voxel_db[p,10].mat == mat):               #W1
        
        UA_s1 = side_samemat(p,11,p,8,dx,G,K[p,11],K[p,8])
        
    if(voxel_db[p,12].mat != mat):               #W3
        
        UA_s2 = side_diffmat(p,11,p,10,dx,G,K[p,11],HTC)
        
    if(voxel_db[p,12].mat == mat):               #W3
        
        UA_s2 = side_samemat(p,11,p,10,dx,G,K[p,11],K[p,10])
        
    if(voxel_db[p,19].mat != mat):               #F2
        
        UA_diag = diag_diffmat(p,11,p,17,dx,G,K[p,11],HTC)
        
    if(voxel_db[p,19].mat == mat):               #F2
        
        UA_diag = diag_samemat(p,11,p,17,dx,G,K[p,11],K[p,17])
    
    i_cor.append(p*24+11) ; j_cor.append((p-ny)*24+13) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+11) ; j_cor.append(p*24+8) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+11) ; j_cor.append(p*24+10) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+11) ; j_cor.append(p*24+17) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+11) ; j_cor.append(p*24+11) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
        
   

# East Tetras ###########################################################

def E1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("E1")
    if(voxel_db[p + ny,10].mat != mat):         #W1
        
        UA_perp = perp_diffmat(p,12,p+ny,8,dx,G,K[p,12],HTC)
        
    if(voxel_db[p + ny,10].mat == mat):         #W1
        
        UA_perp = perp_samemat(p,12,p+ny,8,dx,G,K[p,12],K[p+ny,8])
        
    if(voxel_db[p,15].mat != mat):               #E2
        
        UA_s1 = side_diffmat(p,12,p,13,dx,G,K[p,12],HTC)
        
    if(voxel_db[p,15].mat == mat):               #E2
        
        UA_s1 = side_samemat(p,12,p,13,dx,G,K[p,12],K[p,13])
        
    if(voxel_db[p,17].mat != mat):              #E4
        
        UA_s2 = side_diffmat(p,12,p,15,dx,G,K[p,12],HTC)
        
    if(voxel_db[p,17].mat == mat):              #E4
        
        UA_s2 = side_samemat(p,12,p,15,dx,G,K[p,12],K[p,15])
        
    if(voxel_db[p,5].mat != mat):               #N4
        
        UA_diag = diag_diffmat(p,12,p,3,dx,G,K[p,12],HTC)
        
    if(voxel_db[p,5].mat == mat):               #N4
        
        UA_diag = diag_samemat(p,12,p,3,dx,G,K[p,12],K[p,3])

    i_cor.append(p*24+12) ; j_cor.append((p+ny)*24+8) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+12) ; j_cor.append(p*24+13) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+12) ; j_cor.append(p*24+15) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+12) ; j_cor.append(p*24+3) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+12) ; j_cor.append(p*24+12) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
def E2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("E2")
    if(voxel_db[p + ny,13].mat != mat):         #W4
        
        UA_perp = perp_diffmat(p,13,p+ny,11,dx,G,K[p,13],HTC)
        
    if(voxel_db[p + ny,13].mat == mat):         #W4
        
        UA_perp = perp_samemat(p,13,p+ny,11,dx,G,K[p,13],K[p+ny,11])
        
    if(voxel_db[p,16].mat != mat):              #E3
        
        UA_s1 = side_diffmat(p,13,p,14,dx,G,K[p,13],HTC)
        
    if(voxel_db[p,16].mat == mat):              #E3
        
        UA_s1 = side_samemat(p,13,p,14,dx,G,K[p,13],K[p,14])
        
    if(voxel_db[p,14].mat != mat):              #E1
        
        UA_s2 = side_diffmat(p,13,p,12,dx,G,K[p,13],HTC)
        
    if(voxel_db[p,14].mat == mat):              #E1
        
        UA_s2 = side_samemat(p,13,p,12,dx,G,K[p,13],K[p,12])
                
    if(voxel_db[p,21].mat != mat):              #F4
        
        UA_diag = diag_diffmat(p,13,p,19,dx,G,K[p,13],HTC)
        
    if(voxel_db[p,21].mat == mat):              #F4
        
        UA_diag = diag_samemat(p,13,p,19,dx,G,K[p,13],K[p,19])
    
    i_cor.append(p*24+13) ; j_cor.append((p+ny)*24+11) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+13) ; j_cor.append(p*24+14) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+13) ; j_cor.append(p*24+12) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+13) ; j_cor.append(p*24+19) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+13) ; j_cor.append(p*24+13) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
        
        
def E3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("E3")
    if(voxel_db[p + ny,12].mat != mat):         #W3
        
        UA_perp = perp_diffmat(p,14,p+ny,10,dx,G,K[p,14],HTC)
        
    if(voxel_db[p + ny,12].mat == mat):         #W3
        
        UA_perp = perp_samemat(p,14,p+ny,10,dx,G,K[p,14],K[p+ny,10])
        
    if(voxel_db[p,17].mat != mat):              #E4
        
        UA_s1 = side_diffmat(p,14,p,15,dx,G,K[p,14],HTC)
        
    if(voxel_db[p,17].mat == mat):              #E4
        
        UA_s1 = side_samemat(p,14,p,15,dx,G,K[p,14],K[p,15])
        
    if(voxel_db[p,15].mat != mat):              #E2
        
        UA_s2 = side_diffmat(p,14,p,13,dx,G,K[p,14],HTC)
        
    if(voxel_db[p,15].mat == mat):              #E2
        
        UA_s2 = side_samemat(p,14,p,13,dx,G,K[p,14],K[p,13])
        
    if(voxel_db[p,9].mat != mat):               #S4
        
        UA_diag = diag_diffmat(p,14,p,7,dx,G,K[p,14],HTC)
        
    if(voxel_db[p,9].mat == mat):               #S4
        
        UA_diag = diag_samemat(p,14,p,7,dx,G,K[p,14],K[p,7])
    
    i_cor.append((p*24+14)) ; j_cor.append((p+ny)*24+10) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append((p*24+14)) ; j_cor.append(p*24+15) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append((p*24+14)) ; j_cor.append(p*24+13) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append((p*24+14)) ; j_cor.append(p*24+7) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append((p*24+14)) ; j_cor.append((p*24+14)) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
    
def E4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("E4")
    if(voxel_db[p + ny,11].mat != mat):         #W2
       
        UA_perp = perp_diffmat(p,15,p+ny,9,dx,G,K[p,15],HTC)
        
    if(voxel_db[p + ny,11].mat == mat):         #W2
       
        UA_perp = perp_samemat(p,15,p+ny,9,dx,G,K[p,15],K[p+ny,9])
        
    if(voxel_db[p,14].mat != mat):              #E1
        
        UA_s1 = side_diffmat(p,15,p,12,dx,G,K[p,15],HTC)
        
    if(voxel_db[p,14].mat == mat):              #E1
        
        UA_s1 = side_samemat(p,15,p,12,dx,G,K[p,15],K[p,12])
        
    if(voxel_db[p,16].mat != mat):              #E3
        
        UA_s2 = side_diffmat(p,15,p,14,dx,G,K[p,15],HTC)
        
    if(voxel_db[p,16].mat == mat):              #E3
        
        UA_s2 = side_samemat(p,15,p,14,dx,G,K[p,15],K[p,14])
        
    if(voxel_db[p,23].mat != mat):              #B2
        
        UA_diag = diag_diffmat(p,15,p,21,dx,G,K[p,15],HTC)
        
    if(voxel_db[p,23].mat == mat):              #B2
        
        UA_diag = diag_samemat(p,15,p,21,dx,G,K[p,15],K[p,21])
    
    i_cor.append(p*24+15) ; j_cor.append((p+ny)*24+9) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+15) ; j_cor.append(p*24+12) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+15) ; j_cor.append(p*24+14) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+15) ; j_cor.append(p*24+21) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+15) ; j_cor.append(p*24+15) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
    
        
# Front Tetras ############################################################

def F1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("F1")
    if(voxel_db[p + 1,22].mat != mat):          #B1
        
        UA_perp = perp_diffmat(p,16,p+1,20,dx,G,K[p,16],HTC)
        
    if(voxel_db[p + 1,22].mat == mat):          #B1
        
        UA_perp = perp_samemat(p,16,p+1,20,dx,G,K[p,16],K[p+1,20])
        
    if(voxel_db[p,19].mat != mat):              #F2
        
        UA_s1 = side_diffmat(p,16,p,17,dx,G,K[p,16],HTC)
        
    if(voxel_db[p,19].mat == mat):              #F2
        
        UA_s1 = side_samemat(p,16,p,17,dx,G,K[p,16],K[p,17])
        
    if(voxel_db[p,21].mat != mat):              #F4
        
        UA_s2 = side_diffmat(p,16,p,19,dx,G,K[p,16],HTC)
        
    if(voxel_db[p,21].mat == mat):              #F4
        
        UA_s2 = side_samemat(p,16,p,19,dx,G,K[p,16],K[p,19])
        
    if(voxel_db[p,4].mat != mat):               #N3
        
        UA_diag = diag_diffmat(p,16,p,2,dx,G,K[p,16],HTC)
        
    if(voxel_db[p,4].mat == mat):               #N3
        
        UA_diag = diag_samemat(p,16,p,2,dx,G,K[p,16],K[p,2])
    
    i_cor.append(p*24+16) ; j_cor.append((p+1)*24+20) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+16) ; j_cor.append(p*24+17) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+16) ; j_cor.append(p*24+19) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+16) ; j_cor.append(p*24+2) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+16) ; j_cor.append(p*24+16) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    

def F2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("F2")
    if(voxel_db[p + 1,25].mat != mat):          #B4
        
        UA_perp = perp_diffmat(p,17,p+1,23,dx,G,K[p,17],HTC)
        
    if(voxel_db[p + 1,25].mat == mat):          #B4
        
        UA_perp = perp_samemat(p,17,p+1,23,dx,G,K[p,17],K[p+1,23])
        
    if(voxel_db[p,20].mat != mat):              #F3
        
        UA_s1 = side_diffmat(p,17,p,18,dx,G,K[p,17],HTC)
        
    if(voxel_db[p,20].mat == mat):              #F3
        
        UA_s1 = side_samemat(p,17,p,18,dx,G,K[p,17],K[p,18])
        
    if(voxel_db[p,18].mat != mat):              #F1
        
        UA_s2 = side_diffmat(p,17,p,16,dx,G,K[p,17],HTC)
        
    if(voxel_db[p,18].mat == mat):              #F1
        
        UA_s2 = side_samemat(p,17,p,16,dx,G,K[p,17],K[p,16])
        
    if(voxel_db[p,13].mat != mat):              #W4
        
        UA_diag = diag_diffmat(p,17,p,11,dx,G,K[p,17],HTC)
        
    if(voxel_db[p,13].mat == mat):              #W4
        
        UA_diag = diag_samemat(p,17,p,11,dx,G,K[p,17],K[p,11])
     
    i_cor.append(p*24+17) ; j_cor.append((p+1)*24+23) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+17) ; j_cor.append(p*24+18) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+17) ; j_cor.append(p*24+16) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+17) ; j_cor.append(p*24+11) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+17) ; j_cor.append(p*24+17) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
            
def F3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("F3")
    if(voxel_db[p + 1,24].mat != mat):          #B3
        
        UA_perp = perp_diffmat(p,18,p+1,22,dx,G,K[p,18],HTC)
        
    if(voxel_db[p + 1,24].mat == mat):          #B3
        
        UA_perp = perp_samemat(p,18,p+1,22,dx,G,K[p,18],K[p+1,22])
        
    if(voxel_db[p,21].mat != mat):              #F4
        
        UA_s1 = side_diffmat(p,18,p,19,dx,G,K[p,18],HTC)
        
    if(voxel_db[p,21].mat == mat):              #F4
        
        UA_s1 = side_samemat(p,18,p,19,dx,G,K[p,18],K[p,19])
        
    if(voxel_db[p,19].mat != mat):              #F2
        
        UA_s2 = side_diffmat(p,18,p,17,dx,G,K[p,18],HTC)
        
    if(voxel_db[p,19].mat == mat):              #F2
        
        UA_s2 = side_samemat(p,18,p,17,dx,G,K[p,18],K[p,17])
        
    if(voxel_db[p,8].mat != mat):               #S3
        
        UA_diag = diag_diffmat(p,18,p,6,dx,G,K[p,18],HTC)
        
    if(voxel_db[p,8].mat == mat):               #S3
        
        UA_diag = diag_samemat(p,18,p,6,dx,G,K[p,18],K[p,6])
    
    i_cor.append(p*24+18) ; j_cor.append((p+1)*24+22) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+18) ; j_cor.append(p*24+19) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+18) ; j_cor.append(p*24+17) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+18) ; j_cor.append(p*24+6) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+18) ; j_cor.append(p*24+18) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
            
def F4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("F4")
    if(voxel_db[p + 1,23].mat != mat):          #B2
        
        UA_perp = perp_diffmat(p,19,p+1,21,dx,G,K[p,19],HTC)
        
    if(voxel_db[p + 1,23].mat == mat):          #B2
        
        UA_perp = perp_samemat(p,19,p+1,21,dx,G,K[p,19],K[p+1,21])
        
    if(voxel_db[p,18].mat != mat):              #F1
        
        UA_s1 = side_diffmat(p,19,p,16,dx,G,K[p,19],HTC)
        
    if(voxel_db[p,18].mat == mat):              #F1
        
        UA_s1 = side_samemat(p,19,p,16,dx,G,K[p,19],K[p,16])
        
    if(voxel_db[p,20].mat != mat):              #F3
        
        UA_s2 = side_diffmat(p,19,p,18,dx,G,K[p,19],HTC)
        
    if(voxel_db[p,20].mat == mat):              #F3
        
        UA_s2 = side_samemat(p,19,p,18,dx,G,K[p,19],K[p,18])
        
    if(voxel_db[p,15].mat != mat):
       
       UA_diag = diag_diffmat(p,19,p,13,dx,G,K[p,19],HTC)
       
    if(voxel_db[p,15].mat == mat):
       
       UA_diag = diag_samemat(p,19,p,13,dx,G,K[p,19],K[p,13])
    
    i_cor.append(p*24+19) ; j_cor.append((p+1)*24+21) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+19) ; j_cor.append(p*24+16) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+19) ; j_cor.append(p*24+18) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+19) ; j_cor.append(p*24+13) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+19) ; j_cor.append(p*24+19) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    

# Back Tetras #########################################################   
     
def B1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
 #   print("B1")
    if(voxel_db[p - 1,18].mat != mat):          #F1
        
        UA_perp = perp_diffmat(p,20,p-1,16,dx,G,K[p,20],HTC)
        
    if(voxel_db[p - 1,18].mat == mat):          #F1
        
        UA_perp = perp_samemat(p,20,p-1,16,dx,G,K[p,20],K[p-1,16])
        
    if(voxel_db[p,23].mat != mat):              #B2
        
        UA_s1 = side_diffmat(p,20,p,21,dx,G,K[p,20],HTC)
        
    if(voxel_db[p,23].mat == mat):              #B2
        
        UA_s1 = side_samemat(p,20,p,21,dx,G,K[p,20],K[p,21])
        
    if(voxel_db[p,25].mat != mat):              #B4
        
        UA_s2 = side_diffmat(p,20,p,23,dx,G,K[p,20],HTC)
        
    if(voxel_db[p,25].mat == mat):              #B4
        
        UA_s2 = side_samemat(p,20,p,23,dx,G,K[p,20],K[p,23])
        
    if(voxel_db[p,2].mat != mat):               #N1
        
        UA_diag = diag_diffmat(p,20,p,0,dx,G,K[p,20],HTC)
        
    if(voxel_db[p,2].mat == mat):               #N1
        
        UA_diag = diag_samemat(p,20,p,0,dx,G,K[p,20],K[p,0])
    
    i_cor.append(p*24+20) ; j_cor.append((p-1)*24+16) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+20) ; j_cor.append(p*24+21) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+20) ; j_cor.append(p*24+23) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+20) ; j_cor.append(p*24+0) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+20) ; j_cor.append(p*24+20) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    

def B2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("B2")
    if(voxel_db[p - 1,21].mat != mat):          #F4
        
        UA_perp = perp_diffmat(p,21,p-1,19,dx,G,K[p,21],HTC)
        
    if(voxel_db[p - 1,21].mat == mat):          #F4
        
        UA_perp = perp_samemat(p,21,p-1,19,dx,G,K[p,21],K[p-1,19])
        
    if(voxel_db[p,24].mat != mat):              #B3
        
        UA_s1 = side_diffmat(p,21,p,22,dx,G,K[p,21],HTC)
        
    if(voxel_db[p,24].mat == mat):              #B3
        
        UA_s1 = side_samemat(p,21,p,22,dx,G,K[p,21],K[p,22])
        
    if(voxel_db[p,22].mat != mat):              #B1
        
        UA_s2 = side_diffmat(p,21,p,20,dx,G,K[p,21],HTC)
        
    if(voxel_db[p,22].mat == mat):              #B1
        
        UA_s2 = side_samemat(p,21,p,20,dx,G,K[p,21],K[p,20])
        
    if(voxel_db[p,17].mat != mat):              #E4
        
        UA_diag = diag_diffmat(p,21,p,15,dx,G,K[p,21],HTC)
        
    if(voxel_db[p,17].mat == mat):              #E4
        
        UA_diag = diag_samemat(p,21,p,15,dx,G,K[p,21],K[p,15])
    
    i_cor.append(p*24+21) ; j_cor.append((p-1)*24+19) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+21) ; j_cor.append(p*24+22) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+21) ; j_cor.append(p*24+20) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+21) ; j_cor.append(p*24+15) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+21) ; j_cor.append(p*24+21) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    

def B3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("B3")
    if(voxel_db[p - 1,20].mat != mat):          #F3
        
        UA_perp = perp_diffmat(p,22,p-1,18,dx,G,K[p,22],HTC)
        
    if(voxel_db[p - 1,20].mat == mat):          #F3
        
        UA_perp = perp_samemat(p,22,p-1,18,dx,G,K[p,22],K[p-1,18])
        
    if(voxel_db[p,25].mat != mat):              #B4
        
        UA_s1 = side_diffmat(p,22,p,23,dx,G,K[p,22],HTC)
        
    if(voxel_db[p,25].mat == mat):              #B4
        
        UA_s1 = side_samemat(p,22,p,23,dx,G,K[p,22],K[p,23])
        
    if(voxel_db[p,23].mat != mat):              #B2
        
        UA_s2 = side_diffmat(p,22,p,21,dx,G,K[p,22],HTC)
        
    if(voxel_db[p,23].mat == mat):              #B2
        
        UA_s2 = side_samemat(p,22,p,21,dx,G,K[p,22],K[p,21])
        
    if(voxel_db[p,6].mat != mat):               #S3
        
        UA_diag = diag_diffmat(p,22,p,4,dx,G,K[p,22],HTC)
        
    if(voxel_db[p,6].mat == mat):               #S3
        
        UA_diag = diag_samemat(p,22,p,4,dx,G,K[p,22],K[p,4])
    
    i_cor.append(p*24+22) ; j_cor.append((p-1)*24+18) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+22) ; j_cor.append(p*24+23) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+22) ; j_cor.append(p*24+21) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+22) ; j_cor.append(p*24+4) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+22) ; j_cor.append(p*24+22) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    

def B4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
#    print("B4")
    if(voxel_db[p - 1,19].mat != mat):          #F2
        
        UA_perp = perp_diffmat(p,23,p-1,17,dx,G,K[p,23],HTC)
        
    if(voxel_db[p - 1,19].mat == mat):          #F2
        
        UA_perp = perp_samemat(p,23,p-1,17,dx,G,K[p,23],K[p-1,17])
        
    if(voxel_db[p,22].mat != mat):              #B1
        
        UA_s1 = side_diffmat(p,23,p,20,dx,G,K[p,23],HTC)
        
    if(voxel_db[p,22].mat == mat):              #B1
        
        UA_s1 = side_samemat(p,23,p,20,dx,G,K[p,23],K[p,20])
        
    if(voxel_db[p,24].mat != mat):              #B3
        
        UA_s2 = side_diffmat(p,23,p,22,dx,G,K[p,23],HTC)
        
    if(voxel_db[p,24].mat == mat):              #B3
        
        UA_s2 = side_samemat(p,23,p,22,dx,G,K[p,23],K[p,22])
        
    if(voxel_db[p,11].mat != mat):              #W2
        
        UA_diag = diag_diffmat(p,23,p,9,dx,G,K[p,23],HTC)
        
    if(voxel_db[p,11].mat == mat):              #W2
        
        UA_diag = diag_samemat(p,23,p,9,dx,G,K[p,23],K[p,9])
    
    i_cor.append(p*24+23) ; j_cor.append((p-1)*24+17) ; ua_data.append(UA_perp) # UA_perp
    i_cor.append(p*24+23) ; j_cor.append(p*24+20) ; ua_data.append(UA_s1) # UA_s1    
    i_cor.append(p*24+23) ; j_cor.append(p*24+22) ; ua_data.append(UA_s2) # UA_s2
    i_cor.append(p*24+23) ; j_cor.append(p*24+9) ; ua_data.append(UA_diag) # UA_diag
    i_cor.append(p*24+23) ; j_cor.append(p*24+23) ; ua_data.append(-(UA_perp + UA_s1 + UA_s2 + UA_diag)) # sum of UA
    Q.append(-q)
    
def func(tag,voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data):
    
    if(tag == 'N1'):
        (N1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'N2'):
        (N2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'N3'):
        (N3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'N4'):
        (N4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
 
    elif(tag == 'S1'):
        (S1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'S2'):
        (S2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'S3'):
        (S3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'S4'):
        (S4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))

    elif(tag == 'W1'):
        (W1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'W2'):
        (W2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'W3'):
        (W3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'W4'):
        (W4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
   
    elif(tag == 'E1'):
        (E1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'E2'):
        (E2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'E3'):
        (E3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'E4'):
        (E4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))

    elif(tag == 'F1'):
        (F1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'F2'):
        (F2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'F3'):
        (F3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'F4'):
        (F4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    
    elif(tag == 'B1'):
        (B1(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'B2'):
        (B2(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'B3'):
        (B3(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))
    elif(tag == 'B4'):
        (B4(voxel_db,p,nx,ny,nz,dx,mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data))

# UA matrix generator function

def matrixgenerator(voxel_db,K,HTC,dx,dy,dz,n_voxel,G,nx,ny,nz,q,Q,Tamb):
    i_cor = []
    j_cor = []
    ua_data = []
    Q = []
    vc = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                for t in range(24):
                    poi = 24*i + 24*ny*j + 24*ny*nz*k + t
                    if(voxel_db[vc,t+2].mat == 0):
                        i_cor.append(vc*24 + t)
                        j_cor.append(vc*24 + t)
                        ua_data.append(1.0)
                        Q.append(Tamb)
                        
                    if(voxel_db[vc,t+2].mat == 1):
                        loc = voxel_db[vc,(t+2)].pos
                        func(loc,voxel_db,vc,nx,ny,nz,dx,voxel_db[vc,(t+2)].mat,G,HTC,K,q,Q,poi,i_cor,j_cor,ua_data)
                        
                vc = vc + 1

    UA_csc = sp.csc_matrix((ua_data,(i_cor,j_cor)))
    print(len(i_cor),len(j_cor)); print("memory of i and j",sys.getsizeof(i_cor),sys.getsizeof(j_cor))
    return(UA_csc,np.asarray(Q))

