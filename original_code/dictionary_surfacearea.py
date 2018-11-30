

def N1(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("N1")
    if(voxel_db[p + nx*ny,6].mat != mat):       #S1
        sumarea = sumarea + dx*dx*0.25 # ; print ("perp S1",sumarea)
    if(voxel_db[p,3].mat != mat):               #N2
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 1 N2 ",p)
    if(voxel_db[p,5].mat != mat):               #N4
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 2 N4",p)
    if(voxel_db[p,22].mat != mat):              #B4
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side B4",p)
    return sumarea

def N2(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("N2")
    if(voxel_db[p + nx*ny,7].mat != mat):       #S2
        sumarea = sumarea + dx*dx*0.25 # ; print ("perp S2",sumarea)
    if(voxel_db[p,4].mat != mat):               #N3
        sumarea = sumarea + 0.17677*dx**2 # ; print ("side 1 N3",p)
    if(voxel_db[p,2].mat != mat):               #N1
        sumarea = sumarea + 0.17677*dx**2 # ; print ("side 2 N1",p)
    if(voxel_db[p,10].mat != mat):              #W1
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag W1",p)
    return sumarea

def N3(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("N3")
    if(voxel_db[p + nx*ny,8].mat != mat):       #S3
        sumarea = sumarea + dx*dx*0.25 # ; print("perp S3",sumarea)
    if(voxel_db[p,5].mat != mat):               #N4
       sumarea = sumarea + 0.17677*dx**2 # ; print("side 1 N4",p)
    if(voxel_db[p,3].mat != mat):               #N2
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 2 N2",p)
    if(voxel_db[p,18].mat != mat):              #F1
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag F1",p)
    return sumarea


def N4(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("N4")
    if(voxel_db[p + nx*ny,9].mat != mat):       #S4
        sumarea = sumarea + dx*dx*0.25 # ; print("perp S4",sumarea)
    if(voxel_db[p,2].mat != mat):              #N1
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 1 N1",p)
    if(voxel_db[p,4].mat != mat):              #N3
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 2 N3",p)
    if(voxel_db[p,14].mat != mat):             #E1
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag E1",p)
    return sumarea

    
def S1(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("S1")
    if(voxel_db[p - nx*ny,2].mat != mat):       #N1
        sumarea = sumarea + dx*dx*0.25 # ; print("perp N1",sumarea)
    if(voxel_db[p,9].mat != mat):              #S4
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 1 S4",p)
    if(voxel_db[p,7].mat != mat):               #S2
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 1 S2",p)
    if(voxel_db[p,24].mat != mat):             #B3
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side B3",p)
    return sumarea

def S2(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("S2")
    if(voxel_db[p - nx*ny,3].mat != mat):       #N2
        sumarea = sumarea + dx*dx*0.25 # ; print("perp N2",sumarea)
    if(voxel_db[p,9].mat != mat):                #S4
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 1 S4",p)
    if(voxel_db[p,8].mat != mat):                #S3
        sumarea = sumarea +  0.17677*dx**2 # ; print("side 2 S3",p)
    if(voxel_db[p,12].mat != mat):              #W3
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side W3",p)
    return sumarea
    
def S3(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("S3")
    if(voxel_db[p - nx*ny,4].mat != mat):        #N3
        sumarea = sumarea + dx*dx*0.25 # ; print("perp N3",sumarea)
    if(voxel_db[p,9].mat != mat):               #S4
        sumarea = sumarea +  0.17677*dx**2 # ; print("side 1 S4",p)
    if(voxel_db[p,7].mat != mat):               #S2
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 2 S2",p)
    if(voxel_db[p,20].mat != mat):               #F3
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side F3",p)
    return sumarea

def S4(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("S4")
    if(voxel_db[p - nx*ny,5].mat != mat):        #N4
        sumarea = sumarea + dx*dx*0.25 # ; print("perp N4",sumarea)
    if(voxel_db[p,6].mat != mat):               #S1
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 1 S1",p)
    if(voxel_db[p,8].mat != mat):                #S3
        sumarea = sumarea + 0.17677*dx**2 # ; print(" side 2 S3",p)
    if(voxel_db[p,16].mat != mat):               #E3
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side E3",p)
    return sumarea

def W1(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("W1")
    if(voxel_db[p - ny,14].mat != mat):         #E1
        sumarea = sumarea + dx*dx*0.25 # ; print("perp E1",sumarea)
    if(voxel_db[p,11].mat != mat):              #W2
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 1 W2")
    if(voxel_db[p,13].mat != mat):              #W4
        sumarea = sumarea + 0.17677*dx**2 # ; print("side 2 W4")
    if(voxel_db[p,3].mat != mat):             #N2
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side N2",p)
    return sumarea

def W2(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("W2")
    if(voxel_db[p - ny,17].mat != mat):          #E4
        sumarea = sumarea + dx*dx*0.25 # ; print("perp E4",sumarea)
    if(voxel_db[p,12].mat != mat):              #W3
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,10].mat != mat):               #W1
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,25].mat != mat):               #B4
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side B4",p)
    return sumarea

def W3(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("W3")
    if(voxel_db[p - ny,16].mat != mat):          #E3
        sumarea = sumarea + dx*dx*0.25 # ; print("perp E3",sumarea)
    if(voxel_db[p,13].mat != mat):              #W4
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,11].mat != mat):               #W2
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,7].mat != mat):                #S2
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side S2",p)
    return sumarea

def W4(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("W4")
    if(voxel_db[p - ny,15].mat != mat):          #E2
        sumarea = sumarea + dx*dx*0.25 # ; print("perp E2",sumarea)
    if(voxel_db[p,10].mat != mat):               #W1
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,12].mat != mat):               #W3
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,19].mat != mat):               #F2
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side F2",p)
    return sumarea

def E1(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("E1")
    if(voxel_db[p + ny,10].mat != mat):         #W1
        sumarea = sumarea + dx*dx*0.25 # ; print("perp W1",sumarea)
    if(voxel_db[p,15].mat != mat):               #E2
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,17].mat != mat):              #E4
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,5].mat != mat):               #N4
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side N4",p)
    return sumarea

def E2(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("E2")
    if(voxel_db[p + ny,13].mat != mat):         #W4
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp W4",sumarea)
    if(voxel_db[p,16].mat != mat):              #E3
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,14].mat != mat):              #E1
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,21].mat != mat):              #F4
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side F4",p)
    return sumarea

def E3(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("E3")
    if(voxel_db[p + ny,12].mat != mat):         #W3
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp W3",sumarea)
    if(voxel_db[p,17].mat != mat):              #E4
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,15].mat != mat):              #E2
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,9].mat != mat):               #S4
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side S4",p)
    return sumarea

def E4(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("E4")
    if(voxel_db[p + ny,11].mat != mat):         #W2
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp W2",sumarea)
    if(voxel_db[p,14].mat != mat):              #E1
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,16].mat != mat):              #E3
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,23].mat != mat):              #B2
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side B2",p)
    return sumarea

def F1(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("F1")
    if(voxel_db[p + 1,22].mat != mat):          #B1
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp B1",sumarea)
    if(voxel_db[p,19].mat != mat):              #F2
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,21].mat != mat):              #F4
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,4].mat != mat):               #N3
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side N3",p)
    return sumarea

def F2(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("F2")
    if(voxel_db[p + 1,25].mat != mat):          #B4
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp B4",sumarea)
    if(voxel_db[p,20].mat != mat):              #F3
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,18].mat != mat):              #F1
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,13].mat != mat):              #W4
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side W4",p)
    return sumarea

def F3(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("F3")
    if(voxel_db[p + 1,24].mat != mat):          #B3
        sumarea = sumarea + dx*dx*0.25 # ; print("perp B3",sumarea)
    if(voxel_db[p,21].mat != mat):              #F4
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,19].mat != mat):              #F2
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,8].mat != mat):               #S3
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side S3",p)
    return sumarea

def F4(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("F4")
    if(voxel_db[p + 1,23].mat != mat):          #B2
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp B2",sumarea)
    if(voxel_db[p,18].mat != mat):              #F1
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,20].mat != mat):              #F3
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,15].mat != mat):              #E2
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side E2",p)
    return sumarea

def B1(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("B1")
    if(voxel_db[p - 1,18].mat != mat):          #F1
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp F1",sumarea)
    if(voxel_db[p,23].mat != mat):              #B2
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,25].mat != mat):              #B4
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,2].mat != mat):               #N1
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side N1",p)
    return sumarea

def B2(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("B2")
    if(voxel_db[p - 1,21].mat != mat):          #F4
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp F4",sumarea)
    if(voxel_db[p,24].mat != mat):              #B3
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,22].mat != mat):              #B1
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,17].mat != mat):              #E4
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side E4",p)
    return sumarea

def B3(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("B3")
    if(voxel_db[p - 1,20].mat != mat):          #F3
        sumarea = sumarea + dx*dx*0.25 # ; print("perp F3",sumarea)
    if(voxel_db[p,25].mat != mat):              #B4
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,23].mat != mat):              #B2
        sumarea = sumarea +  0.17677*dx**2
    if(voxel_db[p,6].mat != mat):               #S3
        sumarea = sumarea + 0.3535*dx**2 # ; print("diag side S3",p)
    return sumarea

def B4(voxel_db,p,nx,ny,nz,dx,mat,sumarea):
#    print("B4")
    if(voxel_db[p - 1,19].mat != mat):          #F2
        sumarea = sumarea +  dx*dx*0.25 # ; print("perp F2",sumarea)
    if(voxel_db[p,22].mat != mat):              #B1
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,24].mat != mat):              #B3
        sumarea = sumarea + 0.17677*dx**2
    if(voxel_db[p,11].mat != mat):              #W2
        sumarea = sumarea +  0.3535*dx**2 # ; print("diag side w2",p)
    return sumarea

def func(tag,voxel_db,p,nx,ny,nz,dx,mat):
    sumarea = 0.0
    if(tag == 'N1'):
        sumarea = N1(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'N2'):
        sumarea = N2(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'N3'):
        sumarea = N3(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'N4'):
        sumarea = N4(voxel_db,p,nx,ny,nz,dx,mat,sumarea)        

    elif(tag == 'S1'):
        sumarea =  S1(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'S2'):
        sumarea =  S2(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'S3'):
        sumarea =  S3(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'S4'):
        sumarea =  S4(voxel_db,p,nx,ny,nz,dx,mat,sumarea)

    elif(tag == 'W1'):
        sumarea =  W1(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'W2'):
        sumarea =  W2(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'W3'):
        sumarea =  W3(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'W4'):
        sumarea =  W4(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    
    elif(tag == 'E1'):
        sumarea =  E1(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'E2'):
        sumarea = E2(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'E3'):
        sumarea = E3(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'E4'):
        sumarea = E4(voxel_db,p,nx,ny,nz,dx,mat,sumarea)

    elif(tag == 'F1'):
        sumarea = F1(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'F2'):
        sumarea = F2(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'F3'):
        sumarea = F3(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'F4'):
        sumarea = F4(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    
    elif(tag == 'B1'):
        sumarea = B1(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'B2'):
        sumarea = B2(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'B3'):
        sumarea = B3(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
    elif(tag == 'B4'):
        sumarea = B4(voxel_db,p,nx,ny,nz,dx,mat,sumarea)
  
    return sumarea

