def n1(voxel_db, p, nx, ny, dx, mat, sumarea):
    #    print("n1")
    if voxel_db[p + nx * ny, 6].mat != mat:  # s1
        sumarea = sumarea + dx * dx * 0.25  # ; print ("perp s1",sumarea)
    if voxel_db[p, 3].mat != mat:  # n2
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 n2 ",p)
    if voxel_db[p, 5].mat != mat:  # n4
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 2 n4",p)
    if voxel_db[p, 22].mat != mat:  # b4
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side b4",p)
    return sumarea


def n2(voxel_db, p, nx, ny, dx, mat, sumarea):
    #    print("n2")
    if voxel_db[p + nx * ny, 7].mat != mat:  # s2
        sumarea = sumarea + dx * dx * 0.25  # ; print ("perp s2",sumarea)
    if voxel_db[p, 4].mat != mat:  # n3
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print ("side 1 n3",p)
    if voxel_db[p, 2].mat != mat:  # n1
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print ("side 2 n1",p)
    if voxel_db[p, 10].mat != mat:  # w1
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag w1",p)
    return sumarea


def n3(voxel_db, p, nx, ny, dx, mat, sumarea):
    #    print("n3")
    if voxel_db[p + nx * ny, 8].mat != mat:  # s3
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp s3",sumarea)
    if voxel_db[p, 5].mat != mat:  # n4
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 n4",p)
    if voxel_db[p, 3].mat != mat:  # n2
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 2 n2",p)
    if voxel_db[p, 18].mat != mat:  # f1
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag f1",p)
    return sumarea


def n4(voxel_db, p, nx, ny, dx, mat, sumarea):
    #    print("n4")
    if voxel_db[p + nx * ny, 9].mat != mat:  # s4
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp s4",sumarea)
    if voxel_db[p, 2].mat != mat:  # n1
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 n1",p)
    if voxel_db[p, 4].mat != mat:  # n3
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 2 n3",p)
    if voxel_db[p, 14].mat != mat:  # e1
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag e1",p)
    return sumarea


def s1(voxel_db, p, nx, ny, dx, mat, sumarea):
    #    print("s1")
    if voxel_db[p - nx * ny, 2].mat != mat:  # n1
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp n1",sumarea)
    if voxel_db[p, 9].mat != mat:  # s4
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 s4",p)
    if voxel_db[p, 7].mat != mat:  # s2
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 s2",p)
    if voxel_db[p, 24].mat != mat:  # b3
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side b3",p)
    return sumarea


def s2(voxel_db, p, nx, ny, dx, mat, sumarea):
    #    print("s2")
    if voxel_db[p - nx * ny, 3].mat != mat:  # n2
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp n2",sumarea)
    if voxel_db[p, 9].mat != mat:  # s4
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 s4",p)
    if voxel_db[p, 8].mat != mat:  # s3
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 2 s3",p)
    if voxel_db[p, 12].mat != mat:  # w3
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side w3",p)
    return sumarea


def s3(voxel_db, p, nx, ny, dx, mat, sumarea):
    #    print("s3")
    if voxel_db[p - nx * ny, 4].mat != mat:  # n3
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp n3",sumarea)
    if voxel_db[p, 9].mat != mat:  # s4
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 s4",p)
    if voxel_db[p, 7].mat != mat:  # s2
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 2 s2",p)
    if voxel_db[p, 20].mat != mat:  # f3
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side f3",p)
    return sumarea


def s4(voxel_db, p, nx, ny, dx, mat, sumarea):
    #    print("s4")
    if voxel_db[p - nx * ny, 5].mat != mat:  # n4
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp n4",sumarea)
    if voxel_db[p, 6].mat != mat:  # s1
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 s1",p)
    if voxel_db[p, 8].mat != mat:  # s3
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print(" side 2 s3",p)
    if voxel_db[p, 16].mat != mat:  # e3
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side e3",p)
    return sumarea


def w1(voxel_db, p, ny, dx, mat, sumarea):
    #    print("w1")
    if voxel_db[p - ny, 14].mat != mat:  # e1
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp e1",sumarea)
    if voxel_db[p, 11].mat != mat:  # w2
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 1 w2")
    if voxel_db[p, 13].mat != mat:  # w4
        sumarea = sumarea + 0.17677 * dx ** 2  # ; print("side 2 w4")
    if voxel_db[p, 3].mat != mat:  # n2
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side n2",p)
    return sumarea


def w2(voxel_db, p, ny, dx, mat, sumarea):
    #    print("w2")
    if voxel_db[p - ny, 17].mat != mat:  # e4
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp e4",sumarea)
    if voxel_db[p, 12].mat != mat:  # w3
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 10].mat != mat:  # w1
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 25].mat != mat:  # b4
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side b4",p)
    return sumarea


def w3(voxel_db, p, ny, dx, mat, sumarea):
    #    print("w3")
    if voxel_db[p - ny, 16].mat != mat:  # e3
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp e3",sumarea)
    if voxel_db[p, 13].mat != mat:  # w4
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 11].mat != mat:  # w2
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 7].mat != mat:  # s2
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side s2",p)
    return sumarea


def w4(voxel_db, p, ny, dx, mat, sumarea):
    #    print("w4")
    if voxel_db[p - ny, 15].mat != mat:  # e2
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp e2",sumarea)
    if voxel_db[p, 10].mat != mat:  # w1
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 12].mat != mat:  # w3
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 19].mat != mat:  # f2
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side f2",p)
    return sumarea


def e1(voxel_db, p, ny, dx, mat, sumarea):
    #    print("e1")
    if voxel_db[p + ny, 10].mat != mat:  # w1
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp w1",sumarea)
    if voxel_db[p, 15].mat != mat:  # e2
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 17].mat != mat:  # e4
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 5].mat != mat:  # n4
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side n4",p)
    return sumarea


def e2(voxel_db, p, ny, dx, mat, sumarea):
    #    print("e2")
    if voxel_db[p + ny, 13].mat != mat:  # w4
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp w4",sumarea)
    if voxel_db[p, 16].mat != mat:  # e3
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 14].mat != mat:  # e1
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 21].mat != mat:  # f4
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side f4",p)
    return sumarea


def e3(voxel_db, p, ny, dx, mat, sumarea):
    #    print("e3")
    if voxel_db[p + ny, 12].mat != mat:  # w3
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp w3",sumarea)
    if voxel_db[p, 17].mat != mat:  # e4
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 15].mat != mat:  # e2
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 9].mat != mat:  # s4
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side s4",p)
    return sumarea


def e4(voxel_db, p, ny, dx, mat, sumarea):
    #    print("e4")
    if voxel_db[p + ny, 11].mat != mat:  # w2
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp w2",sumarea)
    if voxel_db[p, 14].mat != mat:  # e1
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 16].mat != mat:  # e3
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 23].mat != mat:  # b2
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side b2",p)
    return sumarea


def f1(voxel_db, p, dx, mat, sumarea):
    #    print("f1")
    if voxel_db[p + 1, 22].mat != mat:  # b1
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp b1",sumarea)
    if voxel_db[p, 19].mat != mat:  # f2
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 21].mat != mat:  # f4
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 4].mat != mat:  # n3
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side n3",p)
    return sumarea


def f2(voxel_db, p, dx, mat, sumarea):
    #    print("f2")
    if voxel_db[p + 1, 25].mat != mat:  # b4
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp b4",sumarea)
    if voxel_db[p, 20].mat != mat:  # f3
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 18].mat != mat:  # f1
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 13].mat != mat:  # w4
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side w4",p)
    return sumarea


def f3(voxel_db, p, dx, mat, sumarea):
    #    print("f3")
    if voxel_db[p + 1, 24].mat != mat:  # b3
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp b3",sumarea)
    if voxel_db[p, 21].mat != mat:  # f4
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 19].mat != mat:  # f2
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 8].mat != mat:  # s3
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side s3",p)
    return sumarea


def f4(voxel_db, p, dx, mat, sumarea):
    #    print("f4")
    if voxel_db[p + 1, 23].mat != mat:  # b2
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp b2",sumarea)
    if voxel_db[p, 18].mat != mat:  # f1
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 20].mat != mat:  # f3
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 15].mat != mat:  # e2
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side e2",p)
    return sumarea


def b1(voxel_db, p, dx, mat, sumarea):
    #    print("b1")
    if voxel_db[p - 1, 18].mat != mat:  # f1
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp f1",sumarea)
    if voxel_db[p, 23].mat != mat:  # b2
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 25].mat != mat:  # b4
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 2].mat != mat:  # n1
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side n1",p)
    return sumarea


def b2(voxel_db, p, dx, mat, sumarea):
    #    print("b2")
    if voxel_db[p - 1, 21].mat != mat:  # f4
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp f4",sumarea)
    if voxel_db[p, 24].mat != mat:  # b3
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 22].mat != mat:  # b1
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 17].mat != mat:  # e4
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side e4",p)
    return sumarea


def b3(voxel_db, p, dx, mat, sumarea):
    #    print("b3")
    if voxel_db[p - 1, 20].mat != mat:  # f3
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp f3",sumarea)
    if voxel_db[p, 25].mat != mat:  # b4
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 23].mat != mat:  # b2
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 6].mat != mat:  # s3
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side s3",p)
    return sumarea


def b4(voxel_db, p, dx, mat, sumarea):
    #    print("b4")
    if voxel_db[p - 1, 19].mat != mat:  # f2
        sumarea = sumarea + dx * dx * 0.25  # ; print("perp f2",sumarea)
    if voxel_db[p, 22].mat != mat:  # b1
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 24].mat != mat:  # b3
        sumarea = sumarea + 0.17677 * dx ** 2
    if voxel_db[p, 11].mat != mat:  # w2
        sumarea = sumarea + 0.3535 * dx ** 2  # ; print("diag side w2",p)
    return sumarea


def func(tag, voxel_db, p, nx, ny, nz, dx, mat):
    sumarea = 0.0
    if tag == 'n1':
        sumarea = n1(voxel_db, p, nx, ny, dx, mat, sumarea)
    elif tag == 'n2':
        sumarea = n2(voxel_db, p, nx, ny, dx, mat, sumarea)
    elif tag == 'n3':
        sumarea = n3(voxel_db, p, nx, ny, dx, mat, sumarea)
    elif tag == 'n4':
        sumarea = n4(voxel_db, p, nx, ny, dx, mat, sumarea)

    elif tag == 's1':
        sumarea = s1(voxel_db, p, nx, ny, dx, mat, sumarea)
    elif tag == 's2':
        sumarea = s2(voxel_db, p, nx, ny, dx, mat, sumarea)
    elif tag == 's3':
        sumarea = s3(voxel_db, p, nx, ny, dx, mat, sumarea)
    elif tag == 's4':
        sumarea = s4(voxel_db, p, nx, ny, dx, mat, sumarea)

    elif tag == 'w1':
        sumarea = w1(voxel_db, p, ny, dx, mat, sumarea)
    elif tag == 'w2':
        sumarea = w2(voxel_db, p, ny, dx, mat, sumarea)
    elif tag == 'w3':
        sumarea = w3(voxel_db, p, ny, dx, mat, sumarea)
    elif tag == 'w4':
        sumarea = w4(voxel_db, p, ny, dx, mat, sumarea)

    elif tag == 'e1':
        sumarea = e1(voxel_db, p, ny, dx, mat, sumarea)
    elif tag == 'e2':
        sumarea = e2(voxel_db, p, ny, dx, mat, sumarea)
    elif tag == 'e3':
        sumarea = e3(voxel_db, p, ny, dx, mat, sumarea)
    elif tag == 'e4':
        sumarea = e4(voxel_db, p, ny, dx, mat, sumarea)

    elif tag == 'f1':
        sumarea = f1(voxel_db, p, dx, mat, sumarea)
    elif tag == 'f2':
        sumarea = f2(voxel_db, p, dx, mat, sumarea)
    elif tag == 'f3':
        sumarea = f3(voxel_db, p, dx, mat, sumarea)
    elif tag == 'f4':
        sumarea = f4(voxel_db, p, dx, mat, sumarea)

    elif tag == 'b1':
        sumarea = b1(voxel_db, p, dx, mat, sumarea)
    elif tag == 'b2':
        sumarea = b2(voxel_db, p, dx, mat, sumarea)
    elif tag == 'b3':
        sumarea = b3(voxel_db, p, dx, mat, sumarea)
    elif tag == 'b4':
        sumarea = b4(voxel_db, p, dx, mat, sumarea)

    return sumarea
