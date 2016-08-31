#-*- coding: utf-8 -*-
# Funtions for matrix
import math
import BasicFunctions as bfrg

# matrix size  (square size = 1A)
SIZE=1.0

# calculate limit size (half diagonal) 3D
SIZE_SIDE=0.5*(math.sqrt((math.pow(SIZE, 2))* 3.))

# get radius from [Lipid_name Atom_name]
def get_radius(radius, res_name, atom_name):
    return radius[res_name+" "+atom_name]

# get aliphatic flag [Lipid_name Atom_name]
def get_aliphatic(aliphatic, res_name, atom_name):
    return aliphatic[res_name+" "+atom_name]

# fill the matrix for each atom
def fill_matrix(matrix,coordtmp,res_number,res_name,atom_name, listX,listY,listZ,
                radius_res, FlagPDtype, aliph_atoms):
    v=5
    # for each lipid get listZ
    xLip = coordtmp[0]
    yLip = coordtmp[1]
    zLip = coordtmp[2]
    # works only around X, Y position of the matrix +/- v cases
    iX,iY = find_X_Y(coordtmp, listX, listY)
    iX = check_edges(iX, v, len(listX))
    iY = check_edges(iY, v, len(listY))
    listXM = listX[iX - v : iX + (v + 1)]
    listYM = listY[iY - v : iY + (v + 1)]

    dist_lim = (SIZE+radius_res) ** 2
    dist_meet = (SIZE_SIDE + radius_res) ** 2

    for indZ, sliceZ in enumerate(listZ):
        distZ_to_sliceZ = bfrg.dist_oneAxis(zLip,sliceZ)
        if distZ_to_sliceZ > dist_lim:
            continue
        for indX,sliceX in enumerate(listXM):
            distX_to_sliceX = bfrg.dist_oneAxis(xLip, sliceX)
            if distX_to_sliceX > dist_lim:
                continue
            for indY,sliceY in enumerate(listYM):
                distY_to_sliceY = bfrg.dist_oneAxis(yLip, sliceY)
                if distY_to_sliceY > dist_lim:
                    continue
                X = indX + (iX - v)
                Y = indY + (iY - v)
                Z = indZ
                coordCenter = [float(sliceX), float(sliceY), float(sliceZ)]
                distance = bfrg.dist(coordCenter, coordtmp)
                if matrix[X][Y]=="NA":
                    matrix[X][Y]=0.
                val_mat = matrix[X][Y]
                if distance <= dist_meet:
                    # shallow PDefects (2) or all PDefects (0)
                    if FlagPDtype is 2 or FlagPDtype is 0:
                        matrix[X][Y] = setDefects(atom_name, aliph_atoms, val_mat)
                    # deep PDefects (1)
                    else:
                        matrix[X][Y] +=1.
    return matrix


# fill matrix cell depending on the defect type (Deep, Shallow)
def setDefects(atom_name, list_aliphatic, val_mat):
    if list_aliphatic == "a":
        val_mat += 0.001
    else:
        val_mat += 1.
    return val_mat

def initialize_matrix3D(val1, val2, val3, defaut):
    Matrix=list(range(0,val1))
    for i in range(0, val1):
        Matrix[i]=list(range(0,val2))
        for j in range(0, val2):
            Matrix[i][j]=list(range(0,val3))
            for k in range(0, val3):
                Matrix[i][j][k]=defaut
    return Matrix

def initialize_matrix2D(val1, val2, defaut):
    Matrix=list(range(0,val1))
    for i in range(0, val1):
        Matrix[i]=list(range(0,val2))
        for j in range(0, val2):
            Matrix[i][j]=defaut
    return Matrix       

def binarize_matrix(matrix, matrix_ini, val=0.):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[0])):
            #print(matrix[i][j])
            if matrix[i][j] == "NA":
                matrix_ini[i][j]=0.
            else:
                if matrix[i][j] > val:
                    matrix_ini[i][j]=1.
                else:
                    matrix_ini[i][j]=0.
    return matrix_ini

def binarize_matrix_whithout0(matrix, matrix_ini, val1=0, val2=0.99):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[0])):
            #print(matrix[i][j])
            if matrix[i][j] == "NA":
                matrix_ini[i][j]=0.
            else:
                if (matrix[i][j] > val1 and  matrix[i][j] < val2):
                    matrix_ini[i][j]=0.
                else:
                    matrix_ini[i][j]=1.
    return matrix_ini

def binarize_matrix_whithout0_bis(matrix, matrix_ini, val=0.99):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[0])):
            #print(matrix[i][j])
            if matrix[i][j] == "NA":
                matrix_ini[i][j]=0.
            else:
                if (matrix[i][j] > 0. and  matrix[i][j] < val):
                    matrix_ini[i][j]=0.
                else:
                    matrix_ini[i][j]=1.
    return matrix_ini

def clean_NA_inside(Matrix_labels, clust_edge, Matrix_ini, total_edge):
    clustPb={}
    for i in range(0, len(Matrix_ini)):
        for j in range(0, len(Matrix_ini[0])):
            if Matrix_ini[i][j] == "NA" and Matrix_labels[i][j] not in clust_edge:
                if Matrix_labels[i][j] in clustPb:
                    clustPb[Matrix_labels[i][j]]+=1
                else:
                    clustPb[Matrix_labels[i][j]]=1
                Matrix_labels[i][j]=clust_edge[0]
                total_edge+=1
    return Matrix_labels, total_edge, clustPb
            
def determinelastNA(matrix):
    for i in range(len(matrix)-1, -1, -1):
        print(i, matrix[i][0][0])

# function that determine from x,y coordinates the central value in matrix
def find_X_Y(coord, matX, matY):
    tmpX = int(coord[0])
    tmpY = int(coord[1])
    iX = matX.index(tmpX)
    iY = matY.index(tmpY)
    #print(tmpX, matX, iX)
    return iX,iY

# compute the difference between the last value in matrixZ and the Z for one atom
def find_Z(coordZ, matZ):
    return matZ[-1]-coordZ

def check_edges(val, val_lim1, val_lim2):
    if val < val_lim1:
        return val_lim1
    elif val >= val_lim2-val_lim1:
        return val_lim2-(val_lim1+1)
    else:
        return val

def find_indexfromvalue(value, liste):
    mini=100.
    mini_i=len(liste)
    for i, val in enumerate(liste):
        if math.fabs(val-value) < mini:
            mini=math.fabs(val-value)
            mini_i=i
    return mini_i
    
def modify_matrix(mat1, mat2, listval):
    for i in range(0, len(mat1)):
        for j in range(0, len(mat1[0])):
            if mat1[i][j] in listval:
                mat2[i][j]= 0.0
    return mat2
    
