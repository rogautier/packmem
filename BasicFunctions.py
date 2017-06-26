#-*- coding: utf-8 -*-
# R. Gautier 2014
# Basic functions for packing defects determination

# load file return list[line(string)]
def read_file(filename):
    try: 
        with open(filename) as f:
            data = f.readlines()
    except : 
        print("ERROR : Something went wrong with the file %s" % filename)
    return data

# open VdW radius -> dict["RES(3L) ATOMname"]=radius in A 
def read_radius(filename):
    dico={}
    data = read_file(filename)
    for i in range(0, len(data)):
        data[i]=data[i].strip()
        if data[i] != "":
            data[i]=data[i].split()
            cle=data[i][0]+" "+data[i][1]
            dico[cle]=float(data[i][2])
    return dico

#new 2016-08-18
# open aliphatic table-> dict["RES(3L) ATOMname"]= p (polar) a (aliphatic) 
def read_aliphatic(filename):
    dico={}
    data = read_file(filename)
    for i in range(0, len(data)):
        data[i]=data[i].strip()
        if data[i] != "":
            data[i]=data[i].split()
            cle=data[i][0]+" "+data[i][1]
            dico[cle]=data[i][3]
    return dico

# Euclidian distance without square root
def dist(coord1, coord2):
    return ((coord1[0]-coord2[0])**2 +
           (coord1[1]-coord2[1])**2 +
           (coord1[2]-coord2[2])**2)

# Euclidian distance without square root
def dist2D(coord1, coord2):
    return ((coord2[0]-coord1[0])**2 +
           (coord2[1]-coord1[1])**2)

# distance for one axis without square root
def dist_oneAxis(axe1,axe2):
    return (axe1-axe2)**2

