#-*- coding: utf-8 -*-
#Function about PDB data and file

FLAG_TO_DEFECT = ["All", "Deep", "Shallow"]

# modify data from PDB with the good lipid name
# Rename lipid from 3 letters to 4 
def modifyPDBdata(datapdb, dico,startID,endID):
    for i, val in enumerate(datapdb):
        if val[0:4] == "ATOM":
            for key in dico:
                if int(val[startID:endID]) in dico[key]:
                    #print(key,len(val[17:21].strip()))
                    datapdb[i]=val[:17]+key+" "+val[21:]
    return datapdb

#determine position of res_identifier in pdb file
def determine_pos_resid(data):
    # The residue number can exceed PDB_limit(9999)
    flagID=0
    i=0
    while flagID == 0:
        if data[i][:4] == "ATOM":
            flagID=1
            try:
                resid=int(data[i][22:26])
            except:
                startID=21
                endID=27
            else:
                startID=22
                endID=26
        i=i+1
    #print(startID,endID)
    return startID,endID

#find num frame in pdb file (line MODEL) or 1 by default
def find_numframe(data):
    num_frame = 1
    for line in data:
        if line[:6].strip() == "MODEL":
            num_frame = int(line.split()[-1])
    return num_frame
    
# create output file of atoms in listres in PDB format
def outputPDB_leaflet(data, listres, outputname, startID, endID, num_frame=1):
    with open(outputname,"w") as f:
        f.write("MODEL      %3d\n"%(num_frame))
        for val in data:
            if val[0:4] == "ATOM":
                res_number=int(val[startID:endID])
                if res_number in listres:
                    f.write(val)
        f.write("ENDMDL\n")

# create output file of atm_name in listres in PDB format
def outputPDB_leafletAtm(data, listres, outputname, startID, endID, atm_name):
    with open(outputname,"w") as f:
        for val in data:
            if val[0:4] == "ATOM":
                res_number=int(val[startID:endID])
                atm=val[12:16].strip()
                if res_number in listres and atm == atm_name:
                    f.write(val)

# create output file of lipids in listres and molecules in listLD in PDB format
def outputPDB_leafletLD(data, listres, outputname, startID, endID, listLD):
    with open(outputname,"w") as f:
        for val in data:
            if val[0:4] == "ATOM":
                res_number=int(val[startID:endID])
                if res_number in listres:
                    f.write(val)
                if res_number in listLD:
                    f.write(val)

def write_a_pdb_line(file, nb_atm, atm_name, nb_res, coords, nb_defect):
    file.write("%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%
                ("ATOM  ",nb_atm,"  H1", atm_name, nb_res ,float(coords[0]),float(coords[1]),
                float(coords[2]),1.0,nb_defect))

# create output file with Total Matrix in pdb format
def outputPDB_Total_matrix(outputname, FlagPDtype, leaflet, num_frame, listX, 
                            listY, Rpos, Matrix_final):
    outputname = outputname + "_Total" + leaflet+ "_" + FLAG_TO_DEFECT[FlagPDtype] + ".pdb"
    with open(outputname,"w") as f:
        f.write("MODEL      %3d\n"%(num_frame))
        nb=0
        for j, sliceX in enumerate(listX):
            nb+=1
            for k, sliceY in enumerate(listY):
                coordtmp = [sliceX, sliceY, Rpos]
                if Matrix_final[j][k] == "NA":
                    write_a_pdb_line(f, nb, "EDG",  nb, coordtmp, -1)
                else:
                    write_a_pdb_line(f, nb, "MAT", nb, coordtmp, Matrix_final[j][k])
        f.write("ENDMDL\n")


# create output file with defects only in pdb format
def outputPDB_defects(outputname, FlagPDtype, leaflet, num_frame, listX, listY, Rpos, 
                        Matrix_fin, cluster_edge):
    outputname = outputname + "_Defect" + leaflet+ "_"+FLAG_TO_DEFECT[FlagPDtype]+".pdb"
    with open(outputname,"w") as f:
        f.write("MODEL      %3d\n"%(num_frame))
        nb=0
        dicoDef={}
        for j, sliceX in enumerate(listX):
            for k, sliceY in enumerate(listY):
                num_def = Matrix_fin[j][k]
                if num_def !=0. and num_def not in cluster_edge:
                    if int(num_def) not in dicoDef:
                        dicoDef[int(num_def)]=[]
                    nb+=1
                    coordtmp=[sliceX,sliceY, Rpos]
                    dicoDef[int(num_def)].append(
                                "%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%
                                ("ATOM  ", nb,"  H1", "DEF", int(num_def),
                                float(coordtmp[0]),float(coordtmp[1]),float(coordtmp[2]),
                                1.0,num_def))
                        #f.write("%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%("ATOM  ", nb,"  H1", "DEF", int(Matrix_fin[j][k]),float(coordtmp[0]),float(coordtmp[1]),float(coordtmp[2]),1.0,Matrix_fin[j][k]))
                    #else:
                    #   nb+=1
                    #  	coordtmp=[sliceX,sliceY, Rpos]
                    # 	dicoDef[int(Matrix_fin[j][k])].append(write_a_pdb_line(nb, 
                    #	"%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%("ATOM  ", nb,"  H1", "DEF", int(Matrix_fin[j][k]),float(coordtmp[0]),float(coordtmp[1]),float(coordtmp[2]),1.0,Matrix_fin[j][k]))
        #print(dicoDef)
        listdef = list(dicoDef.keys())
        listdef.sort()
        nb=0
        for l in listdef:
            for line in dicoDef[l]:
                nb+=1
                f.write("%s%5d%s"%(line[:6],nb,line[11:]))
        f.write("ENDMDL\n")

# create output file with defects only in pdb format
def outputPDB_defects_old(outputname, FlagPDtype, leaflet, num_frame, listX, listY, Rpos, 
                            Matrix_fin, cluster_edge):
    outputname = outputname + "_Defect" + leaflet+ "_"+FLAG_TO_DEFECT[FlagPDtype]+".pdb"
    with open(outputname,"w") as f:
        f.write("MODEL      %3d\n"%(num_frame))
        nb=0
        for j, sliceX in enumerate(listX):
            for k, sliceY in enumerate(listY):
                num_def = Matrix_fin[j][k]
                if num_def !=0. and num_def not in cluster_edge:
                    nb+=1
                    coordtmp=[sliceX, sliceY, Rpos]
                    write_a_pdb_line(f, nb, "DEF", int(num_def), coordtmp, num_def)
                    #f.write("%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%("ATOM  ", nb,"  H1", "DEF", int(Matrix_fin[j][k]),float(coordtmp[0]),float(coordtmp[1]),float(coordtmp[2]),1.0,Matrix_fin[j][k]))
        f.write("ENDMDL\n")

#create output file for packing defects in TXT format
# #defect size Xposition   Yposition
# 1    7     7.00    55.00 
# 2    1     7.00    58.00 
def outputTXT_defects(outputname, FlagPDtype, leaflet, dico_def_area, dico_def_coor, 
                        total_size, total_edge, listX, listY):
    outputname = outputname + "_" + leaflet+ "_"+FLAG_TO_DEFECT[FlagPDtype]+"_result.txt"
    with open(outputname,"w") as f:
        f.write("## MatrixSize %5d %5d \n"%(total_size-total_edge,total_size))
        f.write("## Total %4d %5d %4.2f %5.3f\n"%(len(dico_def_area),sum(dico_def_area.values()), sum(dico_def_area.values())/float(len(dico_def_area)), (sum(dico_def_area.values())*100.)/(total_size-total_edge)))
        i=0
        for key in dico_def_coor:
            # version with matrice position
            #f.write("%3d %4d   %6.2f   %6.2f \n"%(i+1, dico_def_area[key], dico_def_coor[key][0], dico_def_coor[key][1]))
            # version with x, y coordinates
            f.write("%3d %4d   %6.2f   %6.2f \n"%(i+1, dico_def_area[key], listX[dico_def_coor[key][0]], listY[dico_def_coor[key][1]]))
            i+=1

