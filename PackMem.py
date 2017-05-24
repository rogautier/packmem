#!/usr/bin/python
#-*- coding: utf-8 -*-
# Pg compute matrix z and x/y for determine packing defects
# R. Gautier april 2015
# version matrix3D
# modif parameters 2016-08-18 (new aliphatic dictionnary)

import sys
import math
import argparse

import listes as l
import matrix as m
import pdb as pdb
import connected_component as cc
import dico as d
import BasicFunctions as bfrg
import param as p


##########################################################################################
##### main
if __name__ == '__main__':
    
    """
    Python
    Script to compute Packing defect in flat bilayers, Use pdb file (from for example editconf) 
    be careful to PBC to create pdb file
    Lipids parameters adapted to Berger lipid FF (with corrections) Be careful to the atoms 
    name if you use other lipids
    fileRadius.txt example:
    DOPC  C02 1.875 a (aliphatic)
    DOPC  O8  1.48 n (non aliphatic)
    DOPC  C25 1.98 n (non aliphatic)
    Output files:
    outputname_Up/Lo_Shallow/Deep/All_result.txt
    ## MatrixSize  9646  9801         # Membrane matrix size, Total matrix size
    ## Total   51   582 11.41 6.034   # number of packing defects, total area of packing defects, average size, pourcent of membrane (Membrane matrix size)
    1    1    77.00     2.00    # for each packing defect num size x_position (first pixel) y_mean (first pixel)
    2    1    21.00     3.00 
    3    8    59.12     7.00 

    outputname_TotalUp/Lo_Shallow/Deep/All.pdb # PDB format only if verbose
    Matrix x,y with, for each cell, the value of "packing defect" (in B_factor column, the last column) 
    if 0 = Deep packing defects The value increases with the number of atoms in the cell.
    if >0 and < 1 = Shallow defect
    if -1 = Edge
    ATOM      3   H1 EDG     3       4.000  54.000  56.940  1.00 -1.00
    ATOM      3   H1 MAT     3       4.000  54.000  56.940  1.00  0.00
    ATOM      3   H1 MAT     3       4.000  54.000  56.940  1.00  0.00

    outputname+"_DefectUp/Lo_Shallow/Deep/All.pdb # PDB format only if verbose
    Packing defects in pdb format:
    the residue number corresponds to the different packing defects.
    ATOM      5   H1 DEF     2       6.000  22.000  56.940  1.00  2.00
    ATOM      5   H1 DEF     3       6.000  85.000  56.940  1.00  3.00
    ATOM      5   H1 DEF     3       6.000  86.000  56.940  1.00  3.00
    ATOM      6   H1 DEF     3       7.000  86.000  56.940  1.00  3.00
    """
    #######PARAMETRES et INPUT#####
    
    try:
        outputname = "output"
        pdbout = 0
        parser = argparse.ArgumentParser(description = 'Arguments for the app')
        parser.add_argument('-i', action = 'store', dest = 'filename', help = 'PDB File')
        parser.add_argument('-r', action = 'store', dest = 'filesrad',
                            help = 'File for radius (default vdw_radiiFinal2014.txt)',
                            default = 'vdw_radiiFinal2014.txt')
        parser.add_argument('-o', action = 'store', dest = 'outputname',
                            default = 'output', help = 'Name for output file (default output)')
        parser.add_argument('-d', action = 'store', dest = 'dist_suppl_Z', default = 1.0, type = float,
                            help = '(default value 1.0)')
        parser.add_argument('-t', action = 'store', dest = 'pd_type', default = 'all',
                            help = 'Packing Defect Type (all/shallow/deep) default all')
        parser.add_argument('-p', action = 'store', dest = 'paramFile', default = 'param.txt',
                            help = 'File for lipid parameters')
        parser.add_argument('-n', action = 'store', dest = 'indexFile',
                            help = 'File for index (Gromacs ndx style). Only Lower/Upper group accepted')
        parser.add_argument('-v', dest = 'pdbout', action = 'store_true',
                            help = 'Increase the verbosity')
                            
        args = parser.parse_args()
        
        if len(sys.argv) == 1:
            parser.print_usage()
            sys.exit()
        
        DICT_3L, RESNAME_GLYC, LIPID = p.set_params(args.paramFile)

        # determine packing defedct type (all/shallow/deep)
        if args.pd_type == 'all':
            FlagPDtype = 0
        elif args.pd_type == 'shallow':
            FlagPDtype = 2
        elif args.pd_type == 'deep' :
            FlagPDtype = 1
        else: 
            print('ERROR : Packing defect type not known. Please choose : all/deep/shallow')
            sys.exit()
            
        if args.dist_suppl_Z < 0.0 :
            print('ERROR : The distance for the -d option must be > 0')
            sys.exit()
            
        pdblines = bfrg.read_file(args.filename)
        radius = bfrg.read_radius(args.filesrad)
        aliphatic = bfrg.read_aliphatic(args.filesrad)
        

    except:
        print('Command line: PacMem.py -i file.pdb -r fileRadius.txt \
               -p param.txt -o output -d distGlyc -t deep/all/shallow [-n index file] [-v] or -h for help')
        sys.exit()

    #extract num_frame in pdb file (line MODEL) or 1 by default
    num_frame=pdb.find_numframe(pdblines)

    # determine the position of res_identifier ([22:26] or [21:27]) depending the number of residues
    startID,endID = pdb.determine_pos_resid(pdblines)

    ###################Which lipids in bilayer?###################
    zaxis_byatoms = []
    xaxis_byatoms = []
    yaxis_byatoms = []
    dicoMb = {}
    #dict[lipid type]=list of residue number
    for line_atom in pdblines:
        if line_atom[0:4] == "ATOM":
            res_name=line_atom[17:21].strip()
            atm_name=line_atom[12:16].strip()
            zaxis_byatoms.append(float(line_atom[46:54]))
            xaxis_byatoms.append(float(line_atom[30:38]))
            yaxis_byatoms.append(float(line_atom[38:46]))
            key = res_name+  "_" + atm_name
            if res_name in LIPID:
                if key in DICT_3L:
                    if DICT_3L[key] in dicoMb:
                        dicoMb[DICT_3L[key]].append(int(line_atom[startID:endID]))
                    else:
                        dicoMb[DICT_3L[key]] = []
                        dicoMb[DICT_3L[key]].append(int(line_atom[startID:endID]))
            else :
                print("Lipid name %s not found in the parameter file submitted"%(res_name))
                sys.exit()
    xmin, xmax, xmean = l.min_max(xaxis_byatoms)
    ymin, ymax, ymean = l.min_max(yaxis_byatoms)
    zmin, zmax, zmean = l.min_max(zaxis_byatoms)

    # print(dicoMb)
    # modify the lipid name in data (3L->new3L to take account problem POPvsPOPS/POPE)
    pdblines = pdb.modifyPDBdata(pdblines, dicoMb, startID, endID)

    ########## extract UPPER/LOWER leaflet ##########
    if args.indexFile == None:
        lower_leaflet = []
        lower_listZ = {}
        upper_leaflet = []
        upper_listZ = {}
        for atm_line in pdblines :
            if atm_line[0:4] == "ATOM":
                # for the upper leaflet
                atom_name = atm_line[12:16].strip()
                res_name = atm_line[17:21].strip()
                z_coord = float(atm_line[46:54])
                if (atom_name == RESNAME_GLYC[res_name] and z_coord > zmean):
                    res_number = int(atm_line[startID:endID])
                    # append res_id in the list of lipids in upper leaflet
                    upper_leaflet.append(res_number)
                    # build list of Z-level for the Upper leaflet
                    tmp = l.create_list_ascend(z_coord - args.dist_suppl_Z, 
                                                    zmax + 1.0, m.SIZE)
                    tmp.reverse()
                    upper_listZ[res_number] = tmp
                # for the lower leaflet
                if (atom_name == RESNAME_GLYC[res_name] and z_coord < zmean):
                    res_number = int(atm_line[startID:endID])
                    # append res_id in the list of lipids in lower leaflet
                    lower_leaflet.append(res_number)
                    # build list of Z-level for the lower leaflet
                    tmp = l.create_list_descend(z_coord + args.dist_suppl_Z, 
                                                    zmin - 1.0, m.SIZE * -1)
                    tmp.reverse()
                    lower_listZ[res_number] = tmp
    else:
        (lower_leaflet, upper_leaflet) = p.read_ndx(args.indexFile)
        lower_listZ = {}
        upper_listZ = {}
        for atm_line in pdblines :
            if atm_line[0:4] == "ATOM":
                atom_name = atm_line[12:16].strip()
                res_name = atm_line[17:21].strip()
                z_coord = float(atm_line[46:54])
                res_number = int(atm_line[startID:endID])
                if (atom_name == RESNAME_GLYC[res_name] and 
                        res_number in upper_leaflet) :
                    tmp = l.create_list_ascend(z_coord - args.dist_suppl_Z, 
                                                    zmax + 1.0, m.SIZE)
                    tmp.reverse()
                    upper_listZ[res_number] = tmp
                if (atom_name == RESNAME_GLYC[res_name] and
                        res_number in lower_leaflet):
                    tmp = l.create_list_descend(z_coord + args.dist_suppl_Z, 
                                                    zmin - 1.0, m.SIZE * -1)
                    tmp.reverse()
                    lower_listZ[res_number] = tmp
                
                

    
    # create listX[xmin,xmax] and listY[ymin,ymax] 
    listX = l.create_list_ascend(int(xmin - 1.0), int(xmax + 1.0), m.SIZE)
    listY = l.create_list_ascend(int(ymin - 1.0), int(ymax + 1.0), m.SIZE)

    ##################################################  Initialize Matrix   ##############
    # initialize  matrix 2D for Upper and Lower leaflet
    MatrixUp = m.initialize_matrix2D(len(listX), len(listY), "NA")
    MatrixLo = m.initialize_matrix2D(len(listX), len(listY), "NA")
                
    #################################################  Compute Matrix    #################
    #value to limit the search around
    v = 5
    #for each atoms of lipids
    for atm_line in pdblines:
        if atm_line[0:4] == "ATOM":
            coordtmp = []
            atom_name = atm_line[12:16].strip()
            res_name = atm_line[17:21].strip()
            res_id = int(atm_line[startID:endID])
            radius_res=m.get_radius(radius, res_name, atom_name)
            aliph_atom=m.get_aliphatic(aliphatic, res_name, atom_name)
            coordtmp.append(float(atm_line[30:38])) #X
            coordtmp.append(float(atm_line[38:46])) #Y
            coordtmp.append(float(atm_line[46:54])) #Z
            # ##################################### Upper leaflet ########################
            if res_id in upper_leaflet :
                # limit the loop around dfZ
                dfZ = m.find_Z(coordtmp[2], upper_listZ[res_id])
                if dfZ < (1. * v):
                    MatrixUp= m.fill_matrix(MatrixUp, coordtmp, res_id, res_name,
                                            atom_name,listX, listY, upper_listZ[res_id], 
                                            radius_res,FlagPDtype, aliph_atom)

            # #################################### lower leaflet #########################
            if res_id in lower_leaflet :
                # limit the loop around dfZ
                dfZ = m.find_Z(coordtmp[2], lower_listZ[res_id])
                if dfZ > (-1. * v):
                    MatrixLo = m.fill_matrix(MatrixLo, coordtmp, res_id, res_name,
                                            atom_name, listX, listY, lower_listZ[res_id], 
                                            radius_res, FlagPDtype, aliph_atom)

    ########################preliminary process only for shallow defect and problem of the edges 
    ### to eliminate shallow defects on edges first: binarize on all defects and storage edges coord
    # shallow PDefects (2)
    if FlagPDtype == 2:
        # binarization Matrices according defect type [XM][YM]
        Matrix_UpbinM = m.initialize_matrix2D(len(listX),len(listY),0.)
        Matrix_LobinM = m.initialize_matrix2D(len(listX),len(listY),0.)
        Matrix_UpbinM = m.binarize_matrix_whithout0(MatrixUp, Matrix_UpbinM , -0.01, 0.99)
        Matrix_LobinM = m.binarize_matrix_whithout0(MatrixLo, Matrix_LobinM, -0.01, 0.99)
    
        # get temp packing defects 
        Matrix_labels_UpM, root_labels_UpM, area_clusters_UpM, coor_clusters_UpM = \
            cc.get_connected_components(Matrix_UpbinM,len(listX),len(listY))
        Matrix_labels_LoM, root_labels_LoM, area_clusters_LoM, coor_clusters_LoM = \
            cc.get_connected_components(Matrix_LobinM,len(listX),len(listY))

        # get cluster on the edge
        clust_edge_UpM=cc.get_clusters_on_the_edge(Matrix_labels_UpM,len(listX),len(listY))
        clust_edge_LoM=cc.get_clusters_on_the_edge(Matrix_labels_LoM,len(listX),len(listY))
        # print(clust_edge_UpM, clust_edge_LoM)

    ############################################## binarisation matrices ##################
    #print("PackingDefect generation")
    # matrix binarization
    Matrix_Upbin=m.initialize_matrix2D(len(listX),len(listY),0.)
    Matrix_Lobin=m.initialize_matrix2D(len(listX),len(listY),0.)
    # shallow PDefects (2)
    if FlagPDtype == 2:
        Matrix_Upbin = m.binarize_matrix_whithout0(MatrixUp, Matrix_Upbin , 0, 0.99)
        Matrix_Lobin = m.binarize_matrix_whithout0(MatrixLo, Matrix_Lobin, 0, 0.99)
    # all PDefects (2)
    elif FlagPDtype == 0:
        Matrix_Upbin = m.binarize_matrix_whithout0(MatrixUp, Matrix_Upbin , -0.01, 0.99)
        Matrix_Lobin = m.binarize_matrix_whithout0(MatrixLo, Matrix_Lobin, -0.01, 0.99)
    # deep
    else:
        Matrix_Upbin = m.binarize_matrix(MatrixUp, Matrix_Upbin, 0.)
        Matrix_Lobin = m.binarize_matrix(MatrixLo, Matrix_Lobin, 0.)

    # Exception when shallow defects
    if FlagPDtype == 2:
        #modify the binary matrix to take account edges (determined by all packing defects)
        #add the first coonected_component (edge
        Matrix_Upbin = m.modify_matrix(Matrix_labels_UpM, Matrix_Upbin, clust_edge_UpM)
        Matrix_Lobin = m.modify_matrix(Matrix_labels_LoM, Matrix_Lobin, clust_edge_LoM)

    ############################################## packing defects determination  ########
    # get packing defects 
    Matrix_labels_Up, root_labels_Up, area_clusters_Up, coor_clusters_Up = \
        cc.get_connected_components(Matrix_Upbin,len(listX),len(listY))
    Matrix_labels_Lo, root_labels_Lo, area_clusters_Lo, coor_clusters_Lo = \
        cc.get_connected_components(Matrix_Lobin,len(listX),len(listY))

    # get cluster on the edge
    clust_edge_Up = cc.get_clusters_on_the_edge(Matrix_labels_Up,len(listX),len(listY))
    clust_edge_Lo = cc.get_clusters_on_the_edge(Matrix_labels_Lo,len(listX),len(listY))

    # count area of the edge
    total_edge_Up = 0
    for key in clust_edge_Up:
        total_edge_Up += area_clusters_Up[key]
    total_edge_Lo = 0
    for key in clust_edge_Lo:
        total_edge_Lo += area_clusters_Lo[key]
    
    #clean dico defects (without edge)
    area_clusters_Up = d.del_key_dico(area_clusters_Up, clust_edge_Up)
    coor_clusters_Up = d.del_key_dico(coor_clusters_Up, clust_edge_Up)
    area_clusters_Lo = d.del_key_dico(area_clusters_Lo, clust_edge_Lo)
    coor_clusters_Lo = d.del_key_dico(coor_clusters_Lo, clust_edge_Lo)
   
    if FlagPDtype == 2:
        # Eliminate NA inside (deep not shallow defect)
        Matrix_labels_Lo, total_edge_Lo, clustPb_Lo = \
            m.clean_NA_inside(Matrix_labels_Lo, clust_edge_Lo, MatrixLo, total_edge_Lo)
        root_labels_Lo, area_clusters_Lo, coor_clusters_Lo = \
            cc.delete_NApoints_inside(clustPb_Lo, Matrix_labels_Lo, root_labels_Lo, 
                                       area_clusters_Lo, len(listX),len(listY))
        Matrix_labels_Up, total_edge_Up, clustPb_Up = \
            m.clean_NA_inside(Matrix_labels_Up, clust_edge_Up, MatrixUp, total_edge_Up)
        root_labels_Up, area_clusters_Up, coor_clusters_Up = \
            cc.delete_NApoints_inside(clustPb_Up, Matrix_labels_Up, root_labels_Up, 
                                       area_clusters_Up, len(listX),len(listY))
        
        # reclean dico defects (without edge)
        area_clusters_Up = d.del_key_dico(area_clusters_Up, clust_edge_Up)
        coor_clusters_Up = d.del_key_dico(coor_clusters_Up, clust_edge_Up)
        area_clusters_Lo = d.del_key_dico(area_clusters_Lo, clust_edge_Lo)
        coor_clusters_Lo = d.del_key_dico(coor_clusters_Lo, clust_edge_Lo)
    
    # ############################################# output text file  ####################
    total_size = len(listX) * len(listY)
    pdb.outputTXT_defects(args.outputname, FlagPDtype, "Up", area_clusters_Up, 
                          coor_clusters_Up, total_size, total_edge_Up, listX, listY)
    pdb.outputTXT_defects(args.outputname, FlagPDtype, "Lo", area_clusters_Lo, 
                          coor_clusters_Lo, total_size, total_edge_Lo, listX, listY)

    ############ output PDB files ########################################################
    # final matrix values PD (X,Y) with Z cooresponding to the max(Upper/lowerZlevel)
    if args.pdbout :
        valzmax=float(d.max_value_dico(upper_listZ))
        valzmin=float(d.min_value_dico(lower_listZ))
        pdb.outputPDB_leaflet(pdblines, upper_leaflet, args.outputname + "_Upper.pdb", 
                              startID, endID, num_frame)
        pdb.outputPDB_leaflet(pdblines, lower_leaflet, args.outputname + "_Lower.pdb", 
                              startID, endID, num_frame)
                        
        pdb.outputPDB_Total_matrix(args.outputname, FlagPDtype, "Up", num_frame,
                                   listX, listY, valzmax, MatrixUp)
        pdb.outputPDB_Total_matrix(args.outputname, FlagPDtype, "Lo", num_frame,
                                   listX, listY, valzmin, MatrixLo)
                                   
        pdb.outputPDB_defects(args.outputname, FlagPDtype, "Up", num_frame,
                              listX, listY, valzmax, Matrix_labels_Up, clust_edge_Up)
        pdb.outputPDB_defects(args.outputname, FlagPDtype, "Lo", num_frame,
                              listX, listY, valzmin, Matrix_labels_Lo, clust_edge_Lo)
