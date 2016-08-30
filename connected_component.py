# -*- coding: utf-8 -*-

##### Here is implemented the Connected-component_labeling algorithm
##### Source: http://en.wikipedia.org/wiki/Connected-component_labeling
#####
##### P. Fuchs 08/2014
#####

import copy

def get_random_matrix(nrows,ncols):
    import random
    M=[]
    for i in range(nrows):
        line=[]
        for j in range(ncols):
            line.append(random.randint(0,1))
        M.append(line)
    return M

def get_random_matrix2(nrows,ncols):
    import random
    M=[]
    for i in range(nrows):
        line=[]
        for j in range(ncols):
            line.append(random.randint(0,1))
        M.append(line)
    for i in range(nrows):
        M[i][ncols/2]=0
    for j in range(ncols):
        M[nrows/2][j]=0
    return M

    
def read_matrix(filename): # WARNING: no extra empty line at the end
    M=[] ; nrows = 0 ; ncols = 0
    with open(filename) as f:
        for l in f.readlines():
            M.append(l.split())
    nrows = len(M)
    ncols = len(M[0])
    for l in M:
        if len(l) != ncols:
            exit("Pb with matrix, not the same nb of cols!!!")
    for i in range(nrows):
        for j in range(ncols):
            M[i][j] = int(M[i][j])
    return M,nrows,ncols

# each 0 -> 1; each 1 -> 0
def rev_matrix(M,nrows,ncols):
    for i in range(nrows):
        for j in range(ncols):
            if M[i][j] == 0:
                M[i][j] = 1
            elif M[i][j] == 1:
                M[i][j] = 0
            else:
                exit("Pb in rev_matrix, not a binary matrix")
    return M

def get_Romain_matrix():
    M=[] ; nrows = 0 ; ncols = 0
    with open("testmatrix_romain.txt","r") as f:
        for l in f.readlines():
            new_line=[]
            l = l.replace("[","").replace("]","").replace("'","").split(",")
            for elt in l:
                if elt.strip() == 'NA':
                    new_line.append(1)
                else:
                    value = float(elt)
                    if value < 0.000001: #case 0.0
                        new_line.append(1)
                    else:
                        new_line.append(0)
            M.append(new_line)
            nrows += 1
    ncols = len(M[0])
    for l in M:
        if len(l) != ncols:
            exit("Pb with matrix, not the same nb of cols!!!")
    return M,nrows,ncols
    
def get_binary_matrix(filename):
    # read data: M is the initial matrix with 1 & 0
    M=[]
    with open(filename,'r') as f :
        for li in f:
            M.append(list(li[:-1]))
    nrows = len(M)
    ncols = len(M[0])
    # check if ncols consistent
    for i in range(nrows):
        if len(M[i]) != ncols:
            exit("Matrix has not the same number of columns in each row")
    # convert M to int
    for i in range(nrows):
        for j in range(ncols):
            M[i][j] = int(M[i][j])
    return M,nrows,ncols

def print_matrix(M,nrows,ncols):
    s = ""
    for i in range(nrows):
        for j in range(ncols):
            s += "%3i " % M[i][j]
        s += "\n"
    return s

def print_matrix2file(filename,M,nrows,ncols):
    with open(filename,"w") as f:
        f.write(print_matrix(M,nrows,ncols))

def init_matrix_0(nrows,ncols):
    m=[]
    for i in range(nrows):
        line=[]
        for j in range(ncols):
            line.append(0)
        m.append(line)
    return m

def get_uniq_labels(nb_labels,list_equiv_labels):
    # list_equiv_labels looks like that: [[1,2], [3,4,5,8], [6,7]...]
    # At the end we want to get something like this:
    # {1: 1, 2: 1, 3: 3, 4: 3, 5: 3, 6: 6, 7: 6, 8:3}
    dico_uniq_labels = {} ; root_labels = []
    for label in range(1,nb_labels+1):
        # is label in list_equiv_labels?
        is_label_in_list_equiv_labels = 0
        for tmp_list_labels in list_equiv_labels:
            if label in tmp_list_labels:
                # assign the root label to the current label
                dico_uniq_labels[label] = tmp_list_labels[0]
                is_label_in_list_equiv_labels = 1
                # is it a root label (thus at the first position of the list)?
                if label == tmp_list_labels[0]:
                    if label not in root_labels: # in case of big cluster, it avoids writing multiple times the same label!
                        root_labels.append(label)
        if not is_label_in_list_equiv_labels:
            # if we end up overhere, the label is a root label but unique
            dico_uniq_labels[label] = label
            root_labels.append(label)
    DEBUG=0
    if DEBUG:
        print("list_equiv_labels:" , list_equiv_labels)
        print("nb labels=%i" % nb_labels)
        print("Dictionnary of equivalence:")
        print(dico_uniq_labels)
        print("nb of root labels=%i" % len(root_labels))
        print("root labels are:"); print(root_labels)
        #exit()
    return dico_uniq_labels, root_labels

# returns the intersection of 2 lists
def intersect(l1, l2):
    return list(set(l1) & set(l2))

def merge_avoid_duplicate_sort_list(l1,l2):
    new_l = list(set(l1 + l2))
    new_l.sort()
    return new_l

def is_duplicate_in_list_of_lists(list_of_lists, starting_index):
    L = len(list_of_lists)
    for i in range(starting_index+1,L):
        if len(intersect(list_of_lists[starting_index],list_of_lists[i])) > 0:
            return True
    return False

def is_empty_sublist(list_of_lists):
    L = len(list_of_lists)
    for i in range(L):
        if len(list_of_lists[i]) == 0:
            return True
    return False
    
# merge a list of list
def merge_list_of_lists(list_of_lists):
    # we want to merge the sublists that intersect in a list_of_lists, e.g.
    # init list_of_lists:  [[1, 2, 3, 9], [4, 5, 8], [1, 7, 9], [8, 9, 10], [13, 16], [12, 44], [16, 54]]
    # final list_of_lists: [[1, 2, 3, 4, 5, 7, 8, 9, 10], [13, 16, 54], [12, 44]]
    #print "Intial list_of_lists", list_of_lists
    L = len(list_of_lists)
    starting_index = 0
    while starting_index <= L - 1:
        # we try to merge sublist starting_index with all other sublists
        while is_duplicate_in_list_of_lists(list_of_lists,starting_index):
            for i in range(starting_index+1,L):
                if len(intersect(list_of_lists[starting_index],list_of_lists[i])) > 0:
                    # OK we merge sublist starting_index with sublist i
                    list_of_lists[starting_index] = merge_avoid_duplicate_sort_list(list_of_lists[starting_index], list_of_lists[i])
                    # we delete sublist i
                    list_of_lists[i] = []
        # we increment starting index
        starting_index += 1
    # cleanup empty sublists
    while is_empty_sublist(list_of_lists):
        L = len(list_of_lists)
        for i in range(L):
            if len(list_of_lists[i]) == 0:
                list_of_lists.pop(i)
                break
    #print "Final list_of_lists" , list_of_lists
    return list_of_lists


# this is the important fct 
def get_connected_components(M,nrows,ncols,val_bin=0):
##    print_matrix2file("matrix_init.txt",M,nrows,ncols)
    # Initialize a new matrix (for labels)
    M_labels = init_matrix_0(nrows,ncols) #set each cell to 0 (integer)
    # counter for the number of labels
    nb_labels = 0 
    # list for storing equivalent labels (see below)
    list_equiv_labels = []

    # First pass
    for i in range(nrows):
        for j in range(ncols):
            if M[i][j] == val_bin: # this is a positive cell
                # create a list that contains the positive neighbors
                # e.g.: [[1, 2]] means the neighboring cell [1,2] (of the current cell) is positive
                # e.g.: [[6, 6], [6, 4]] means the neighboring cells (of the current cell)
                #       [6,6] and [6,4] are positive
                positive_neighbors = []
                # look for connectivity of the current positive cell
                if i > 0:
                    if j > 0 and M[i-1][j-1] == val_bin: # North-West
                        positive_neighbors.append([i-1,j-1])                
                    if M[i-1][j] == val_bin: # North
                        positive_neighbors.append([i-1,j])
                    if j < ncols - 1 and M[i-1][j+1] == val_bin: # North-East
                        positive_neighbors.append([i-1,j+1])
                if j > 0 and M[i][j-1] == val_bin: # West
                    positive_neighbors.append([i,j-1])
                # if no neighbor is positive, this cell is a new label
                if len(positive_neighbors) == 0:
                    nb_labels += 1
                    M_labels[i][j] = nb_labels
                # if only 1 neighbor is positive, assign the neighbor label to current cell
                if len(positive_neighbors) == 1:
                    row_uniq_neighbor = positive_neighbors[0][0]
                    col_uniq_neighbor = positive_neighbors[0][1]
                    M_labels[i][j] = M_labels[row_uniq_neighbor][col_uniq_neighbor]
                # if multiple neighbors are positive, assign one of their label
                #  --> it doesn't matter which one: arbitrarily, we choose the first one
                # in the list positive_neighbors
                if len(positive_neighbors) > 1:
                    #print("positive_neighbors:", positive_neighbors)
                    # loop over positive_neighbors and store their label in a list called label_neighbors
                    label_neighbors = []
                    for neighbor in positive_neighbors: # neighbor is a list (e.g. [11,81]) containing coor of the neighbor in M
                        row_neighbor = neighbor[0]
                        col_neighbor = neighbor[1]
                        label_neighbors.append(M_labels[row_neighbor][col_neighbor])
                    # eliminate duplicate & sort
                    label_neighbors = list(set(label_neighbors))
                    label_neighbors.sort()
                    # assign the label (min of label_neighbors) to that cell
                    # I tested, any label label_neighbors work as well!
                    M_labels[i][j] = min(label_neighbors)
                    
                    # in case of multiple labels, keep track they are equivalent in list_equiv_labels:
                    # e.g. [[1,2],[3,4,5,7],[6,8], ...] (each sublist is ordered)
                    #  -> means 1&2 are equiv, 3,..,7 are equiv, etc
                    if len(label_neighbors) > 1:
                        #print("positive_neighbors:", positive_neighbors)
                        #print("cell %3i-%3i: label %3i, label_neighbors: %s" %(i,j,M_labels[i][j],label_neighbors))
                        #print()
                        # is list_equiv_labels empty? -> easy to fill in
                        if len(list_equiv_labels) == 0:
                            list_equiv_labels.append(label_neighbors)
                        # otherwise, it's not empty! -> and we have to fill it in
                        else:
                            # is one of the label_neighbors present in list_equiv_labels?
                            is_present_in_list_equiv_labels = 0
                            for k,tmplist in enumerate(list_equiv_labels):
                                if len(intersect(tmplist,label_neighbors)) > 0:
                                    # found it! So we fuse the lists, avoid duplicates and sort!
                                    list_equiv_labels[k] = merge_avoid_duplicate_sort_list(list_equiv_labels[k] , label_neighbors)
                                    is_present_in_list_equiv_labels = 1
                            # label is not present, so we create a new sublist
                            if not is_present_in_list_equiv_labels:
                                list_equiv_labels.append(label_neighbors)

    # There are redundancies in list_equiv_labels -> clean them up
    list_equiv_labels = merge_list_of_lists(list_equiv_labels)
    # uncomment the following for an output
##    print("First pass done!")
##    print_matrix2file("matrix_1after_1stpass.txt",M_labels,nrows,ncols)
##    print("Intermediary matrix of labels printed to file")

    # Second pass
    # Build a dictionnary for assigning a uniq label
    # should look like this: {1: 1, 2: 1, 3: 3, 4: 4}
    # --> label 1 & 2 are equiv, thus we set label 2 to label 1
    dico_of_uniq_labels, root_labels = get_uniq_labels(nb_labels,list_equiv_labels)
    area_clusters = {} # to store the area of each clusteer
    coor_clusters = {} # to store the coor (in the matrix) of the 1st element of each cluster
    for root_label in root_labels:
        area_clusters[root_label] = 0
    for i in range(nrows):
        for j in range(ncols):
            if M_labels[i][j]:
                M_labels[i][j] = dico_of_uniq_labels[M_labels[i][j]]
                area_clusters[M_labels[i][j]] += 1
                if M_labels[i][j] not in coor_clusters:
                    coor_clusters[M_labels[i][j]]=[i,j]
                
    # uncomment the following for an output                
##    print("Second pass done!")
##    print_matrix2file("matrix_2after_2ndpass.txt",M_labels,nrows,ncols)
##    print("Final matrix of labels printed to file")
    return M_labels, root_labels, area_clusters, coor_clusters

def get_clusters_on_the_edge(M_labels,nrows,ncols):
    clusters_on_the_edge = []
    # check clusters on 1st & last row
    #print(M_labels[0])
    #print(M_labels[nrows-1])
    clusters_on_the_edge += list(set([i for i in M_labels[0] if i != 0]))
    clusters_on_the_edge += list(set([i for i in M_labels[nrows-1] if i != 0]))
    # check clusters on 1st & last col
    for i in range(nrows):
        if M_labels[i][0] not in clusters_on_the_edge and M_labels[i][0] != 0:
            clusters_on_the_edge.append(M_labels[i][0])
        if M_labels[i][ncols-1] not in clusters_on_the_edge and M_labels[i][ncols-1] != 0:
            clusters_on_the_edge.append(M_labels[i][ncols-1])
    # get rid of duplicates & sort
    clusters_on_the_edge = sorted(list(set(clusters_on_the_edge)))
    return clusters_on_the_edge

def delete_NApoints_inside(clustPb, Matrix_labels, root_labels, area_clusters, nrows, ncols):
    #### supprimer les defauts totalament NA
    clust_2_delete=[]
    for key in clustPb:
        if area_clusters[key] == clustPb[key]:
            clust_2_delete.append(key)
    if len(clust_2_delete) != 0:
        for num in clust_2_delete:
            tmp=root_labels.index(num)
            del root_labels[tmp]

    area_clusters = {} # to store the area of each cluster
    coor_clusters = {} # to store the coor (in the matrix) of the 1st element of each cluster
    for root_label in root_labels:
        area_clusters[root_label] = 0
    for i in range(nrows):
        for j in range(ncols):
            if Matrix_labels[i][j]:
                area_clusters[Matrix_labels[i][j]] += 1
                if Matrix_labels[i][j] not in coor_clusters:
                    coor_clusters[Matrix_labels[i][j]]=[i,j]
    return root_labels, area_clusters, coor_clusters
    
