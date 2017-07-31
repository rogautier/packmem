#-*- coding: utf-8 -*-
# Functions about dictionaries
# R. Gautier 2015

# remove in dict keys from list_key
def del_key_dico(dico, list_key):
    if len(list_key) != 0:
        for m in list_key:
            del dico[m]
    return dico

# determine if duplicate in residue number
def detect_duplicate(dico):
    flagD=0
    dicotmp={}
    for key in dico:
        for val in dico[key]:
            if val in dicotmp:
                flagD=1
            else:
                dicotmp[val]=1
    return flagD

# return the maximum value from dico[key]=list
def max_value_dico(dico):
    return max([max(data2) for data2 in [dico[data] for data in dico.keys()]])

# return the minimum value from dico[key]=list
def min_value_dico(dico):
    return min([min(data2) for data2 in [dico[data] for data in dico.keys()]])

# from dico[liste] determine max length
def determine_lenMax2(dicoTot):
    return max([len(dicoTot[key]) for key in dicoTot])
