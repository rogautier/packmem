#-*- coding: utf-8 -*-
# Functions to create, analyse list variables
# R. Gautier 2015

# create list[v1 to v2, step]
def create_list_ascend(v1, v2, step):
    liste = []
    i = v1
    while i <= v2:
        liste.append(i)
        i = i+step
    return liste

# create list[
def create_list_descend(v1, v2, step):
    liste = []
    i = v1
    while i >= v2:
        liste.append(i)
        i = i+step
    return liste

#from list[list] determine max length
def determine_lenMax(listeTot):
    return max([len(data) for data in listeTot])

# return the maximum value from list[list]
def max_value_list(data):
    return max(data)

# return the minimum value from list[list]
def min_value_list(data):
    return min(data)

# return min, max, mean from list[values]
def min_max(data):
    miniz = min(data)
    maxiz = max(data)
    meanz=(maxiz+miniz)/2.
    return miniz, maxiz, meanz
