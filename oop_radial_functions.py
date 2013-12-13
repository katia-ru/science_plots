# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 04:07:29 2013

@author: katia
"""
from matplotlib import pyplot
import matplotlib.pyplot as plt
import numpy
import math
from scipy.interpolate import interp1d

N_COM = (180)
R_CUT = float(6)
R_BOND_oo = 3.0
R_BOND_oh = 2.3

pattern = "TIP3W"
#THRESHOLD = 5
#ширина окошка для радиальной функции
width = 0.1
width_angle = 1
N_COM = 1000

width_angle = 1

class WaterMolecule(object):
    def __init__(self, coords):
        self.o_x, self.o_y, self.o_z = coords[0]
        self.h1_x, self.h1_y, self.h1_z = coords[1]
        self.h2_x, self.h2_y, self.h2_z = coords[2]
       

def for_read(input_file):
    in_file = open(input_file, "r")
    oo = []
    h_1 = []
    h_2 = [] 
    for ii in in_file:
        word = ii.split()
        if pattern in ii:
            if "OH2" in ii:
                #пишем все координаты
                oo.append([float(q) for q in word[4:7]])
            elif "H1" in ii:
                h_1.append([float(q) for q in word[4:7]])
            elif "H2" in ii:
                h_2.append([float(q) for q in word[4:7]])
    return zip(oo, h_1, h_2)
    
def classFactory(input_file):
    all_coord = for_read(input_file)
    counter = 0
    wat_mols = []
    for ii in for_read(input_file):
        #wat_mols - массив объектов класса WaterMolecular
        wat_mols.append(WaterMolecule(ii))
    return wat_mols
        
def Dist_OO(wat_1, wat_2):
    square_dist = (wat_1.o_x - wat_2.o_x)**2 + (wat_1.o_y - wat_2.o_y)**2 \
    + (wat_1.o_z - wat_2.o_z)**2
    return sqrt(square_dist)
    
def Dist(coord_1, coord_2):
    square_dist = 0
    for i in xrange(3):
        square_dist += (coord_1[i] - coord_2[i])**2
    return sqrt(square_dist)
    

def HydrogeBondExist(wat_1, wat_2):
    dist_oo = Dist_OO(wat_1, wat_2)
    dist = [[0,0],[0,0]]
    marker = 0
    dist[0][0] = Dist([wat_1.h1_x, wat_1.h1_y, wat_1.h1_z], [wat_2.o_x, wat_2.o_y, wat_2.o_z])
    dist[0][1] = Dist([wat_1.h2_x, wat_1.h2_y, wat_1.h2_z], [wat_2.o_x, wat_2.o_y, wat_2.o_z])
    dist[1][0] = Dist([wat_2.h1_x, wat_2.h1_y, wat_2.h1_z], [wat_1.o_x, wat_1.o_y, wat_1.o_z])
    dist[1][1] = Dist([wat_2.h2_x, wat_2.h2_y, wat_2.h2_z], [wat_1.o_x, wat_1.o_y, wat_1.o_z])
    for i in xrange(2):
        for j in xrange(2):
            if dist[i][j] < R_BOND_oh:
                return 1
    return 0


def RadialDistrFuncAll(array_mol):
#    array_mol = []
    my_dict = {}
    x_axis, y_axis = [], []
 
    for ii in array_mol:
        for jj in array_mol[array_mol.index(ii):]:
            dist_oo = Dist_OO(ii, jj)
            dist_oo = int(dist_oo/width)*width
            if dist_oo < R_CUT:
                    if my_dict.has_key(dist_oo):
                        my_dict[dist_oo] += 1
                    else:
                        my_dict[dist_oo] = 1
    for ii in sorted(my_dict):
        if ii:
            my_dict[ii] = my_dict[ii]/(4*pi*ii**2/(R_CUT/N_COM))
            x_axis.append(ii)
            y_axis.append(my_dict[ii])
    return x_axis, y_axis

def DictOfBonds(array_mol):
    bonds_dict = []
    c_1 = -1
    c_2 = -1
    for ii in xrange(len(array_mol)-1):
        bonds_dict.append([])
        for jj in xrange(ii+1,len(array_mol)):
            bond_exist = 0
            bond_exist = HydrogeBondExist(array_mol[ii], array_mol[jj])
            if bond_exist:
                bonds_dict[ii].append(jj)
    return bonds_dict
                
def napr(wat_1, wat_2):
    return [wat_1.o_x-wat_2.o_x, wat_1.o_y-wat_2.o_y, wat_1.o_z-wat_2.o_z]

                
def SetOfBonds(array_mol, N_BONDS):
    new_set_mol = []
    bonds_list = DictOfBonds(array_mol)
    
    for i in xrange(len(bonds_list)):
        if len(bonds_list[i]) == N_BONDS:
            new_set_mol.append(array_mol[i])
    

                
#def AngleValFunc(array_mol, N_BONDS):
def AngleValFunc(array_mol):
    bonds_list = DictOfBonds(array_mol)
    angle_list = []
    my_dict_angle = {}
    for i in xrange(len(bonds_list)):
#        if len(bonds_list[i]) == N_BONDS:
#        if len(bonds_list[i]) <100:
            for k in bonds_list[i]:
                vect_ik = napr(array_mol[i], array_mol[k])
                num = bonds_list[i].index(k)+1
                for l in bonds_list[i][num:]:
                    vect_il = napr(array_mol[i], array_mol[l])
                    current = AngleBetweenVectors(vect_ik, vect_il)
                    cur_angle = int(current/width_angle)*width_angle
                    if my_dict_angle.has_key(cur_angle):
                         my_dict_angle[cur_angle] += 1
                    else:
                        my_dict_angle[cur_angle] = 1
#                    angle_list.append(current_angle)
        
#    for ii in sorted(my_dict_angle):
#        if ii:
#            my_dict[ii] /= 2
#            my_dict[ii] = my_dict[ii]/(4*pi*ii**2/(R_CUT/N_COM))
#            x_axis.append(ii)
#            y_axis.append(my_dict[ii]) 
    
    x_axis_angle = sorted(my_dict_angle)
    y_axis_angle = [my_dict_angle[i] for i in x_axis_angle]
    return x_axis_angle, y_axis_angle
                    
def moving_average(x, y, N=1):
    y_new = y
    for jj in xrange(N,len(y)-N):
        a = jj - N
        b = jj+N
        y_new[jj] = sum(y[a:b])/float(len(y[a:b]))
    return x, y_new
                    
                    
def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def AngleBetweenVectors(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))*180/pi
    
            
def SomePlot(x_axis, y_axis, legend):
#        plt.bar(x_axis, y_axis,  width=x_axis[1]-x_axis[0])
        plt.plot(x_axis, y_axis, 'r-')
        plt.legend(legend)
        #labels = [str(i/60) for i in x_axis]
        #plt.xticks(x_axis, labels)
        #my_title = "Radial Distribution Function"
        #plt.title(my_title)
        #plt.savefig(my_title+".png")
    

input_file = "/home/katia/Science/hemoglobin/ver2_1MBN/fixed_prot_1MBN_10ps.pdb"
array_mol = classFactory(input_file)
#x_axis, y_axis = RadialDistrFuncAll(array_mol)
#SomePlot(x_axis, y_axis)

#N_BONDS = 5
#x_ax, y_ax = AngleValFunc(array_mol, N_BONDS)
#если хотим не просто весь array_mol - то надо,например, выделять молекулы, 
#у которых по N_BONDS исходящих связей:

x_dict, y_dict = {}, {}
for N_BONDS in [3,4]:
    new_array_mol = SetOfBonds(array_mol, N_BONDS)

    x_ax, y_ax = AngleValFunc(array_mol)
    x_smooth, y_smooth =  x_ax, y_ax
    for i in xrange(2):
        x_smooth, y_smooth = moving_average(x_smooth, y_smooth, N = 2)
    x_dict[N_BONDS], y_dict[N_BONDS] = x_smooth, y_smooth
SomePlot(x_dict[3], y_dict[3],str(3))
SomePlot(x_dict[4], y_dict[4],str(4))        

#SomePlot(x_smooth, y_smooth)
#plt.figure()
#plt.plot(, '--')


plt.show()
   
