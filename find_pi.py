#!/usr/bin/env python

# find face face pi stacks (the pi-danza motif, in pdbs)
# vers 0.1

import sys
import numpy as np
import itertools
import math
import collections

headless = False
testfile = sys.argv[1]
if '--h' in sys.argv:
    headless = True
vers = 0.1
dthresh = 6
angthresh = 20
def read_pdb_find_aromatics(file):
    
    '''find all the aromatics get the coordinates of all six carbons
    make sure all six are there '''
    aromatics = {}      #dic {name:{CG:array(coords),CD1:array(coords), and etc...}} 
    data = open(file,'r').readlines()
    for line in data:
        if line[0:4] == 'ATOM' and line[17:20] in ['PHE','TYR'] and line[13:16] in ['CG ','CD1','CD2','CE1','CE2','CZ ']:
            aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
            if aaname not in aromatics:
                aromatics[aaname] = {}
            aromatics[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
        if line[0:4] == 'ATOM' and line[17:20] == 'TRP' and line[13:16] in ['CD2','CE3','CZ3','CH2','CZ2','CE2']:
            aaname ='{0}.{1}{2}'.format(line[21],line[17:20],line[22:26])
            if aaname not in aromatics:
                aromatics[aaname] = {}
            aromatics[aaname][line[13:16].strip(' ')] = np.array([float(line[31:38]),float(line[39:46]),float(line[47:54])])
    
    baddies = []
    for i in aromatics:                     # make sure all six carbons are present
        if len(aromatics[i]) != 6:
            baddies.append(i)
    for i in baddies:
        if headless == False:
            print('dropped {0}: missing carbon(s)'.format(i))
        del aromatics[i]
    return(aromatics)
            

def find_ring_center(ring):
    '''find the center of an aromatic ring input dic entry in format aaname:{CG:array(coords),CD1:array(coords), and etc...}'''

    if 'CG' in ring:   # for phe-tyr
        center = (ring['CG']+ring['CD1']+ring['CD2']+ring['CE1']+ring['CE2']+ring['CZ'])/6

    else:               # for trp
        center = (ring['CD2']+ring['CE2'])/2
    
    ''' calculate face of the ring by vector normal to center - CG and Center - CD1 for phe/tyr and center - CD2 and Center - CE2 for trp'''
    if 'CG ' in ring:
        vec1 = np.array([center,ring['CG']])
        vec2 = np.array([center,ring['CD1']])
    else:
        vec1 = np.array([center,ring['CD1']])
        vec2 = np.array([center,ring['CD2']])
    
    v1 = np.subtract(vec1[1,...],vec1[0,...])
    v2 = np.subtract(vec2[1,...],vec2[0,...])
    face_norm = np.cross(v1,v2)

    return(center,face_norm)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    if vector.all() == 0:
        vector = np.array([0.001,0.001,0.001])
    return (vector/np.linalg.norm(vector))

def angle_between(v1, v2):
    """ Return the angle in degrees between two vectors"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    if angle == 0:
        angle = 0.000001
    if angle < 90:
        return(angle)
    else:
        return(180.0 - angle)


data = read_pdb_find_aromatics(testfile)

# get the centers and face normal vectors of all aromatics
centers = {} # center {coords:aaname}
face_norms = {}
for i in data:
    centers[i] = (find_ring_center(data[i])[0])
    face_norms[i] = (find_ring_center(data[i])[1])

# find the different possible combinations of rings
combinations = list(itertools.combinations(centers,2))

# calculate the distances between all of the centers
distances = {}
for i in combinations:
    combo = '{0}-{1}'.format(i[0],i[1])
    distances[combo] = np.linalg.norm(centers[i[0]]-centers[i[1]])

## find the distances less than threshold and find chains of multiple within threshold pairs
threshed = []
threshed_pairs = []
threshed_aas = {}

## get the interactions within thresholds
# find all interactions within distance threshold
for i in distances:
    if distances[i] < dthresh:
        threshed.append(i)
        threshed_pairs.append(i.split('-'))


# calculate the angles between the threshed pairs
final_pairs = []
if headless == False:
    print('AllPairs: Pair,Distance,Angle')
for pair in threshed_pairs:
    pair_ang = (angle_between(face_norms[pair[0]], face_norms[pair[1]]))    
    combo = '{0}-{1}'.format(pair[0],pair[1])   # make nmae to retreve distance
    pair_dist = distances[combo]
    if headless == False:
        print(pair,pair_dist,pair_ang)
    pair_screen = (pair,pair_ang)
    if pair_screen[1] < angthresh:
        final_pairs.append(pair_screen[0])

# make a bild file to show the pi stacks

bildfile = open('bildfile.bld','w')
for i in final_pairs:
    bildfile.write('.color red \n.arrow {0} {1} {2} {3} {4} {5} 0.1 0.1 0.1\n'.format(centers[i[0]].item(0),centers[i[0]].item(1),centers[i[0]].item(2),centers[i[1]].item(0),centers[i[1]].item(1),centers[i[1]].item(2)))
bildfile.close()

# get the individual aas
for i in final_pairs:        
    for j in i:
        if j not in threshed_aas:
            threshed_aas[j] = 0

# find how may hits each aa has
for aa in threshed_aas:
    for pair in final_pairs:
        if aa in pair:
            threshed_aas[aa] +=1

#make the intial chain nodes
chain_nodes = []
for i in threshed_aas:
    if threshed_aas[i] > 1 :
        chain_nodes.append(i)

chains = {}
for i in chain_nodes:
    if i not in chains:
        chains[i] = []
    for pair in final_pairs:
        if i in pair:
            for aa in pair:
                if aa not in chains[i]:
                    chains[i].append(aa)
chain_seeds = []
for i in chains:
    chain_seeds.append(chains[i])

# iterate over list of seeds and link connected chains. 
len_l = len(chain_seeds)
i = 0
while i < (len_l - 1):
    for j in range(i + 1, len_l):
        i_set = set(chain_seeds[i])
        j_set = set(chain_seeds[j])
        if len(i_set.intersection(j_set)) > 0:
            chain_seeds.pop(j)
            chain_seeds.pop(i)
            ij_union = list(i_set.union(j_set))
            chain_seeds.append(ij_union)
            len_l -= 1
            i -= 1
            break

    i += 1

if headless == False:
    print('\nMultiple stacks')
chain_counts = []
for i in chain_seeds:
    if headless == False:
        print(i)
    chain_counts.append(len(i))
if len(chain_counts) < 1:
    if headless == False:
        print('None found')

# stats on multi-stacks
count = collections.Counter(chain_counts)
if headless == False:
    if len(count) > 0:
            print('Length\tcount')
    for i in count:
        print('{0: 3}\t{1: 2}'.format(i,count[i]))

def get_numbers(list):
    ''' gets residue numbers for writing chimera commands '''
    aanos = []
    for i in list:
        aanos.append('{0}.{1}'.format(i[-4:].strip(' '),i[0]))
    return (','.join(aanos))
## write out command for chimera 
if headless == False:
    print('\n**** chimera command ****\n')

chimera_commands = []
colors = ('red','yellow','blue','green','purple','pink')
colorn = 0
for i in final_pairs:
    chimera_commands.append(('color {1} #0:{0} ; disp #0:{0}'.format(get_numbers(i),colors[colorn])))
    colorn+=1
    if colorn == 6:
        colorn = 0
for i in chain_seeds:
    chimera_commands.append(('color {1} #0:{0} ; disp #0:{0}'.format(get_numbers(i),colors[colorn])))
    colorn+=1
    if colorn == 6:
        colorn = 0
if headless == False:
    print ('color white ; ~disp ; {0}'.format(' ; '.join(chimera_commands)))
if headless == True and len(count) >= 1:
    outname = '{}.txt'.format(testfile.split('.')[0].replace('/','-'))
    output = open('hits/{0}'.format(outname),'w')
    output.write('color white ; ~disp ; {0}'.format(' ; '.join(chimera_commands)))
