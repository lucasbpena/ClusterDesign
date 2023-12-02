import kmeans1d
from scipy.spatial import distance

import os
import sys

from os.path import join as ospj
from DesignFunctions import * 


# Get Inputs - 
# ~$ python3 sys.argv[0] (DesignPrinciples.py) sys.argv[1] (filename) ...

filename = str(sys.argv[1]) # .xyz file in which permutations with dopant atoms will be conducted

dopant = str(sys.argv[2])  # Dopant atom to insert over unary structure (.xyz file)

dopant_n = int(sys.argv[3]) # Number of dopant atoms

permutations = int(sys.argv[4])  # How many permutations to generate for each design principle

## Show
print(f"\n=== Showing Inputs ===")
print(f"\n    .xyz File  =  {filename}")
print(f"    Dopant  =  {dopant}")
print(f"    Number of dopant atoms =  {str(dopant_n)}")
print(f"    Permutations per design =  {permutations}")
input("\n    press return...")

# Reading .xyz basis file for permutations
diretorio = os.path.dirname(os.path.realpath(__file__)) # This script directory

## Get structure dict 
with open(ospj(diretorio, filename), 'r') as xyz:

    atom_n = int(xyz.readline())
    empty = xyz.readline()

    coordinates = {}

    for c, line in enumerate(xyz):
        
        atom,x,y,z = line.split()
        
        atom_og = atom
        
        coordinates[str(atom) + str(c)] = [float(x), float(y), float(z)] # Gets dict with Atom1, Atom2, Atom3 ... : Coordinates

# Makes directory for generated structures
try:
    os.mkdir("generated_structures")
except:
    pass

# Atom distances to cluster center

## Calculating structure geometric center
coordslist = list(coordinates.values())
x = []
y = []
z = []
for coords in coordslist:
    x.append(coords[0])
    y.append(coords[1])
    z.append(coords[2])

media_coord = [sum(x)/len(x), sum(y)/len(y), sum(z)/len(z)] # Average of coordinates (structure geometric center)
distances = {}

## Calculating atom distances to geometric center
for atom in coordinates.keys():
    distances[atom] = distance.euclidean(media_coord, coordinates[atom]) # Dictionary for each atom euclidean distance to structure center


# kmeans1d package - Cluster layers    
print("\n=== Generated Structure Layers === ")
print("    (kmeans1d package)")

kmeans_k = 3 # k number - defines how many layers the cluster has 

distances_clean = list(distances.values())
clusters, centroids = kmeans1d.cluster(distances_clean, kmeans_k) # Running 1D k-means for grouping atoms based on their distances to structure center

## Assign each coord to each layer (camada)
camada0 = {}
camada1 = {}
camada2 = {}

camadas = [camada0, camada1, camada2]

coordinates_keys = list(coordinates.keys())

for c, cluster in enumerate(clusters):

    if cluster == 0:
        camada0[coordinates_keys[c]] = coordinates[coordinates_keys[c]]

    elif cluster == 1: 
        camada1[coordinates_keys[c]] = coordinates[coordinates_keys[c]]

    elif cluster == 2: 
        camada2[coordinates_keys[c]] = coordinates[coordinates_keys[c]]

## View generated k-means clusters (structure layers)
print(f"\n    Atoms per layer:\n")
print(f"    Layer 0: {len(camada0)}")
print(f"    Layer 1: {len(camada1)}")
print(f"    Layer 2: {len(camada2)}")
print(f"\n    Centroids distances: {centroids}")
input("\n    press return...")

# Design Principles Generation
## Core-shell
dicios_core = {}

for x in list(range(permutations)):

    dicios_core[f'{dopant}{dopant_n}_core_{x}'] = CoreShell(dopant, dopant_n, camadas)

for dicio in dicios_core:

    WriteDicio(ospj('generated_structures', dicio + '.xyz'), dicios_core[dicio], atom_n)

## Shell-core
dicios_shell = {}

for x in list(range(permutations)):

    dicios_shell[f'{dopant}{dopant_n}_shell_{x}'] = CoreShell(dopant, dopant_n, camadas, reverse=True)

for dicio in dicios_shell:

    WriteDicio(ospj('generated_structures',  dicio + '.xyz'), dicios_shell[dicio], atom_n)

## Onion
### Center Dopant
dicios_onion_c ={}

for x in list(range(permutations)):

    dicios_onion_c[f'{dopant}{dopant_n}_Onion_c_{x}'] = Onion(dopant, dopant_n, camadas, center_dopant=True)

for dicio in dicios_onion_c:

    WriteDicio(ospj('generated_structures',  dicio + '.xyz'), dicios_onion_c[dicio], atom_n)

### No center Dopant
dicios_onion ={}

for x in list(range(permutations)):

    dicios_onion[f'{dopant}{dopant_n}_Onion_nc_{x}'] = Onion(dopant, dopant_n, camadas)

for dicio in dicios_onion:

    WriteDicio(ospj('generated_structures',  dicio + '.xyz'), dicios_onion[dicio], atom_n)


## Segmented
folga = 5 # Number of additional atoms to include in permutations for segmented generation
### Segmented generation permutates the closest "dopant_n" atoms randomly. It's possible to consider more atoms with this parameter

for x in list(range(permutations)):

    sample = list(coordinates.keys())[random.randint(0, len(coordinates) - 1)]
    
    distances_collections = {}
    distances_seg = {}

    for atomd in coordinates:    
        
        distances_seg[atomd] = distance.euclidean(coordinates[atom], coordinates[atomd])

    distances_collections[c] = distances_seg


    for dist_dicio in distances_collections:
        
        small, dicio_rest = GetMin(distances_collections[dist_dicio], dopant_n + folga)
        small_values_sorted = sorted(small.values())

        small_sorted = {}

        for i in small_values_sorted:
            for k in small.keys():
                if small[k] == i:
                    small_sorted[k] = small[k]
                    break
        
        small_res ={}
        rest_res = {}

        for key in small_sorted:
            small_res[key] = coordinates[key]
        
        for key in dicio_rest:
            rest_res[key] = coordinates[key]

        res_dicio = Segmented(dopant, dopant_n, small_res, rest_res)

        WriteDicio(ospj('generated_structures',  f'{dopant}{dopant_n}_segmented_{x}' + '.xyz'), res_dicio, atom_n)