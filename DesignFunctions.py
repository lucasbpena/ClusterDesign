
import random
import copy

# Functions for Design_Principles.py 

def CoreShell(dopant, dopant_n, camadas_dicios, reverse=False):
    '''
    Randomly dopes cluster shells from 
    inside to outside. Can be reversed 
    with 'reverse'
    '''
    composition_c = dopant_n  # Number of dopant atoms
    res_dicio = {} # Resulting dict with {Atom: xyz coordinates}

    camadas_clean = copy.copy(camadas_dicios)
    
    def SubAtom(dicio, key):

        nonlocal composition_c
         
        if composition_c > 0:

            res_dicio[str(dopant + key[2:])] = dicio[key]
            composition_c -= 1
        else:
            res_dicio[str(key)] = dicio[key]

        return 
    
    ###

    if reverse==True:
        
        camadas_clean.reverse()

        for camada in camadas_clean:

            if composition_c >= len(camada):
        
                for key in camada:
                    
                    SubAtom(camada, key)
    
            else:
                
                camada_shuffled = list(camada.keys())

                random.shuffle(camada_shuffled)

                for x in list(range(len(camada))):
                    
                    SubAtom(camada, camada_shuffled[x])
    
    elif reverse==False:
        
        for camada in camadas_clean: 

            if composition_c >= len(camada):
                
                for key in camada:
                    
                    SubAtom(camada, key)
        
            else:
                
                camada_shuffled = list(camada.keys())
                random.shuffle(camada_shuffled)

                for x in list(range(len(camada))):
                    
                    SubAtom(camada, camada_shuffled[x])        

    return res_dicio


def Onion(dopant, dopant_n, camadas_dicios, center_dopant=False):

    '''
    Randomly dopes cluster layers 0 and 2 (center_dopant=True)
    or ony layer 1 (center_dopant=False)
    '''

    composition_c = copy.copy(dopant_n)  # Number of dopant atoms
    res_dicio = {}

    camadas_clean = copy.copy(camadas_dicios)

    for camada in camadas_dicios:
        
        if len(camada) > 0:
            pass
        
        else:
            camadas_clean.remove(camada)

    
    def SubAtom(dicio, key):

        nonlocal composition_c
         
        if composition_c > 0:

            res_dicio[str(dopant + key[2:])] = dicio[key]
            composition_c -= 1

        else:
            res_dicio[str(key)] = dicio[key]
 
        return 


    camada_shuffled1 = list(camadas_clean[1].keys())
    random.shuffle(camada_shuffled1)
    camada_shuffled2 = list(camadas_clean[2].keys())
    random.shuffle(camada_shuffled2)

    
    if center_dopant==True:
        
        for key in camadas_clean[0]:

            SubAtom(camadas_clean[0], key)
        
        for key in camada_shuffled2:

            SubAtom(camadas_clean[2], key)
        
        for key in camada_shuffled1:

            SubAtom(camadas_clean[1], key)
        

    if center_dopant==False:

        for key in camada_shuffled1:

            SubAtom(camadas_clean[1], key)
        
        for key in camada_shuffled2:

            SubAtom(camadas_clean[2], key)

        for key in camadas_clean[0]:

            SubAtom(camadas_clean[0], key)

    return res_dicio


def Segmented(dopant, dopant_n, closest_dist, rest_dist):
    '''
    Selects one atom from entire cluster and replaces the closest ones
    '''

    composition_c = copy.copy(dopant_n)
    res_dicio = {}

    def SubAtom(dicio, key):

        nonlocal composition_c
         
        if composition_c > 0:

            res_dicio[str(dopant + key[2:])] = dicio[key]
            composition_c -= 1

        else:
            res_dicio[str(key)] = dicio[key]
 
        return 

    def WAtom(dicio, key):

        res_dicio[str(key)] = dicio[key]

        return
    
    for key in closest_dist:

        SubAtom(closest_dist, key)
    
    for key in rest_dist:

        WAtom(rest_dist, key)

    return res_dicio

def GetMin(dicio, n):

    '''
    Gets 'n' index minimum distance from dict
    and returns dicts with the minimuns
    and another with non minimuns
    '''
        
    small = {}
    dicio_c = copy.copy(dicio)

    for x in range(n):

        small[min(dicio_c, key=dicio_c.get)] = dicio_c[min(dicio_c, key=dicio_c.get)]
        del dicio_c[min(dicio_c, key=dicio_c.get)]
    
    return small, dicio_c


def GetMax(dicio, n):

    '''
    Gets 'n' maximum distance from dict
    and returns dicts with the maximuns and
    another with non maximums
    '''
        
    large = {}
    dicio_c = dicio

    for x in range(n):

        large[max(dicio_c, key=dicio_c.get)] = dicio_c[max(dicio_c, key=dicio_c.get)]
        del dicio_c[max(dicio_c, key=dicio_c.get)]
    
    return large, dicio_c


def WriteDicio(name, dicio, atom_n):

    '''
    Writes dict (Atoms : xyz coords) to .xyz file.
    '''

    with open(name, 'w') as w:

        w.write(str(atom_n) + '\n')
        w.write('\n')

        for key in dicio:
            
            w.write(key[:2] + "  " + "{:.7f}".format(dicio[key][0]) + "  " + "{:.7f}".format(dicio[key][1]) + "  " + "{:.7f}".format(dicio[key][2]) + "\n")
    
    return


