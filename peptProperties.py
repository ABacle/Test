import mdtraj as mdtraj
import numpy as np

STRUCTURE = ["C", "H", "B", "E", "G", "I", "T", "S", " "]

def deg_to_rad(deg):
    return deg * np.pi / 180

def rad_to_deg(rad):
    return rad * 180 / np.pi

def compute_structure(res_dssp):
    list_structure = []
    for i in range(res_dssp.shape[1]):
        dico_structure = count_structure(res_dssp[:,i])
        list_structure.append(dico_structure)
    return list_structure


def count_structure(arr):
    dic_count = {}
    for c in arr:
        if c not in dic_count.keys():
            dic_count[c] = 1
        else:
            dic_count[c] += 1
    return dic_count

def get_consensus_structure(dic_count_structure):
    consensus_struct = []
    for i in dic_count_structure:
        if len(i.keys()) == 1:
            consensus_struct.append(i.keys()[0])
        else :
            keyMax = keywithmaxval(i)
            if consensus_struct[-1] in keyMax:
                consensus_struct.append(consensus_struct[-1])
    return consensus_struct


def keys_with_max_val(d):
    """ a) create a list of the dict's keys and values;
        b) return the key with the max value"""
    v=list(d.values())
    nbMax =  v.count(max(v))
    k=list(d.keys())
    if nbMax == 1:
        return k[v.index(max(v))]
    else:
        keys = []
        for i in range(nbMax):
            v=list(d.values())
            k=list(d.keys())
            firstKey = k[v.index(max(v))]
            keys.append(firstKey)
            del d[firstKey]
        return keys

def count_structure(consensus_structure):
    dic_count= {}
    for struct in STRUCTURE :
        dic_count[struct] = consensus_structure.count(struct)
    return dic_count

def get_helicity(consensus_structure):
    dic_count_struct = count_structure(consensus_structure)
    len_prot = sum(dic_count_struct.values())
    aa_in_helix = dic_count_struct["H"]
    helicity = aa_in_helix * 100. / len_prot
    return round(helicity,0)

def get_helices(struct):
    nb_h = 0
    #after this threshold the helix is taken into account
    threshold = 4
    dic_struct = {}
    for i in range(len(struct)):
        if struct[i] != struct[i-1]:
            if struct[i-1] in ["H", "G", "I"]:
                dic_struct[struct[i-1]][-1][1] = i-1
            nb_h =1
        else :
            nb_h +=1
            if nb_h == threshold :
                if struct[i] not in dic_struct:
                    dic_struct[struct[i]] = [[i-3, 0]]
                else:
                    dic_struct[struct[i]].append([i-3, 0])
    return dic_struct

def mean(list_angle):
    #Faire la somme de tous les angles
    #puis faire une moyenne
    sum_tot = 0
    length = 0
    for ang in list_angle:
        sum_tot += sum(ang)
        length += len(ang)
    return sum_tot / length

# en degré
def get_phi_psi(traj, dic_helix):
    dic_phsi ={}
    for key in dic_helix.keys():
        for j in range(len(dic_helix[key])):
            start = dic_helix[key][j][0]
            end = dic_helix[key][j][1]
            selection = 'resid %i to %i' % (start, end)
            t_helix = t_prot.atom_slice(topology.select(selection))
            print t_helix
            phis = md.compute_phi(t_helix)[1]
            psis = md.compute_psi(t_helix)[1]
            mean_phi = rad_to_deg(mean(psis))
            mean_psi = rad_to_deg(mean(phis))
            if key in dic_phsi:
                dic_phsi[key].append([mean_phi, mean_psi])
            else :
                dic_phsi[key] = [[mean_phi, mean_psi]]
    return dic_phsi


def get_mean_xy(pos):
    #to change to test the real data
    pos = np.array([
                    [[ 4.81400013 , 2.69099998 , 4.91400003],
                     [ 5.03299999 ,3.15199995 , 5.10699987],
                     [ 5.20100021 , 2.8039999  , 5.51100016],
                     [ 5.24499989 , 2.44300008 , 5.09899998]],
                    [[ 4.67       , 2.56       ,      0    ] ,       #frame 2
                     [ 5.04       , 3.5        ,     0     ],
                     [ 5.4        , 3.0        ,    0      ],
                     [ 5.5        , 2.5        ,   0       ]]
                    ])
                    
                    
                    mean_xy =   sum(pos)/len(pos)
                    #only X/Y columns for all residues
    return mean_xy[:,:2]

def get_xy_helix(t_prot, dic_helix):
    dic_pos_xy = {}
    for key in dic_helix.keys():
        for j in range(len(dic_helix[key])):
            start = dic_helix[key][j][0]
            end = dic_helix[key][j][1]
            selection = 'resid %i to %i and name CB' % (start, end)
            t_helix = t_prot.atom_slice(topology.select(selection))
            pos = t_helix.xyz
            mean_xy = get_mean_xy(pos)
            if key not in dic_pos_xy.keys():
                dic_pos_xy[key] = [mean_xy]
            else:
                dic_pos_xy[key].append([mean_xy])
    return dic_pos_xy


def check_11aa_helix(list_xy):
    #Function qui check si les résidues éloignés de 11 ont des X/Y proche
    # i.e detection d'hélice 11-3
    return 1


def ind_2_leaflets(t_POPC):
    pos_POPC = t_POPC.xyz
    P8_mean = pos_POPC[:,:,2].mean()
    #np.where = array[0] = rows, array[1] = columns
    ind_up_POPC =  np.where(pos_POPC[:,:,2] > P8_mean)[1]
    ind_low_POPC =  np.where(pos_POPC[:,:,2] < P8_mean)[1]
    return (ind_up_POPC, ind_low_POPC)


def mean_z_CA_prot(t_CA):
    z = t_protCA.xyz[:,:,2]
    return z.mean()

def mean_z_leaflets_POPC(t_POPC):
    (ind_up_POPC, ind_low_POPC) = ind_2_leaflets(t_POPC)
    up_POPC_z = pos_POPC[:,ind_up_POPC,2]
    mean_up =  up_POPC.mean()
    low_POPC_z = pos_POPC[:,ind_low_POPC,2]
    mean_low =  low_POPC.mean()
    return (mean_up, mean_low)

def dist_prot_lip(t_CA, t_POPC):
    mean_z_prot = mean_z_CA_prot(t_CA)
    (mean_up, mean_low) = mean_z_leaflets_POPC(t_POPC)
    dist_CA_low = np.abs(mean_low - mean_z_prot)
    dist_CA_up = np.abs(mean_up - mean_z_prot)
    return min(dist_CA_low, dist_CA_up)


def dist_helix_lip(t_prot, t_POPC, dic_helix):
    dic_dist_helix_lip = {}
    for key in dic_helix.keys():
        for j in range(len(dic_helix[key])):
            start = dic_helix[key][j][0]
            end = dic_helix[key][j][1]
            print start, end
            selection = 'resid %i to %i and name CA' % (start, end)
            t_helix = t_prot.atom_slice(topology.select(selection))
            dist = dist_prot_lip(t_helix, t_POPC)
            if key not in dic_dist_helix_lip.keys():
                dic_dist_helix_lip[key] = [dist]
            else:
                dic_dist_helix_lip[key].append([dist])
    return dic_dist_helix_lip


pdb = md.load('last.pdb')
res_dssp = md.compute_dssp(pdb)
prot = topology.select("protein and name CA")
t_prot = pdb.atom_slice(topology.select("protein"))
res_dssp = md.compute_dssp(t_prot)

#Test
dssp_test = np.array(['C', 'C', 'H','H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
                      'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
                      'H', 'H', 'H', 'H', 'H', 'H', 'C'])

dssp_tests = np.vstack((res_dssp[0], dssp_test))
dic_count = compute_structure(dssp_tests)

#Secondary structure and helicity
consensus_struct = get_consensus_structure(dic_count)
helicity = get_helicity(consensus_struct)

#Residu position of the helices and phi/psi mean by helices
helices = get_helices(consensus_struct)
dic_phsi = get_phi_psi(t_prot, helices)

#Attempt to detect 11-mer based on X/Y position superposition
#hyp : if each 11 residues, the X/Y a closed => maybe 11-helix
list_xy = get_xy_helix(t_prot, helices)


#Position of the peptide compared to the bilayer
t_protCA = pdb.atom_slice(topology.select("protein and name CA"))
t_POPC = pdb.atom_slice(topology.select("resname POP and name P8"))
dist_prot_lip(t_protCA, t_POPC)


#Position of each each helix compared to the bilayer (here just one helix)
dist_helix_lip(t_prot, t_POPC, helices)


