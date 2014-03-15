import MDAnalysis.analysis.hbonds as hydbond
import MDAnalysis
import math
import matplotlib.pyplot as plt
import os
import pickle
import re

plt.ion()

# Import .dcd files

# Separate into individual trajectories
# Orient trajectories (know if going from A to B or B to A)
# Align midpoints of reactions

#To investigate:
#   What is the hydroxyethyl group interacting with?
#       KCX?, beta5-beta6 loop?
#   What is balancing the oxygen in the intermediate?
#   What attacks the ester? OH or H2O?
#   Map of proton moving from water
#   
#Can I investigate?:
#   Tautomerization

#Case 1: Doripenem
#   car is 81 CAJ
#   esox is 81 OAI
#   og is 81 OG
#   kcx is 84 OH2
#   prot is W 277 H1
#   wat is W 277 OH2

#Case 2: Imipenem
#   car is 81 C7
#   esox is 81 O62
#   og is 81 OG
#   kcx is 84 OH2
#   prot is W 277 H1
#   wat is W 277 OH2

regex = '[A-Z]+[0-9]+(?=:)'

def getKey(item):
    return item[0]

class BigTraj:
    def __init__(self, universe, is_dori, name, r1, r2):
        self.universe = universe
        self.small_trajs = []
        self.is_dori = is_dori
        self.name = name
        self.r1 = r1
        self.r2 = r2
    def add_traj(self, newtraj):
        self.small_trajs.append(newtraj)
    def analyze_hbonds(self):
        if self.is_dori:
            selestr = "(atom A 81 CAJ) or (atom A 81 OG) or (atom W 277 OH2) or (atom W 277 H1) or " + \
                      "(atom W 277 H2) or (atom A 81 OAI) or (atom A 81 HOI) or (atom A 84 OH1) or " + \
                      "(atom A 84 OH2) or (atom A 81 OAD)"
            OAIstr = "(atom A 81 OAI)"
            OADstr = "(atom A 81 OAD)"
            #Only OAI can be a donor and acceptor, so only create warnings for OAI
            warningOAIacceptor = 'SDR81:OAI'
            warningOAIdonor = 'SDR81:HOI'
            new_donors = ['OAI', 'NAO', 'NAC']
            new_acceptors = ['OAD', 'OAI', 'OAH', 'OAE', 'OH2', 'OH1']
        else:
            selestr = "(atom A 81 C7) or (atom A 81 OG) or (atom W 277 OH2) or (atom W 277 H1) or " + \
                      "(atom W 277 H2) or (atom A 81 O62) or (atom A 81 HO6) or (atom A 84 OH1) or " + \
                      "(atom A 84 OH2) or (atom A 81 O7)"
            OAIstr = "(atom A 81 O62)"
            OADstr = "(atom A 81 O7)"
            warningOAIacceptor = 'SIM81:O62'
            warningOAIdonor = 'SIM81:HO6'
            new_donors = ['O62', 'N24', 'N26']
            new_acceptors = ['O7', 'O62', 'O31', 'O32', 'OH2', 'OH1']
            
        KCXstr = "(atom A 84 OH1) or (atom A 84 OH2)"
        
        #################################
        #GET OAI INFORMATION
        #################################
        hana = hydbond.HydrogenBondAnalysis(self.universe, selection1=OAIstr,
                                            selection2='protein',
                                            donors=new_donors,
                                            acceptors=new_acceptors)
        hana.run()
        h_bond_results = hana.timeseries

        OAIhbonds = {x:{'total':0} for x in range(-700, 700)}

        #Fill OAI dictionary; 0 if bond is not present; 1 if it is
        for r1, frame in zip(self.r1, h_bond_results):
            r1_ind = round(round(r1, 2)*100)
            OAIhbonds[r1_ind]['total'] += 1
            for hbond in frame:
                m = re.search(regex, hbond[2])
                donorstr = m.group(0)
                m = re.search(regex, hbond[3])
                acceptorstr = m.group(0)
                if hbond[2] != warningOAIdonor:
                    if donorstr in OAIhbonds[r1_ind]:
                        OAIhbonds[r1_ind][donorstr] += 1
                    else:
                        OAIhbonds[r1_ind][donorstr] = 1
                if hbond[3] != warningOAIacceptor:
                    if acceptorstr in OAIhbonds[r1_ind]:
                        OAIhbonds[r1_ind][acceptorstr] += 1
                    else:
                        OAIhbonds[r1_ind][acceptorstr] = 1



        #################################
        # GET OAD INFORMATION
        #################################
        hana = hydbond.HydrogenBondAnalysis(self.universe, selection1=OADstr,
                                            selection2='protein',
                                            donors=new_donors,
                                            acceptors=new_acceptors)

        hana.run()
        h_bond_results = hana.timeseries

        OADhbonds = {x:{'total':0} for x in range(-700, 700)}

        #Fill OAD dictionary; 0 if bond is not present; 1 if it is
        for r1, frame in zip(self.r1, h_bond_results):
            r1_ind = round(round(r1, 2)*100)
            OADhbonds[r1_ind]['total'] += 1
            for hbond in frame:
                m = re.search(regex, hbond[2])
                donorstr = m.group(0)
                #OAD can't be a donor, so I don't need to check acceptors
                if donorstr in OADhbonds[r1_ind]:
                    OADhbonds[r1_ind][donorstr] += 1
                else:
                    OADhbonds[r1_ind][donorstr] = 1

        #################################
        # GET KCX INFORMATION
        #
        # Note: This works for now, but once a proton is placed on KCX, this will not count the proton as part of KCX
        # therefore, any interactions that group has will be uncounted
        #################################
        hana = hydbond.HydrogenBondAnalysis(self.universe, selection1=KCXstr,
                                            selection2='protein',
                                            donors=new_donors,
                                            acceptors=new_acceptors)
        hana.run()
        h_bond_results = hana.timeseries

        KCXhbonds = {x:{'total':0} for x in range(-700, 700)}

        #Fill KCX dictionary; 0 if bond is not present; 1 if it is
        for r1, frame in zip(self.r1, h_bond_results):
            r1_ind = round(round(r1, 2)*100)
            KCXhbonds[r1_ind]['total'] += 1
            for hbond in frame:
                m = re.search(regex, hbond[2])
                donorstr = m.group(0)
                #KCX can be a donor if a proton is on it from the water!
                #I am fairly certain this program won't count it, though. Need a better trace...
                if donorstr in KCXhbonds[r1_ind]:
                    KCXhbonds[r1_ind][donorstr] += 1
                else:
                    KCXhbonds[r1_ind][donorstr] = 1
                
        return OAIhbonds, OADhbonds, KCXhbonds


class Trajectory:
    def __init__(self, bigtraj, universe_r1, universe_r2, rxn_beg, rev):
        self.bigtraj = bigtraj
        self.rxn_beg = rxn_beg
        self.rxn_end = rxn_beg + 501
        self.rev = rev
        self.r1 = universe_r1[self.rxn_beg:self.rxn_end]
        self.r2 = universe_r2[self.rxn_beg:self.rxn_end]
        if self.rev:
            self.r1.reverse()
            self.r2.reverse()
        
        j = 0
        i = self.r1[0]
        while i<0:
            j += 1
            i = self.r1[j]
        self.middle = (i, j)
        
    def plot(self):
        x = range(len(self.r1))
        plt.plot(x, self.r1, 'go')
        plt.draw()
    

def get_dist(atom1, atom2):
    return math.sqrt(sum((atom1[i] - atom2[i])**2 for i in range(3)))

def analyze_big_traj(dcd_filepath, psf_filepath, is_dori):
    universe = MDAnalysis.Universe(psf_filepath, dcd_filepath)
    #Dori order is:
        ##0 < Atom 777: name 'OG' of type '74' of resname 'SDR', resid 81 and segid 'A'>
        ##1 < Atom 778: name 'CAJ' of type '32' of resname 'SDR', resid 81 and segid 'A'>
        ##2 < Atom 885: name 'OH1' of type '72' of resname 'KCX', resid 84 and segid 'A'>
        ##3 < Atom 886: name 'OH2' of type '72' of resname 'KCX', resid 84 and segid 'A'>
        ##4 < Atom 3955: name 'OH2' of type '75' of resname 'OH2', resid 277 and segid 'W'>
        ##5 < Atom 3956: name 'H1' of type '4' of resname 'OH2', resid 277 and segid 'W'>
        ##6 < Atom 3957: name 'H2' of type '4' of resname 'OH2', resid 277 and segid 'W'>
        ##
    #Imi order is:
        ##0 < Atom 775: name 'OG' of type '73' of resname 'SIM', resid 81 and segid 'A'>
        ##1 < Atom 778: name 'C7' of type '32' of resname 'SIM', resid 81 and segid 'A'>
        ##2 < Atom 871: name 'OH1' of type '72' of resname 'KCX', resid 84 and segid 'A'>
        ##3 < Atom 872: name 'OH2' of type '72' of resname 'KCX', resid 84 and segid 'A'>
        ##4 < Atom 3941: name 'OH2' of type '75' of resname 'OH2', resid 277 and segid 'W'>
        ##5 < Atom 3942: name 'H1' of type '4' of resname 'OH2', resid 277 and segid 'W'>
        ##6 < Atom 3943: name 'H2' of type '4' of resname 'OH2', resid 277 and segid 'W'>

    if is_dori:
        aoi = universe.selectAtoms("(atom A 81 CAJ) or (atom A 81 OG) or " +
                                   "(atom W 277 OH2) or (atom W 277 H1) or " +
                                   "(atom W 277 H2) or (atom A 84 OH2) or " +
                                   "(atom A 84 OH1)")
    else:
        aoi = universe.selectAtoms("(atom A 81 C7) or (atom A 81 OG) or " +
                                   "(atom W 277 OH2) or (atom W 277 H1) or " +
                                   "(atom W 277 H2) or (atom A 84 OH2) or " +
                                   "(atom A 84 OH1)")
    r1 = []
    #r2 is list of lists:
    #   [dist(H1 - OH1), dist(H1 - OH2), dist(H1-watOH2), dist(H2 - OH1), dist(H2 - OH2), dist(H2-watOH2)]
    r2 = []
    rxnatoms = universe.trajectory.timeseries(aoi)
    for t in range(rxnatoms.shape[1]):
        m = get_dist(rxnatoms[4][t], rxnatoms[1][t])
        b = get_dist(rxnatoms[0][t], rxnatoms[1][t])
        r1.append(b - m)
        h1 = get_dist(rxnatoms[5][t], rxnatoms[2][t])
        h2 = get_dist(rxnatoms[5][t], rxnatoms[3][t])
        h3 = get_dist(rxnatoms[5][t], rxnatoms[4][t])
        h4 = get_dist(rxnatoms[6][t], rxnatoms[2][t])
        h5 = get_dist(rxnatoms[6][t], rxnatoms[3][t])
        h6 = get_dist(rxnatoms[6][t], rxnatoms[4][t])        
        r2.append([h1, h2, h3, h4, h5, h6])
    #501 steps per transition
    big_traj = BigTraj(universe, is_dori, dcd_filepath, r1, r2)
    for i in range(10):
        if r1[i*501] > 0:
            is_rev = True
        else:
            is_rev = False
        t = Trajectory(big_traj, r1, r2, i*501, is_rev)
        big_traj.add_traj(t)
    return big_traj

psfpath = "/data/sguthrie/imivsdori/dori_sim/sim1/template/dori.psf"
rootpath = "/data/sguthrie/imivsdori/dori_sim/sim1/tps_getv2"
isdori=True

##psfpath = "/data/sguthrie/imivsdori/imi_sim/sim1/template/imi.psf"
##rootpath = "/data/sguthrie/imivsdori/imi_sim/sim1/tps_getv/set_r1"
##isdori=False

#cumul_hbonds structure:
#   overall: dictionary indexed by r1 value*100
#       dictionary hashes to list.
#       First element is integer representing the number of simulations that saved at that r1 value
#       Second element is dictionary of OAI interactions
#       Third element is dictionary of OAD interactions
#       Fourth element is dictionary of KCX interactions
cumul_hbonds = {x:[0, {}, {}, {}] for x in range(-700, 700)}
OAIkeys = []
OADkeys = []
KCXkeys = []
trajinfo = []
for root, dirs, files in os.walk(rootpath):
    if len(dirs) == 0:
        print root
        if "tpsv_trajs.dcd" in files:
            dcdpath = root + "/tpsv_trajs.dcd"
            traj = analyze_big_traj(dcdpath, psfpath, isdori)
            OAIhbonds, OADhbonds, KCXhbonds = traj.analyze_hbonds()
            trajinfo.append(traj)
            for r1_val in OAIhbonds:
                assert OAIhbonds[r1_val]['total'] == OADhbonds[r1_val]['total']
                assert OAIhbonds[r1_val]['total'] == KCXhbonds[r1_val]['total']
                cumul_hbonds[r1_val][0] += OAIhbonds[r1_val]['total']
                for key in OAIhbonds[r1_val]:
                    if key != 'total':
                        if key not in OAIkeys:
                            OAIkeys.append(key)
                        if key in cumul_hbonds[r1_val][1]:
                            cumul_hbonds[r1_val][1][key] += OAIhbonds[r1_val][key]
                        else:
                            cumul_hbonds[r1_val][1][key] = OAIhbonds[r1_val][key]
                            
                for key in OADhbonds[r1_val]:
                    if key != 'total':
                        if key not in OADkeys:
                            OADkeys.append(key)
                        if key in cumul_hbonds[r1_val][2]:
                            cumul_hbonds[r1_val][2][key] += OADhbonds[r1_val][key]
                        else:
                            cumul_hbonds[r1_val][2][key] = OADhbonds[r1_val][key]

                for key in KCXhbonds[r1_val]:
                    if key != 'total':
                        if key not in KCXkeys:
                            KCXkeys.append(key)
                        if key in cumul_hbonds[r1_val][3]:
                            cumul_hbonds[r1_val][3][key] += KCXhbonds[r1_val][key]
                        else:
                            cumul_hbonds[r1_val][3][key] = KCXhbonds[r1_val][key]


OAIlist = []
OADlist = []
KCXlist = []
for r1 in cumul_hbonds:
    total = cumul_hbonds[r1][0]
    if total != 0:
        total = float(total)
        #Making sure I have some simulations in that r1 bracket
        mini_list = [r1/100.]
        for key in OAIkeys:
            if key in cumul_hbonds[r1][1]:
                prob = cumul_hbonds[r1][1][key]/total
            else:
                prob = 0
            mini_list.append(prob)
            #TODO: check that this is correct variance
            mini_list.append(prob*(1-prob))
        OAIlist.append(mini_list)
    
        mini_list = [r1/100.]
        for key in OADkeys:
            if key in cumul_hbonds[r1][2]:
                prob = cumul_hbonds[r1][2][key]/total
            else:
                prob = 0
            mini_list.append(prob)
            #TODO: check that this is correct variance
            mini_list.append(prob*(1-prob))
        OADlist.append(mini_list)

        mini_list = [r1/100.]
        for key in KCXkeys:
            if key in cumul_hbonds[r1][3]:  
                prob = cumul_hbonds[r1][3][key]/total
            else:
                prob = 0
            mini_list.append(prob)
            #TODO: check that this is correct variance
            mini_list.append(prob*(1-prob))
        KCXlist.append(mini_list)
        

sorted_OAI = sorted(OAIlist, key=getKey)
sorted_OAD = sorted(OADlist, key=getKey)
sorted_KCX = sorted(KCXlist, key=getKey)

def plot_2D (sorted_list, key_list, intoTS, outofTS, figname):
    """ intoTS and outofTS are TS basin definitions. Should be real numbers
        figname should be string to write figure name to. Should include .png"""
    annoying = zip(*(sorted_list))
    r1_vals, annoying = annoying[0], annoying[1:]
    legend_handles = []
    plt.figure(figsize=(12,10), dpi=100)
    for key in key_list:
        avg, err = annoying[0], annoying[1]
        if len(annoying) > 2:
            annoying = annoying[2:]
        h = plt.errorbar(r1_vals, avg, yerr=err, fmt='-o')
        legend_handles.append(h)
    plt.axvline(x=intoTS, color='k')
    plt.axvline(x=outofTS, color = 'k')
    plt.grid(True)
    plt.axis([-7, 7, -0.5, 1.5])
    plt.legend(legend_handles, key_list)
    plt.draw()
    plt.savefig(figname)

plot_2D(sorted_OAI, OAIkeys, -1, 1.5, 'DoriOAI_interactions.png')
plot_2D(sorted_OAD, OADkeys, -1, 1.5, 'DoriOAD_interactions.png')
plot_2D(sorted_KCX, KCXkeys, -1, 1.5, 'DoriKCX_interactions.png')


if rootpath == "/data/sguthrie/imivsdori/dori_sim/sim1/tps_getv2":
    output = open('Dori_interactions.pkl', 'wb')
    #Note: can't pickle AtomGroups
    #pickle.dump(trajinfo, output, -1)
    pickle.dump(sorted_OAI, output, -1)
    pickle.dump(sorted_OAD, output, -1)
    pickle.dump(sorted_KCX, output, -1)
    output.close()

##inp = open('DoriOAI_interactions.pkl', 'rb')
##whee = pickle.load(inp)
##whee2 = pickle.load(inp)
##whee3 = pickle.load(inp)
##inp.close()


    


