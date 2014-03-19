import MDAnalysis.analysis.hbonds as hydbond
import MDAnalysis
import math
import matplotlib.pyplot as plt
import os
import pickle
import re

plt.ion()


global debug
debug = False

global testing
testing = False

def print_hbonds(b, inter, r1):
    if debug:
        print b
        print "Is classified to be interacting with " + inter
        print "At r1 = ", r1

def incr_dict(dictionary, key):
    dictionary[key] = 1
##    if key in dictionary:
##        dictionary[key] += 1
##    else:
##        dictionary[key] = 1

    
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
    def __init__(self, universe, is_dori, name, r1, r2, moltype):
        self.universe = universe
        self.small_trajs = []
        self.is_dori = is_dori
        self.name = name
        self.r1 = r1
        self.r2 = r2
        self.moltype = moltype
        self.kcx_attached, self.og_attached = self.analyze_protons()
    def add_traj(self, newtraj):
        self.small_trajs.append(newtraj)
        
    def analyze_protons(self):
        #r2extra is list of lists:
        #   [dist(H1 - OH1), dist(H1 - OH2), dist(H1-watOH2),
        #    dist(H2 - OH1), dist(H2 - OH2), dist(H2-watOH2),
        #    dist(H1-OG), dist(H2-OG)]
        prot_kcx = []
        prot_og = []
        r_cov = 1.31
        for h11, h12, h1w, h21, h22, h2w, h1o, h2o in self.r2:
            kcx = []
            kcx_attached = False
            og = []
            og_attached = False
            if h11 <= r_cov:
                kcx_attached = True
                kcx.append('(H1-OH1)')
            if h12 <= r_cov:
                kcx_attached = True
                kcx.append('(H1-OH2)')
            if h21 <= r_cov:
                kcx_attached = True
                kcx.append('(H2-OH1)')
            if h22 <= r_cov:
                kcx_attached = True
                kcx.append('(H2-OH2)')
            if h1o <= r_cov:
                og_attached = True
                og.append('(H1-OG)')
            if h2o <= r_cov:
                og_attached = True
                og.append('(H2-OG)')
            toappend = [kcx_attached] + [x for x in kcx]
            prot_kcx.append(toappend)
            toappend = [og_attached] + [x for x in og]
            prot_og.append(toappend)
        foo = [prot_kcx, prot_og]
        return foo
    def analyze_hbonds(self):
        if self.is_dori:
            selestr = "(atom A 81 CAJ) or (atom A 81 OG) or (atom W 277 OH2) or (atom W 277 H1) or " + \
                      "(atom W 277 H2) or (atom A 81 OAI) or (atom A 81 HOI) or (atom A 84 OH1) or " + \
                      "(atom A 84 OH2) or (atom A 81 OAD)"
            examinestr = "(atom A 81 OAI) or (atom A 81 OAD) or (atom A 84 OH1) or (atom A 84 OH1) " + \
                         "or (atom A 84 OH2) or (atom W 277 OH2) or (atom A 81 OG)"
            #OG is already a donor and acceptor
            new_donors = ['OAI', 'NAO', 'NAC', 'OH1', 'OH2'] 
            new_acceptors = ['OAD', 'OAI', 'OAH', 'OAE', 'OH2', 'OH1']
            OAIacceptor = self.moltype + '81:OAI'
            OAIdonor = self.moltype + '81:HOI'
            OADacceptor = self.moltype + '81:OAD'
        else:
            selestr = "(atom A 81 C7) or (atom A 81 OG) or (atom W 277 OH2) or (atom W 277 H1) or " + \
                      "(atom W 277 H2) or (atom A 81 O62) or (atom A 81 HO6) or (atom A 84 OH1) or " + \
                      "(atom A 84 OH2) or (atom A 81 O7)"
            examinestr = "(atom A 81 O62) or (atom A 81 O7) or (atom A 84 OH1) or (atom A 84 OH1) " + \
                         "or (atom A 84 OH2) or (atom W 277 OH2) or (atom A 81 OG)"
            new_donors = ['O62', 'N24', 'N26', 'OH1', 'OH2']
            new_acceptors = ['O7', 'O62', 'O31', 'O32', 'OH2', 'OH1']
            OAIacceptor = self.moltype + '81:O62'
            OAIdonor = self.moltype + '81:HO6'
            OADacceptor = self.moltype + '81:O7'
        KCX1acceptor = 'KCX84:OH1'
        KCX2acceptor = 'KCX84:OH2'
        WATacceptor = 'OH2277:OH2'
        WATdonor1 = 'OH2277:H1'
        WATdonor2 = 'OH2277:H2'
        OGacceptor = self.moltype + '81:OG'
        #WATdonors can be on KCX or OG, not guaranteed to be on water. 
        
        #################################
        #GET HYDROGEN BONDING INFORMATION
        #################################
        hana = hydbond.HydrogenBondAnalysis(self.universe, selection1=examinestr,
                                            selection2='all',
                                            donors=new_donors,
                                            acceptors=new_acceptors, angle=150.0)
        hana.run()
        h_bond_results = hana.timeseries

        total = {t:[0, 0, 0] for t in range(-700, 700)}
        OAIhbonds = {t:{} for t in range(-700, 700)}
        OADhbonds = {t:{} for t in range(-700, 700)}
        KCXhbonds = {t:{} for t in range(-700, 700)}
        WAThbonds = {t:{} for t in range(-700, 700)}
        OGhbonds = {t:{} for t in range(-700, 700)}

        r_cov = 1.31
        #Fill dictionaries; 0 if bond is not present; 1 if it is
        for tstep, foo in enumerate(zip(self.kcx_attached, self.og_attached, self.r1, h_bond_results)):
            prot_kcx, prot_og, r1, frame = foo
            r1 = int(round(r1, 2)*100)
            total[r1][0] += 1
            if prot_kcx[0]:
                total[r1][1] += 1
            if prot_og[0]:
                total[r1][2] += 1
            #r2 is list of lists:
            #   [dist(H1 - OH1), dist(H1 - OH2), dist(H1-watOH2),
            #    dist(H2 - OH1), dist(H2 - OH2), dist(H2-watOH2),
            #    dist(H1-OG), dist(H2-OG)]
            for hbond in frame:
                m = re.search(regex, hbond[2])
                donorstr = m.group(0)
                m = re.search(regex, hbond[3])
                acceptorstr = m.group(0)
                #Go through donors and add acceptor to correct dictionary
                if hbond[2] == OAIdonor:
                    print_hbonds(hbond, "OAIdonor", r1/100.)
                    incr_dict(OAIhbonds[r1], acceptorstr)
                if hbond[2] == WATdonor1:
                    if self.r2[tstep][0] <= r_cov or self.r2[tstep][1] <= r_cov:
                        print_hbonds(hbond, "KCXdonor", r1/100.)
                        incr_dict(KCXhbonds[r1], acceptorstr)
                    if self.r2[tstep][2] <= r_cov:
                        print_hbonds(hbond, "WATdonor", r1/100.)
                        incr_dict(WAThbonds[r1], acceptorstr)
                    if self.r2[tstep][6] <= r_cov:
                        print_hbonds(hbond, "OGdonor", r1/100.)
                        incr_dict(OGhbonds[r1], acceptorstr)
                    #Should not ever not be in the case where it's flying in the ether, since it wouldn't be
                    #counted by the hbonds.py algorithm
                        
                if hbond[2] == WATdonor2:
                    if self.r2[tstep][3] <= r_cov or self.r2[tstep][4] <= r_cov:
                        print_hbonds(hbond, "KCXdonor", r1/100.)
                        incr_dict(KCXhbonds[r1], acceptorstr)
                    if self.r2[tstep][5] <= r_cov:
                        print_hbonds(hbond, "WATdonor", r1/100.)
                        incr_dict(WAThbonds[r1], acceptorstr)
                    if self.r2[tstep][7] <= r_cov:
                        print_hbonds(hbond, "OGdonor", r1/100.)
                        incr_dict(OGhbonds[r1], acceptorstr)
                    #Should not ever not be in the case where it's flying in the ether, since it wouldn't be
                    #counted by the hbonds.py algorithm
                        
                #Go through acceptors and add donor to correct dictionary
                if hbond[3] == OAIacceptor:
                    print_hbonds(hbond, "OAIacceptor", r1/100.)
                    incr_dict(OAIhbonds[r1], donorstr)
                if hbond[3] == OADacceptor:
                    print_hbonds(hbond, "OADacceptor", r1/100.)
                    incr_dict(OADhbonds[r1], donorstr)
                if (hbond[3] == KCX1acceptor) or (hbond[3] == KCX2acceptor):
                    print_hbonds(hbond, "KCXacceptor", r1/100.)
                    incr_dict(KCXhbonds[r1], donorstr)
                if hbond[3] == WATacceptor:
                    print_hbonds(hbond, "WATacceptor", r1/100.)
                    incr_dict(WAThbonds[r1], donorstr)
                if hbond[3] == OGacceptor:
                    print_hbonds(hbond, "OGacceptor", r1/100.)
                    incr_dict(OGhbonds[r1], donorstr)
        return total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds


class Trajectory:
    '''BUG: assumes path is 501 steps long! '''
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

def analyze_big_traj(dcd_filepath, psf_filepath, is_dori, moltype):
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
    #   [dist(H1 - OH1), dist(H1 - OH2), dist(H1-watOH2), dist(H2 - OH1), dist(H2 - OH2), dist(H2-watOH2), dist(H1-OG), dist(H2-OG)]
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
        h7 = get_dist(rxnatoms[5][t], rxnatoms[0][t])
        h8 = get_dist(rxnatoms[6][t], rxnatoms[0][t])
        r2.append([h1, h2, h3, h4, h5, h6, h7, h8])
    #501 steps per transition
    big_traj = BigTraj(universe, is_dori, dcd_filepath, r1, r2, moltype)
    # Don't really need to work on small trajectories. Just taking up memory
##    for i in range(10):
##        if r1[i*501] > 0:
##            is_rev = True
##        else:
##            is_rev = False
##        t = Trajectory(big_traj, r1, r2, i*501, is_rev)
##        big_traj.add_traj(t)
    return big_traj


def add_to_list(listtoaddto, keylist, rpoint, cumul_hbonds, total, i):
    '''Iterates through keylist, so keylist needs to be ordered before calling add_to_list'''
    z = 1.96
    zsq = 3.84
    mini_list = [rpoint/100.]
    for key in keylist:
        if key in cumul_hbonds[rpoint][i]:
            if total == 0:
                prob = 0
            else:
                prob = cumul_hbonds[rpoint][i][key]/total
        else:
            prob = 0
        #mini_list = [probability, wilson score upper, wilson score lower]
        mini_list.append(prob)
        #mini_list.append(z*math.sqrt(prob*(1-prob)/total)) # Normal approx interval. Not good when prob near 0 or 1
        if total != 0:
            wscorepm = z*math.sqrt(prob*(1-prob)/total + zsq/(4*total**2))
            wmult = 1/(1+zsq/total)
            mini_list.append(wmult*(prob + zsq/(2*total) + wscorepm))
            mini_list.append(wmult*(prob + zsq/(2*total) - wscorepm))
        else:
            mini_list.append(0)
            mini_list.append(0)
        
    listtoaddto.append(mini_list)

def plot_allkeys (sorted_list, key_list, intoTS, outofTS, figname):
    """ intoTS and outofTS are TS basin definitions. Should be real numbers
        figname should be string to write figure name to. Should include .png"""
    annoying = zip(*(sorted_list))
    r1_vals, annoying = annoying[0], annoying[1:]
    legend_handles = []
    plt.figure(figsize=(12,10), dpi=100)
    for key in key_list:
        avg, uperr, lowerr = annoying[0], annoying[1], annoying[2]
        asymmerr = [lowerr, uperr]
        if len(annoying) > 3:
            annoying = annoying[3:]
        h = plt.errorbar(r1_vals, avg, yerr=asymmerr, fmt='-o')
        legend_handles.append(h)
    plt.axvline(x=intoTS, color='k')
    plt.axvline(x=outofTS, color = 'k')
    plt.grid(True)
    plt.axis([-7, 7, -0.5, 1.5])
    plt.legend(legend_handles, key_list)
    plt.title(figname)
    plt.draw()
    plt.savefig(figname)


infolist = []

psfpaths = ["/data/sguthrie/imivsdori/dori_sim/sim1/template/dori.psf", 
            "/data/sguthrie/imivsdori/imi_sim/sim1/template/imi.psf"]
rootpaths = ["/data/sguthrie/imivsdori/dori_sim/sim1/tps_getv2", 
             "/data/sguthrie/imivsdori/imi_sim/sim1/tps_getv"]
isdoris = [True, False]
moltypes = ['SDR', 'SIM']
picklefiles = ['TPSDori_interactions.pkl', 'TPSImi_interactions.pkl']
for i in range(2):
    tmp = [psfpaths[i], rootpaths[i], isdoris[i], moltypes[i],  picklefiles[i]]
    infolist.append(tmp)

#cumul_hbonds structure:
#   overall: dictionary indexed by r1 value*10
#       dictionary hashes to list.
#       First element is integer representing the number of simulations that saved at that r1 value
#       Second element is dictionary of OAI interactions
#       Third element is dictionary of OAD interactions
#       Fourth element is dictionary of KCX interactions
#       Fifth element is dictionary of Water interactions
#       Sixth element is dictionary of OG interactions
#       Seventh element is integer representing the number of simulations that had a proton on KCX
#       Eigth element is integer representing the number of simulations that had a proton on OG

if debug or testing:
    psfpath, rootpath, isdori, moltype, pfile = infolist[1]
    keylist = [[] for x in range(5)]
    cumul_hbonds = {t:[0, {}, {}, {}, {}, {}, 0, 0] for t in range(-700,700)}
    rootpath = rootpath + "/set_r1"
    for root, dirs, files in os.walk(rootpath):
        if len(dirs) == 0:
            print root
            if "tpsv_trajs.dcd" in files:
                dcdpath = root + "/tpsv_trajs.dcd"
                traj = analyze_big_traj(dcdpath, psfpath, isdori, moltype)
                total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds = traj.analyze_hbonds()
                hbondlists = [OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds]
                for r1_val in total:
                    if total[r1_val][0] != 0:
                        cumul_hbonds[r1_val][0] += total[r1_val][0]
                        cumul_hbonds[r1_val][6] += total[r1_val][1]
                        cumul_hbonds[r1_val][7] += total[r1_val][2]
                        for i in range(5):
                            for inter in hbondlists[i][r1_val]:
                                if inter not in keylist[i]:
                                    keylist[i].append(inter)
                                block_incr_dict(hbondlists[i][r1_val][inter], cumul_hbonds[r1_val][i+1], inter)
                            keylist[i].sort()
    # cumul_lists: OAI, OAD, KCX, WAT, OG
    cumul_lists = [[] for x in range(5)]
    
    for rpoint in cumul_hbonds:
        total = cumul_hbonds[rpoint][0]
        #will want an indicator to know if total was zero...
        total = float(total)
        for i in range(5):
            add_to_list(cumul_lists[i], keylist[i], rpoint, cumul_hbonds, total,i+1)

    sorted_cumul_lists = []
    for i in range(5):
        sorted_cumul_lists.append(sorted(cumul_lists[i]))
else:
    all_cumul_lists = []
    all_key_lists = []
    for foo in range(2):
        psfpath, rootpath, isdori, moltype, pfile = infolist[foo]
        print rootpath
        try:
            inp = open(pfile, 'rb')
            cumul_hbonds = pickle.load(inp)
            keylist = pickle.load(inp)
            all_key_lists.append(keylist)
            inp.close()
        except IOError:
            # keylists: OAI, OAD, KCX, WAT, OG
            keylist = [[] for x in range(5)]

            cumul_hbonds = {t:[0, {}, {}, {}, {}, {}, 0, 0] for t in range(-700,700)}
            for root, dirs, files in os.walk(rootpath):
                if len(dirs) == 0:
                    print root
                    if "tpsv_trajs.dcd" in files:
                        dcdpath = root + "/tpsv_trajs.dcd"
                        traj = analyze_big_traj(dcdpath, psfpath, isdori, moltype)
                        total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds = traj.analyze_hbonds()
                        hbondlists = [OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds]
                        for r1_val in total:
                            if total[r1_val][0] != 0:
                                cumul_hbonds[r1_val][0] += total[r1_val][0]
                                cumul_hbonds[r1_val][6] += total[r1_val][1]
                                cumul_hbonds[r1_val][7] += total[r1_val][2]
                                for i in range(5):
                                    for inter in hbondlists[i][r1_val]:
                                        if inter not in keylist[i]:
                                            keylist[i].append(inter)
                                        block_incr_dict(hbondlists[i][r1_val][inter], cumul_hbonds[r1_val][i+1], inter)
                                    keylist[i].sort()
            output = open(pfile, 'wb') 
            pickle.dump(cumul_hbonds, output, -1)
            pickle.dump(keylist, output, -1)
            output.close()

            all_key_lists.append(keylist)
            
        # cumul_lists: OAI, OAD, KCX, WAT, OG
        cumul_lists = [[] for x in range(5)]
        
        for rpoint in cumul_hbonds:
            total = cumul_hbonds[rpoint][0]
            #will want an indicator to know if total was zero...
            total = float(total)
            for i in range(5):
                add_to_list(cumul_lists[i], keylist[i], rpoint, cumul_hbonds, total,i+1)

        sorted_cumul_lists = []
        for i in range(5):
            sorted_cumul_lists.append(sorted(cumul_lists[i]))
        all_cumul_lists.append(sorted_cumul_lists)

    who = ['OAI', 'OAD', 'KCX', 'WAT', 'OG']

    for x in range(5):
        print who[x] + " Interacts with:"
        for y in range(2):
            print all_key_lists[y][x]
            plot_allkeys(all_cumul_lists[y][x], all_key_lists[y][x], -1, 1.5, "TPS: " + moltypes[y] + " : " + who[x])


