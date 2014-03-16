#umbsamp_ensemble_analysis with numpy
#Umbrella Sampling Ensemble Analysis

import MDAnalysis.analysis.hbonds as hydbond
import MDAnalysis
import math
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np
import os
import pickle
import re


plt.ion()

global debug
debug = False

global testing
testing = True

def print_hbonds(b, inter, r1, r2):
    if debug:
        print b
        print "Is classified to be interacting with " + inter
        print "At r1 = ", r1, "And r2 = ", r2
        
#To investigate:
#   What is the hydroxyethyl group interacting with?
#       KCX?, beta5-beta6 loop?
#   What is balancing the oxygen in the intermediate?
#   What attacks the ester? OH or H2O?
#   Map of proton moving from water
#   
#Can I investigate?:
#   Tautomerization

def incr_dict(dictionary, key):
    if key in dictionary:
        dictionary[key] += 1
    else:
        dictionary[key] = 1

def block_incr_dict(amount, d, k):
    if k in d:
        d[k] += amount
    else:
        d[k] = amount

regex = '[A-Z]+[0-9]+(?=:)'

def getKey(item):
    return item[0], item[1]

class UmbEnsemble:
    def __init__(self, universe, is_dori, moltype, foldername, r1, r2, r1trace, r2trace, r2extra):
        """ if is_dori:
                moltype = SDR or SSD
            else:
                moltype = SIM or SMI"""
        self.universe = universe
        self.small_trajs = []
        self.is_dori = is_dori
        self.moltype = moltype
        self.name = foldername
        self.r1 = r1
        self.r2 = r2
        self.r1trace = r1trace
        self.r2trace = r2trace
        self.r2extra = r2extra
        self.attached = self.analyze_protons()
    def analyze_protons(self):
        #r2extra is list of lists:
        #   [dist(H1 - OH1), dist(H1 - OH2), dist(H1-watOH2), dist(H2 - OH1), dist(H2 - OH2), dist(H2-watOH2)]
        prot_kcx = []
        r_cov = 1.31
        for h11, h12, h1w, h21, h22, h2w in self.r2extra:
            s = []
            attached = False
            if h11 <= r_cov:
                attached = True
                s.append('(H1-OH1)')
            elif h12 <= r_cov:
                attached = True
                s.append('(H1-OH2)')
            elif h21 <= r_cov:
                attached = True
                s.append('(H2-OH1)')
            elif h22 <= r_cov:
                attached = True
                s.append('(H2-OH2)')
            toappend = [attached] + [x for x in s]
            prot_kcx.append(toappend)
         self.attached = prot_kcx
         return prot_kcx
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
        r1lower = int((self.r1 - 2)*100)
        r1higher = int((self.r1 + 2)*100)
        r2lower = int((self.r2 - 2)*100)
        r2higher = int((self.r2 + 2)*100)

        matrix = [] 
        for y in range(r2lower, r2higher):
            for x in range(r1lower, r1higher):
                matrix.append((x,y))

        total = {key:0 for key in matrix}
        OAIhbonds = {key:{} for key in matrix}
        OADhbonds = {key:{} for key in matrix}
        KCX1hbonds = {key:{} for key in matrix}
        KCX2hbonds = {key:{} for key in matrix}
        WAThbonds = {key:{} for key in matrix}

        frame_num = 0
        #Fill dictionaries; 0 if bond is not present; 1 if it is
        for prot_kcx, r1, r2, frame in zip(self.attached, self.r1trace, self.r2trace, h_bond_results):
            r1 = int(round(r1, 2)*100)
            r2 = int(round(r2, 2)*100)
            total[(r1,r2)] += 1
            frame_num += 1
            for hbond in frame:
                m = re.search(regex, hbond[2])
                donorstr = m.group(0)
                m = re.search(regex, hbond[3])
                acceptorstr = m.group(0)
                #Go through donors and add acceptor to correct dictionary
                if hbond[2] == OAIdonor:
                    print_hbonds(hbond, "OAIdonor", r1/100., r2/100.)
                    incr_dict(OAIhbonds[(r1,r2)], acceptorstr)
                if (hbond[2] == WATdonor1) or (hbond[2] == WATdonor2):
                    print_hbonds(hbond, "Watdonor", r1/100., r2/100.)
                    incr_dict(WAThbonds[(r1,r2)], acceptorstr)
                    
                #Go through acceptors and add donor to correct dictionary
                if hbond[3] == OAIacceptor:
                    print_hbonds(hbond, "OAIacceptor", r1/100., r2/100.)
                    incr_dict(OAIhbonds[(r1,r2)], donorstr)
                if hbond[3] == OADacceptor:
                    print_hbonds(hbond, "OADacceptor", r1/100., r2/100.)
                    incr_dict(OADhbonds[(r1,r2)], donorstr)
                if hbond[3] == KCX1acceptor:
                    print_hbonds(hbond, "KCX1acceptor", r1/100., r2/100.)
                    incr_dict(KCX1hbonds[(r1,r2)], donorstr)
                if hbond[3] == KCX2acceptor:
                    print_hbonds(hbond, "KCX2acceptor", r1/100., r2/100.)
                    incr_dict(KCX2hbonds[(r1,r2)], donorstr)
                if hbond[3] == WATacceptor:
                    print_hbonds(hbond, "Watacceptor", r1/100., r2/100.)
                    incr_dict(WAThbonds[(r1,r2)], donorstr)       
        return total, OAIhbonds, OADhbonds, KCX1hbonds, KCX2hbonds, WAThbonds
    
    def plot(self):
        plt.figure()
        plt.plot(self.r1trace, self.r2trace, 'go')
        plt.draw()
    

def get_dist(atom1, atom2):
    return math.sqrt(sum((atom1[i] - atom2[i])**2 for i in range(3)))

def analyze_umbsamp(dcd_filepath, psf_filepath, is_dori, moltype, r1, r2):
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
    r1trace = []
    r2trace = []
    r2extra = []
    #r2 is list of lists:
    #   [dist(H1 - OH1), dist(H1 - OH2), dist(H1-watOH2), dist(H2 - OH1), dist(H2 - OH2), dist(H2-watOH2)]
    rxnatoms = universe.trajectory.timeseries(aoi)
    for t in range(rxnatoms.shape[1]):
        m = get_dist(rxnatoms[4][t], rxnatoms[1][t])
        b = get_dist(rxnatoms[0][t], rxnatoms[1][t])
        r1trace.append(b - m)
        h1 = get_dist(rxnatoms[5][t], rxnatoms[2][t])
        h2 = get_dist(rxnatoms[5][t], rxnatoms[3][t])
        h3 = get_dist(rxnatoms[5][t], rxnatoms[4][t])
        h4 = get_dist(rxnatoms[6][t], rxnatoms[2][t])
        h5 = get_dist(rxnatoms[6][t], rxnatoms[3][t])
        h6 = get_dist(rxnatoms[6][t], rxnatoms[4][t])
        if is_dori:
            r2trace.append(-h2)
        else:
            r2trace.append(-h1)
        r2extra.append([h1, h2, h3, h4, h5, h6])
    ensemble = UmbEnsemble(universe, is_dori, moltype, dcd_filepath, r1, r2, r1trace, r2trace, r2extra)
    return ensemble

def add_to_list(listtoaddto, keylist, rpoint, cumul_hbonds, total, i):
    mini_list = [rpoint[0]/10., rpoint[1]/10.]
    
    for key in keylist:
        if key in cumul_hbonds[rpoint][i]:
            if total == 0:
                prob = 0
            else:
                prob = cumul_hbonds[rpoint][i][key]/total
        else:
            prob = 0
        mini_list.append(prob)
        #TODO: check that this is correct variance
        mini_list.append(prob*(1-prob))
    listtoaddto.append(mini_list)

def plot_2Dhist (sorted_list, key_list, figname, size):
    """ intoTS and outofTS are TS basin definitions. Should be real numbers
        figname should be string to write figure name to. Should include .png"""
    annoying = zip(*(sorted_list))
    r1_vals, r2_vals, annoying = annoying[0], annoying[1], annoying[2:]
    legend_handles = []
    
    new_r1_vals = [[] for i in range(len(r1_vals)/size)]
    new_r2_vals = [[] for i in range(len(r2_vals)/size)]
    l_of_np_r1 = []
    l_of_np_r2 = []
    for i, x in enumerate(r1_vals):
        new_r1_vals[i/size].append(x)
    for i, x in enumerate(r2_vals):
        new_r2_vals[i/size].append(x)
    for l in new_r1_vals:
        l_of_np_r1.append(np.array(l))
    for l in new_r2_vals:
        l_of_np_r2.append(np.array(l))
    r1test = np.array(l_of_np_r1)
    r2test = np.array(l_of_np_r2)

    for i,key in enumerate(key_list):
        avg, err = annoying[0], annoying[1]
        if len(annoying) > 2:
            annoying = annoying[2:]
        new_avg = [[] for i in range(len(avg)/size)]
        l = []
        for j,x in enumerate(avg):
            new_avg[j/size].append(x)
        for foo in new_avg:
            l.append(np.array(foo))
        
        avgtest = np.array(l)
        plt.figure()
        plt.pcolor(r1test, r2test, avgtest, cmap='jet', vmin=0, vmax=2)
        plt.xlabel('R1')
        plt.ylabel('R2')
        plt.title(figname + " : " + key)
        plt.colorbar()
        #plt.draw()
        plt.savefig(figname + " : " + key)
        plt.close('all')

infolist = []

psfpaths = ["/data/sguthrie/imivsdori/dori_sim/sim2/template/dori.psf", "/data/sguthrie/imivsdori/m_imi_sim/from_imi/sim1/template/mimi.psf",
            "/data/sguthrie/imivsdori/sdori_sim/from_dori/sim1/template/sdori.psf", "/data/sguthrie/imivsdori/imi_sim/sim2/template/imi.psf"]
rootpaths = ["/data/sguthrie/imivsdori/dori_sim/sim2", "/data/sguthrie/imivsdori/m_imi_sim/from_imi/sim1",
            "/data/sguthrie/imivsdori/sdori_sim/from_dori/sim1", "/data/sguthrie/imivsdori/imi_sim/sim2"]
isdoris = [True, False, True, False]
moltypes = ['SDR', 'SMI', 'SSD', 'SIM']
times = [65, 65, 65, 65]
picklefiles = ['UmbDori_interactions.pkl', 'UmbMimi_interactions.pkl', 'UmbSDori_interactions.pkl', 'UmbImi_interactions.pkl']
for i in range(4):
    tmp = [psfpaths[i], rootpaths[i], isdoris[i], moltypes[i], times[i], picklefiles[i]]
    infolist.append(tmp)

#cumul_hbonds structure:
#   overall: dictionary indexed by (r1 value*100, r2 value*100)
#       dictionary hashes to list.
#       First element is integer representing the number of simulations that saved at that (r1,r2) value
#       Second element is dictionary of OAI interactions
#       Third element is dictionary of OAD interactions
#       Fourth element is dictionary of KCX1 interactions
#       Fifth element is dictionary of KCX2 interactions
#       Sixth element is dictionary of Water interactions


logregex1 = '(-?[0-9]\.[0-9]+) (-?[0-9]\.[0-9]+)'
logregex2 = '(?<=#)(/[-[0-9\.\w]+)+'
if debug or testing:
    psfpath, rootpath, isdori, moltype, time, pfile = infolist[1]
    keylist = [[] for x in range(5)]
    matrix = []
    for y in range(-35, -5):
        for x in range(-30, 40):
            matrix.append((x,y))

    cumul_hbonds = {key:[0, {}, {}, {}, {}, {}] for key in matrix}
    logfile = open(rootpath + "/umbsamp.log", "r")
    num_sims = 0
    for line in logfile:
        m = re.search(logregex1, line)
        if m != None:
            if num_sims >= 1:
                break
            r1 = m.group(1)
            r2 = m.group(2)
            x = re.search(logregex2, line)
            dcdpath = x.group(0) + "/cap_production.dcd"
            try:
                #Check if dcd exists
                f = open(dcdpath, 'r')
                f.close()
                ensemble = analyze_umbsamp(dcdpath, psfpath, isdori, moltype, float(r1), float(r2))
                num_sims += 1
                total, OAIhbonds, OADhbonds, KCX1hbonds, KCX2hbonds, WAThbonds = ensemble.analyze_hbonds()
                hbondlists = [OAIhbonds, OADhbonds, KCX1hbonds, KCX2hbonds, WAThbonds]
                for rvals in total:
                    if total[rvals] != 0:
                        sm_r1 = round(rvals[0]/10.)
                        sm_r2 = round(rvals[1]/10.)
                        r_rvals = (sm_r1, sm_r2)
                        cumul_hbonds[r_rvals][0] += total[rvals]
                        for i in range(5):
                            for inter in hbondlists[i][rvals]:
                                if inter not in keylist[i]:
                                    keylist[i].append(inter)
                                block_incr_dict(hbondlists[i][rvals][inter], cumul_hbonds[r_rvals][i+1], inter)
            except IOError:
                pass
    logfile.close()
    # cumul_lists: OAI, OAD, KCX1, KCX2, WAT
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
    for foo in range(4):
        psfpath, rootpath, isdori, moltype, time, pfile = infolist[foo]
        print rootpath
        try:
            inp = open(pfile, 'rb')
            sorted_cumul_lists = pickle.load(inp)
            all_cumul_lists.append(sorted_cumul_lists)
            keylist = pickle.load(inp)
            all_key_lists.append(keylist)
            inp.close()
        except IOError:
            # keylists: OAI, OAD, KCX1, KCX2, WAT
            keylist = [[] for x in range(5)]
            matrix = []
            for y in range(-70, -5):
                for x in range(-30, 40):
                    matrix.append((x,y))

            cumul_hbonds = {key:[0, {}, {}, {}, {}, {}] for key in matrix}
            logfile = open(rootpath + "/umbsamp.log", "r")
            i = 0
            for line in logfile:
                m = re.search(logregex1, line)
                if m != None:
                    r1 = m.group(1)
                    r2 = m.group(2)
                    x = re.search(logregex2, line)
                    dcdpath = x.group(0) + "/cap_production.dcd"
                    try:
                        #Check if dcd exists
                        f = open(dcdpath, 'r')
                        f.close()
                        ensemble = analyze_umbsamp(dcdpath, psfpath, isdori, moltype, float(r1), float(r2))
                        i += 1
                        total, OAIhbonds, OADhbonds, KCX1hbonds, KCX2hbonds, WAThbonds = ensemble.analyze_hbonds()
                        hbondlists = [OAIhbonds, OADhbonds, KCX1hbonds, KCX2hbonds, WAThbonds]
                        for rvals in total:
                            if total[rvals] != 0:
                                sm_r1 = round(rvals[0]/10.)
                                sm_r2 = round(rvals[1]/10.)
                                r_rvals = (sm_r1, sm_r2)
                                cumul_hbonds[r_rvals][0] += total[rvals]
                                for i in range(5):
                                    for inter in hbondlists[i][rvals]:
                                        if inter not in keylist[i]:
                                            keylist[i].append(inter)
                                        block_incr_dict(hbondlists[i][rvals][inter], cumul_hbonds[r_rvals][i+1], inter)
                    except IOError:
                        pass
            logfile.close()
            # cumul_lists: OAI, OAD, KCX1, KCX2, WAT
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
            
            output = open(pfile, 'wb') 
            pickle.dump(sorted_cumul_lists, output, -1)
            pickle.dump(keylist, output, -1)
            output.close()
            
            all_cumul_lists.append(sorted_cumul_lists)
            all_key_lists.append(keylist)

    who = ['OAI', 'OAD', 'KCX1', 'KCX2', 'WAT']

    for x in range(5):
        print who[x] + " Interacts with:"
        for y in range(4):
            print all_key_lists[y][x]
            
            plot_2Dhist(all_cumul_lists[y][x], all_key_lists[y][x], moltypes[y] + " : " + who[x], 65)


