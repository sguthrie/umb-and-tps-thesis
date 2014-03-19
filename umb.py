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
testing = False

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
    dictionary[key] = 1
##    if key in dictionary:
##        dictionary[key] += 1
##    else:
##        dictionary[key] = 1

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
        self.kcx_attached, self.og_attached = self.analyze_protons()
    def analyze_protons(self):
        #r2extra is list of lists:
        #   [dist(H1 - OH1), dist(H1 - OH2), dist(H1-watOH2),
        #    dist(H2 - OH1), dist(H2 - OH2), dist(H2-watOH2),
        #    dist(H1-OG), dist(H2-OG)]
        prot_kcx = []
        prot_og = []
        r_cov = 1.31
        for h11, h12, h1w, h21, h22, h2w, h1o, h2o in self.r2extra:
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
        r1lower = int((self.r1 - 2)*100)
        r1higher = int((self.r1 + 2)*100)
        r2lower = int((self.r2 - 2)*100)
        r2higher = int((self.r2 + 2)*100)

        matrix = [] 
        for y in range(r2lower, r2higher):
            for x in range(r1lower, r1higher):
                matrix.append((x,y))

        total = {key:[0, 0, 0] for key in matrix}
        OAIhbonds = {key:{} for key in matrix}
        OADhbonds = {key:{} for key in matrix}
        KCXhbonds = {key:{} for key in matrix}
        WAThbonds = {key:{} for key in matrix}
        OGhbonds = {key:{} for key in matrix}

        r_cov = 1.31
        #Fill dictionaries; 0 if bond is not present; 1 if it is
        for tstep, foo in enumerate(zip(self.kcx_attached, self.og_attached, self.r1trace, self.r2trace, h_bond_results)):
            prot_kcx, prot_og, r1, r2, frame = foo
            r1 = int(round(r1, 2)*100)
            r2 = int(round(r2, 2)*100)
            total[(r1,r2)][0] += 1
            if prot_kcx[0]:
                total[(r1,r2)][1] += 1
            if prot_og[0]:
                total[(r1,r2)][2] += 1
            #r2extra is list of lists:
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
                    print_hbonds(hbond, "OAIdonor", r1/100., r2/100.)
                    incr_dict(OAIhbonds[(r1,r2)], acceptorstr)
                if hbond[2] == WATdonor1:
                    if self.r2extra[tstep][0] <= r_cov or self.r2extra[tstep][1] <= r_cov:
                        print_hbonds(hbond, "KCXdonor", r1/100., r2/100.)
                        incr_dict(KCXhbonds[(r1,r2)], acceptorstr)
                    if self.r2extra[tstep][2] <= r_cov:
                        print_hbonds(hbond, "WATdonor", r1/100., r2/100.)
                        incr_dict(WAThbonds[(r1,r2)], acceptorstr)
                    if self.r2extra[tstep][6] <= r_cov:
                        print_hbonds(hbond, "OGdonor", r1/100., r2/100.)
                        incr_dict(OGhbonds[(r1,r2)], acceptorstr)
                    #Should not ever not be in the case where it's flying in the ether, since it wouldn't be
                    #counted by the hbonds.py algorithm
                        
                if hbond[2] == WATdonor2:
                    if self.r2extra[tstep][3] <= r_cov or self.r2extra[tstep][4] <= r_cov:
                        print_hbonds(hbond, "KCXdonor", r1/100., r2/100.)
                        incr_dict(KCXhbonds[(r1,r2)], acceptorstr)
                    if self.r2extra[tstep][5] <= r_cov:
                        print_hbonds(hbond, "WATdonor", r1/100., r2/100.)
                        incr_dict(WAThbonds[(r1,r2)], acceptorstr)
                    if self.r2extra[tstep][7] <= r_cov:
                        print_hbonds(hbond, "OGdonor", r1/100., r2/100.)
                        incr_dict(OGhbonds[(r1,r2)], acceptorstr)
                    #Should not ever not be in the case where it's flying in the ether, since it wouldn't be
                    #counted by the hbonds.py algorithm
                        
                #Go through acceptors and add donor to correct dictionary
                if hbond[3] == OAIacceptor:
                    print_hbonds(hbond, "OAIacceptor", r1/100., r2/100.)
                    incr_dict(OAIhbonds[(r1,r2)], donorstr)
                if hbond[3] == OADacceptor:
                    print_hbonds(hbond, "OADacceptor", r1/100., r2/100.)
                    incr_dict(OADhbonds[(r1,r2)], donorstr)
                if (hbond[3] == KCX1acceptor) or (hbond[3] == KCX2acceptor):
                    print_hbonds(hbond, "KCXacceptor", r1/100., r2/100.)
                    incr_dict(KCXhbonds[(r1,r2)], donorstr)
                if hbond[3] == WATacceptor:
                    print_hbonds(hbond, "WATacceptor", r1/100., r2/100.)
                    incr_dict(WAThbonds[(r1,r2)], donorstr)
                if hbond[3] == OGacceptor:
                    print_hbonds(hbond, "OGacceptor", r1/100., r2/100.)
                    incr_dict(OGhbonds[(r1,r2)], donorstr)
        return total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds
    
    def plot(self):
        plt.figure()
        plt.plot(self.r1trace, self.r2trace, 'go')
        plt.draw()
    
def make_key_data(size, keyname, c_hbonds, inter):
    '''Calculates probability and wilson error bars at each point. Returns np masked arrays '''
    z = 1.96
    zsq = 3.84
    probvals = np.ma.masked_array(np.zeros((len(c_hbonds)/size, size)), mask=np.zeros((len(c_hbonds)/size, size)))
    upper_err = np.ma.masked_array(np.zeros((len(c_hbonds)/size, size)), mask=np.zeros((len(c_hbonds)/size, size)))
    lower_err = np.ma.masked_array(np.zeros((len(c_hbonds)/size, size)), mask=np.zeros((len(c_hbonds)/size, size)))
##    probvals = np.array(np.zeros((len(c_hbonds)/size, size)))
##    upper_err = np.array(np.zeros((len(c_hbonds)/size, size)))
##    lower_err = np.array(np.zeros((len(c_hbonds)/size, size)))
    rpointlist = c_hbonds.keys()
    rpointlist.sort()
    for i, rpoint in enumerate(rpointlist):
        total = float(c_hbonds[rpoint][0])
        if total == 0:
            #pass
            probvals[i/size][i % size] = np.ma.masked
            upper_err[i/size][i % size] = np.ma.masked
            lower_err[i/size][i % size] = np.ma.masked
        else:
            if keyname in c_hbonds[rpoint][inter+1]:
                prob = c_hbonds[rpoint][inter+1][keyname]/total
            else:
                prob = 0
            wscorepm = z*math.sqrt(prob*(1-prob)/total + zsq/(4*total**2))
            wmult = 1/(1+zsq/total)
            probvals[i/size][i % size] = prob
            upper_err[i/size][i % size] = wmult*(prob + zsq/(2*total) + wscorepm)
            lower_err[i/size][i % size] = wmult*(prob + zsq/(2*total) - wscorepm)
    return probvals, upper_err, lower_err

class Umbdata:
    def __init__(self, moltype, cumul_hbonds, keylist, size):
        self.moltype = moltype
        self.size = size
        self.keylist = keylist
        self.r1points = np.zeros((len(cumul_hbonds)/size, size))
        self.r2points = np.zeros((len(cumul_hbonds)/size, size))
        self.nsims = np.zeros((len(cumul_hbonds)/size, size))
        self.kcx_attached = np.zeros((len(cumul_hbonds)/size, size))
        self.og_attached = np.zeros((len(cumul_hbonds)/size, size))
        rpointlist = cumul_hbonds.keys()
        rpointlist.sort()
        for i, rpoint in enumerate(rpointlist):
            self.r1points[i/size][i % size] = rpoint[0]/10.
            self.r2points[i/size][i % size] = rpoint[1]/10.
            self.nsims[i/size][i % size] = cumul_hbonds[rpoint][0]
            if cumul_hbonds[rpoint][0] != 0:
                self.kcx_attached[i/size][i % size] = cumul_hbonds[rpoint][6]/float(cumul_hbonds[rpoint][0])
                self.og_attached[i/size][i % size] = cumul_hbonds[rpoint][7]/float(cumul_hbonds[rpoint][0])
            
        self.OAI = {}
        self.OAD = {}
        self.KCX = {}
        self.WAT = {}
        self.OG = {}
        self.interdicts = [self.OAI, self.OAD, self.KCX, self.WAT, self.OG]

        for i in range(5):
            for keyname in self.keylist[i]:
                probvals, upper_err, lower_err = make_key_data(self.size, keyname, cumul_hbonds, i)
                self.interdicts[i][keyname] = [probvals, upper_err, lower_err]

    def get_plot_data(self, inter_i, keyname):
        if keyname in self.interdicts[inter_i]:
            return self.interdicts[inter_i][keyname]
        else:
            return [None, None, None]


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
    #   [dist(H1 - OH1), dist(H1 - OH2), dist(H1-watOH2), dist(H2 - OH1), dist(H2 - OH2), dist(H2-watOH2), dist(H1-OG), dist(H2-OG)]
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
        h7 = get_dist(rxnatoms[5][t], rxnatoms[0][t])
        h8 = get_dist(rxnatoms[6][t], rxnatoms[0][t])
        if is_dori:
            r2trace.append(-h2)
        else:
            r2trace.append(-h1)
        r2extra.append([h1, h2, h3, h4, h5, h6, h7, h8])
    ensemble = UmbEnsemble(universe, is_dori, moltype, dcd_filepath, r1, r2, r1trace, r2trace, r2extra)
    return ensemble

def plot_2Dhist(bigstruct, inter_i, keyname, figname):
    plt.figure(figsize=(12,12))
    pnum = [[1,3,5],[2,4,6], [7,9,11], [8,10,12]]
    for i in range(4):
        plt.subplot(6, 2, pnum[i][0])
        [avg, uperr, lowerr] = bigstruct[i].get_plot_data(inter_i, keyname)
        if avg != None:
            plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, uperr, vmin=0, vmax=1)
            plt.xlabel('R1')
            plt.ylabel('R2')
            plt.title("UpperError for " + bigstruct[i].moltype + " : " + figname)
            plt.colorbar()
            plt.subplot(6,2, pnum[i][1])
            plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, avg, vmin=0, vmax=1)
            plt.xlabel('R1')
            plt.ylabel('R2')
            plt.title(bigstruct[i].moltype + " : " + figname)
            plt.colorbar()
            plt.subplot(6,2, pnum[i][2])
            plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, lowerr, vmin=0, vmax=1)
            plt.xlabel('R1')
            plt.ylabel('R2')
            plt.title("LowerError for " + bigstruct[i].moltype + " : " + figname)
            plt.colorbar()
        else:
            plt.title(bigstruct[i].moltype + " had no " + figname + " interactions")
    plt.draw()
    plt.savefig("Umb" + figname + ".png")
    plt.close()


def plot_2Dhist_noerror(bigstruct, inter_i, keyname, figname):
    plt.figure(figsize=(12,10))
    for i in range(4):
        [avg, uperr, lowerr] = bigstruct[i].get_plot_data(inter_i, keyname)
        if avg != None:
            plt.subplot(2, 2, 1+i)
            plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, avg, vmin=0, vmax=1)
            plt.xlabel('R1')
            plt.ylabel('R2')
            plt.title(bigstruct[i].moltype + " : " + figname)
            plt.colorbar()
        else:
            plt.text(0, 0, bigstruct[i].moltype + " had no " + figname + " interactions")
    plt.draw()

infolist = []
psfpaths = ["/data/sguthrie/imivsdori/dori_sim/sim2/template/dori.psf", "/data/sguthrie/imivsdori/m_imi_sim/from_imi/sim1/template/mimi.psf",
            "/data/sguthrie/imivsdori/sdori_sim/from_dori/sim1/template/sdori.psf", "/data/sguthrie/imivsdori/imi_sim/sim2/template/imi.psf"]
rootpaths = ["/data/sguthrie/imivsdori/dori_sim/sim2", "/data/sguthrie/imivsdori/m_imi_sim/from_imi/sim1",
            "/data/sguthrie/imivsdori/sdori_sim/from_dori/sim1", "/data/sguthrie/imivsdori/imi_sim/sim2"]
isdoris = [True, False, True, False]
moltypes = ['SDR', 'SMI', 'SSD', 'SIM']
times = [65, 65, 65, 65]
picklefiles = ['UmbDori_interactions2.pkl', 'UmbMimi_interactions2.pkl', 'UmbSDori_interactions2.pkl', 'UmbImi_interactions2.pkl']
for i in range(4):
    tmp = [psfpaths[i], rootpaths[i], isdoris[i], moltypes[i], times[i], picklefiles[i]]
    infolist.append(tmp)

#cumul_hbonds structure:
#   overall: dictionary indexed by (r1 value*100, r2 value*100)
#       dictionary hashes to list.
#       First element is integer representing the number of simulations that saved at that (r1,r2) value
#       Second element is dictionary of OAI interactions
#       Third element is dictionary of OAD interactions
#       Fourth element is dictionary of KCX interactions
#       Fifth element is dictionary of Water interactions
#       Sixth element is dictionary of OG interactions
#       Seventh element is integer representing the number of simulations that had a proton on KCX
#       Eigth element is integer representing the number of simulations that had a proton on OG
##        try:
##            inp = open(pfile, 'rb')
##            bigst = pickle.load(inp)
##            keylist = pickle.load(inp)
##            inp.close()
##        except IOError:
            # keylists: OAI, OAD, KCX, WAT, OG

logregex1 = '(-?[0-9]\.[0-9]+) (-?[0-9]\.[0-9]+)'
logregex2 = '(?<=#)(/[-[0-9\.\w]+)+'
if debug or testing:
    psfpath, rootpath, isdori, moltype, time, pfile = infolist[1]
    keylist = [[] for x in range(5)]
    matrix = []
    for y in range(-30, -5):
        for x in range(-30, 30):
            matrix.append((x,y))

    cumul_hbonds = {key:[0, {}, {}, {}, {}, {}, 0, 0] for key in matrix}
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
                total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds = ensemble.analyze_hbonds()
                hbondlists = [OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds]
                for rvals in total:
                    if total[rvals][0] != 0:
                        sm_r1 = round(rvals[0]/10.)
                        sm_r2 = round(rvals[1]/10.)
                        r_rvals = (sm_r1, sm_r2)
                        cumul_hbonds[r_rvals][0] += total[rvals][0]
                        cumul_hbonds[r_rvals][6] += total[rvals][1]
                        cumul_hbonds[r_rvals][7] += total[rvals][2]
                        for i in range(5):
                            for inter in hbondlists[i][rvals]:
                                if inter not in keylist[i]:
                                    keylist[i].append(inter)
                                block_incr_dict(hbondlists[i][rvals][inter], cumul_hbonds[r_rvals][i+1], inter)
                            keylist[i].sort()
            except IOError:
                pass
    logfile.close()
    
else:
    try:
        inp = open('UmbInteractions.pkl', 'rb')
        bigstruct = pickle.load(inp)
        all_key_lists = pickle.load(inp)
        inp.close()
    except IOError:
        bigstruct = []
        all_key_lists = []
        for foo in range(4):
            psfpath, rootpath, isdori, moltype, time, pfile = infolist[foo]
            print rootpath
            keylist = [[] for x in range(5)]
            matrix = []
            for y in range(-30, -5):
                for x in range(-30, 30):
                    matrix.append((x,y))
            cumul_hbonds = {key:[0, {}, {}, {}, {}, {}, 0, 0] for key in matrix}
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
                        total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds = ensemble.analyze_hbonds()
                        hbondlists = [OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds]
                        for rvals in total:
                            if total[rvals][0] != 0:
                                sm_r1 = round(rvals[0]/10.)
                                sm_r2 = round(rvals[1]/10.)
                                r_rvals = (sm_r1, sm_r2)
                                cumul_hbonds[r_rvals][0] += total[rvals][0]
                                cumul_hbonds[r_rvals][6] += total[rvals][1]
                                cumul_hbonds[r_rvals][7] += total[rvals][2]
                                for i in range(5):
                                    for inter in hbondlists[i][rvals]:
                                        if inter not in keylist[i]:
                                            keylist[i].append(inter)
                                        block_incr_dict(hbondlists[i][rvals][inter], cumul_hbonds[r_rvals][i+1], inter)
                                    keylist[i].sort()
                    except IOError:
                        pass
            logfile.close()
            output = open(pfile, 'wb')
            pickle.dump(cumul_hbonds, output, -1)
            pickle.dump(keylist, output, -1)
            output.close()
            
##            all_key_lists.append(keylist)
##            data = Umbdata(moltype, cumul_hbonds, keylist, time)
##            bigstruct.append(data)
##
##        output = open('UmbInteractions.pkl', 'wb') 
##        pickle.dump(bigstruct, output, -1)
##        pickle.dump(all_key_lists, output, -1)
##        output.close()

    who = ['OAI', 'OAD', 'KCX', 'WAT', 'OG']
    keysets = []
    for x in range(5):
        #print who[x], "interacts with: "
        foo = set()
        for y in range(4):
            foo |= set(all_key_lists[y][x])
        keysets.append(foo)
        #print foo
        
    for x in range(5):
        for keyname in keysets[x]:
            plot_2Dhist(bigstruct, x, keyname, who[x] + " : " + keyname)
#plot_2Dhist(bigstruct, 0, 'SER128', who[0] + " : SER128")

##
##def plot_foo(bigstruct, figname):
##    plt.figure(figsize=(12,10), dpi=100)
##    for i in range(4):
##        plt.subplot(2, 2, 1+i)
##        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, bigstruct[i].kcx_attached, vmin=0, vmax=1)
##        plt.xlabel('R1')
##        plt.ylabel('R2')
##        plt.title(bigstruct[i].moltype + " : " + figname)
##        plt.colorbar()
##
##    plt.draw()

#plot_foo(bigstruct, "testing kcx interactions")



