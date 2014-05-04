#TPS plotted on 2 reaction coordinates
#Umbrella Sampling Ensemble Analysis

import MDAnalysis.analysis.hbonds as hydbond
import MDAnalysis.analysis.align as ali
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

class TrajEnsemble:
    def __init__(self, universe, is_dori, moltype, foldername, r1trace, r2trace, r2traces):
        """ if is_dori:
                moltype = SDR or SSD
            else:
                moltype = SIM or SMI"""
        self.universe = universe
        self.is_dori = is_dori
        self.moltype = moltype
        self.name = foldername
        self.r1min = min(r1trace)
        self.r1max = max(r1trace)
        self.r2min = min(r2trace)
        self.r2max = max(r2trace)
        self.r1trace = r1trace
        self.r2trace = r2trace
        self.r2extra = r2traces
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
            examinestr = "(atom A 81 OAI) or (atom A 81 OAD) or (atom A 84 OH1) or (atom A 84 OH1) " + \
                         "or (atom A 84 OH2) or (atom W 277 OH2) or (atom A 81 OG) " + \
                         "or (atom A 81 NAO) or (atom A 81 OAG) or (atom A 81 OAF) or (atom A 81 NAC) or (atom A 81 NAP)"
            #Not sure that these tail atoms can strictly be said to hydrogen bond
            #OG is already a donor and acceptor
            new_donors = ['OAI', 'NAO', 'NAC', 'OH1', 'OH2', 'NAP'] 
            new_acceptors = ['OAD', 'OAI', 'OAH', 'OAE', 'OH2', 'OH1']
            OAIacceptor = self.moltype + '81:OAI'
            OAIdonor = self.moltype + '81:HOI'
            OADacceptor = self.moltype + '81:OAD'
            taildonors = [self.moltype + '81:HN1', self.moltype + '81:HN2', self.moltype + '81:HN3', self.moltype + '81:HN4', self.moltype + '81:HNP']
            tailacceptors = [self.moltype + '81:OAG', self.moltype + '81:OAF']
        else:
            examinestr = "(atom A 81 O62) or (atom A 81 O7) or (atom A 84 OH1) or (atom A 84 OH1) " + \
                         "or (atom A 84 OH2) or (atom W 277 OH2) or (atom A 81 OG) " + \
                         "or (atom A 81 N24) or (atom A 81 N26)"
            new_donors = ['O62', 'N24', 'N26', 'OH1', 'OH2']
            new_acceptors = ['O7', 'O62', 'O31', 'O32', 'OH2', 'OH1']
            OAIacceptor = self.moltype + '81:O62'
            OAIdonor = self.moltype + '81:HO6'
            OADacceptor = self.moltype + '81:O7'
            taildonors = [self.moltype + '81:HN24', self.moltype + '81:HN61', self.moltype + '81:HN62']
            tailacceptors = []
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
                                            acceptors=new_acceptors) #, angle=150.0
        hana.run()
        h_bond_results = hana.timeseries
        r1lower = int((self.r1min - 1)*100)
        r1higher = int((self.r1max + 1)*100)
        r2lower = int((self.r2min - 1)*100)
        r2higher = int((self.r2max + 1)*100)

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
        tailhbonds = {key:{} for key in matrix}

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
                if hbond[2] in taildonors:
                    print_hbonds(hbond, "TailDonor", r1/100., r2/100.)
                    incr_dict(tailhbonds[(r1,r2)], acceptorstr)
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
                if hbond[3] in tailacceptors:
                    print_hbonds(hbond, "TailAcceptor", r1/100., r2/100.)
                    incr_dict(tailhbonds[(r1,r2)], donorstr)
        return total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds, tailhbonds
    
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
    rpointlist = c_hbonds.keys()
    rpointlist.sort()
    for i, rpoint in enumerate(rpointlist):
        total = float(c_hbonds[rpoint][0])
        if total == 0:
            probvals[i/size][i % size] = np.ma.masked
            upper_err[i/size][i % size] = np.ma.masked
            lower_err[i/size][i % size] = np.ma.masked
        else:
            if keyname in c_hbonds[rpoint][inter+3]:
                prob = c_hbonds[rpoint][inter+3][keyname]/total
            else:
                prob = 0
            wscorepm = z*math.sqrt(prob*(1-prob)/total + zsq/(4*total**2))
            wmult = 1/(1+zsq/total)
            probvals[i/size][i % size] = prob
            upper_err[i/size][i % size] = wmult*(prob + zsq/(2*total) + wscorepm)
            lower_err[i/size][i % size] = wmult*(prob + zsq/(2*total) - wscorepm)
    return probvals, upper_err, lower_err

#cumul_hbonds structure:
#   overall: dictionary indexed by (r1 value*100, r2 value*100)
#       dictionary hashes to list.
#       0: integer - number of simulations that saved at that (r1,r2) value
#       1: integer - number of simulations that had a proton on KCX
#       2: integer - number of simulations that had a proton on OG
#       3: dictionary - OAI interactions
#       4: dictionary - OAD interactions
#       5: dictionary - KCX interactions
#       6: dictionary - Water interactions
#       7: dictionary - OG interactions
#       8: dictionary - Tail interactions  


# keylists: OAI, OAD, KCX, WAT, OG, Tail

class Trajdata:
    def __init__(self, moltype, cumul_hbonds, keylist, size):
        z = 1.96
        zsq = 3.84
        self.moltype = moltype
        self.size = size
        self.keylist = keylist
        self.r1points = np.zeros((len(cumul_hbonds)/size, size))
        self.r2points = np.zeros((len(cumul_hbonds)/size, size))
        self.nsims = np.zeros((len(cumul_hbonds)/size, size))
        self.kcx_attached = np.ma.masked_array(np.zeros((len(cumul_hbonds)/size, size)), mask=np.zeros((len(cumul_hbonds)/size, size)))
        self.kcx_upper_err = np.ma.masked_array(np.zeros((len(cumul_hbonds)/size, size)), mask=np.zeros((len(cumul_hbonds)/size, size)))
        self.kcx_lower_err = np.ma.masked_array(np.zeros((len(cumul_hbonds)/size, size)), mask=np.zeros((len(cumul_hbonds)/size, size)))
        self.og_attached = np.ma.masked_array(np.zeros((len(cumul_hbonds)/size, size)), mask=np.zeros((len(cumul_hbonds)/size, size)))
        self.og_upper_err = np.ma.masked_array(np.zeros((len(cumul_hbonds)/size, size)), mask=np.zeros((len(cumul_hbonds)/size, size)))
        self.og_lower_err = np.ma.masked_array(np.zeros((len(cumul_hbonds)/size, size)), mask=np.zeros((len(cumul_hbonds)/size, size)))
        rpointlist = cumul_hbonds.keys()
        rpointlist.sort()
        for i, rpoint in enumerate(rpointlist):
            self.r1points[i/size][i % size] = rpoint[0]/10.
            self.r2points[i/size][i % size] = rpoint[1]/10.
            total = cumul_hbonds[rpoint][0]
            self.nsims[i/size][i % size] = total
            if cumul_hbonds[rpoint][0] == 0:
                self.kcx_attached[i/size][i % size] = np.ma.masked
                self.og_attached[i/size][i % size] = np.ma.masked
                self.kcx_upper_err[i/size][i % size] = np.ma.masked
                self.kcx_lower_err[i/size][i % size] = np.ma.masked
                self.og_upper_err[i/size][i % size] = np.ma.masked
                self.og_lower_err[i/size][i % size] = np.ma.masked
            else:
                total = float(total)
                probkcx = cumul_hbonds[rpoint][1]/total
                probog = cumul_hbonds[rpoint][2]/total
                kcx_wscorepm = z*math.sqrt(probkcx*(1-probkcx)/total + zsq/(4*total**2))
                og_wscorepm = z*math.sqrt(probog*(1-probog)/total + zsq/(4*total**2))
                wmult = 1/(1+zsq/total)
            
                self.kcx_attached[i/size][i % size] = probkcx
                self.og_attached[i/size][i % size] = probog
                self.kcx_upper_err[i/size][i % size] = wmult*(probkcx + zsq/(2*total) + kcx_wscorepm)
                self.kcx_lower_err[i/size][i % size] = wmult*(probkcx + zsq/(2*total) - kcx_wscorepm)
                self.og_upper_err[i/size][i % size] = wmult*(probog + zsq/(2*total) + og_wscorepm)
                self.og_lower_err[i/size][i % size] = wmult*(probog + zsq/(2*total) - og_wscorepm)

        self.OAI = {}
        self.OAD = {}
        self.KCX = {}
        self.WAT = {}
        self.OG = {}
        self.tail = {}
        self.interdicts = [self.OAI, self.OAD, self.KCX, self.WAT, self.OG, self.tail]

        for i in range(6):
            for keyname in self.keylist[i]:
                probvals, upper_err, lower_err = make_key_data(self.size, keyname, cumul_hbonds, i)
                self.interdicts[i][keyname] = [probvals, upper_err, lower_err]

    def get_plot_data(self, inter_i, keyname):
        if keyname in self.interdicts[inter_i]:
            return self.interdicts[inter_i][keyname]
        elif keyname in ['SDR81', 'SMI81', 'SSD81', 'SIM81']:
            return self.interdicts[inter_i][self.moltype + '81']
        else:
            return [None, None, None]
    def get_kcx_data(self):
        return self.kcx_attached, self.kcx_upper_err, self.kcx_lower_err
    def get_og_data(self):
        return self.og_attached, self.og_upper_err, self.og_lower_err


def get_dist(atom1, atom2):
    return math.sqrt(sum((atom1[i] - atom2[i])**2 for i in range(3)))

def analyze_big_traj(dcd_filepath, psf_filepath, is_dori, moltype):
    universe = MDAnalysis.Universe(psf_filepath, dcd_filepath)
##    smallref_u = MDAnalysis.Universe(refstruc_filepath)
##    bigref_u = MDAnalysis.Universe(bigref_filepath)
    
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
    #r2extra is list of lists:
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
        r2extra.append([h1, h2, h3, h4, h5, h6, h7, h8])
##    smtrajdata = [0 for i in range(10)]
##    for i in range(10):
##        if r1trace[i*501] > 0:
##            smtrajdata[i] = 1
    for i in range(10):
        #Looking for difference in water-proton distance: 2, 5
        #Saving proton-kcx difference: index 0, 1, 3, or 4
        smalltraj = r1trace[i*501:(i+1)*501]
        for indx, j in enumerate(smalltraj):
            if j > 0 and j < 0.1:
                tstep = indx
        #print smtrajdata[i], i*501 + tstep
        #print r2extra[i*501 + tstep][0], r2extra[i*501 + tstep][1], r2extra[i*501 + tstep][3], r2extra[i*501 + tstep][4]
        ## So, unfortunately, the proton moving isn't always the proton that ends up on the OG at the end...
        ## I want the one that interacts with KCX so comment out #if r2extra[(i+1)*501-smtrajdata[i]*501][2] > 1.2:
        if r2extra[i*501 + tstep][0] < 2.5:
            smextra = [-h1 for h1, h2, h3, h4, h5, h6, h7, h8 in r2extra[i*501:(i+1)*501]]
            r2trace.extend(smextra)
        elif r2extra[i*501 + tstep][1] < 2.5:
            smextra = [-h2 for h1, h2, h3, h4, h5, h6, h7, h8 in r2extra[i*501:(i+1)*501]]
            r2trace.extend(smextra)
        elif r2extra[i*501 + tstep][3] < 2.5:
            smextra = [-h4 for h1, h2, h3, h4, h5, h6, h7, h8 in r2extra[i*501:(i+1)*501]]
            r2trace.extend(smextra)
        elif r2extra[i*501 + tstep][4] < 2.5:
            smextra = [-h5 for h1, h2, h3, h4, h5, h6, h7, h8 in r2extra[i*501:(i+1)*501]]
            r2trace.extend(smextra)
        else:
            #print r2extra[(i+1)*501-smtrajdata[i]*501][2], r2extra[(i+1)*501-smtrajdata[i]*501][5]
            print r2extra[i*501 + tstep][0], r2extra[i*501 + tstep][1], r2extra[i*501 + tstep][3], r2extra[i*501 + tstep][4]
            print "Uh oh"

    
    #print len(r2trace)
    ensemble = TrajEnsemble(universe, is_dori, moltype, dcd_filepath, r1trace, r2trace, r2extra)
##    filename1 = 'Smallr1_' + str(r1) + '_r2_' + str(r2) + '_' + moltype
##    filename2 = 'Bigr1_' + str(r1) + '_r2_' + str(r2) + '_' + moltype
##    selection = 'name CA or name O or name N or name C'
##    ali.rms_fit_trj(universe, smallref_u, select=selection, rmsdfile=filename1)
##    ali.rms_fit_trj(universe, bigref_u, select=selection, rmsdfile=filename2)
    return ensemble

def smooth_avg(avg):
    smoothed = np.ma.copy(avg)
    for i, small in enumerate(avg):
        maskarr = np.ma.getmask(small)
        data = np.ma.getdata(small)
        for j in range(len(np.ma.getdata(small))):
            if not maskarr[j]:
                tmp = []
                for up in range(i-1, i+2):
                    for across in range(j-1, j+2):
                        if up != i or across != j:
                            if not np.ma.getmask(avg[up][across]):
                                tmp.append(avg[up][across])
                smoothed[i][j] = sum(tmp)/float(len(tmp))
                if np.ma.getmask(smoothed[i][j]):
                    print smoothed[i][j]
    return smoothed
    
            

def plot_2Dhist(bigstruct, inter_i, keyname, figname, ploterr=True):
    if ploterr:
        d = (22, 22)
    else:
        d = (14,8)
    plt.figure(figsize=d)
    pnum = [[1,3,5],[2,4,6], [7,9,11], [8,10,12]]
    for i in range(4):
        [avg, uperr, lowerr] = bigstruct[i].get_plot_data(inter_i, keyname)
        if ploterr:
            plt.subplot(6, 2, pnum[i][0])
            if avg != None:
                plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, uperr, vmin=0, vmax=1)
                plt.title("UpperError for " + bigstruct[i].moltype + " : " + figname)
                plt.colorbar()
                plt.subplot(6,2, pnum[i][1])
                plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, avg, vmin=0, vmax=1)
                plt.title(bigstruct[i].moltype + " : " + figname)
                plt.colorbar()
                plt.subplot(6,2, pnum[i][2])
                plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, lowerr, vmin=0, vmax=1)
                plt.title("LowerError for " + bigstruct[i].moltype + " : " + figname)
                plt.colorbar()
            else:
                plt.title(bigstruct[i].moltype + " had no " + figname + " interactions")
        else:
            plt.subplot(2,2,1+i)
            if avg != None:
                #avg = smooth_avg(avg)
                plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, avg, vmin=0, vmax=1)
                plt.xlabel('R1')
                plt.ylabel('R2')
                plt.title(bigstruct[i].moltype + " : " + figname)
                plt.colorbar()
            else:
                plt.title(bigstruct[i].moltype + " had no " + figname + " interactions")
    plt.draw()
    if ploterr:
        plt.savefig("TPSpmf" + figname + "WithError.png")
    else:
        plt.savefig("TPSpmf" + figname + ".png")
    plt.close()

def plot_prot_movement(bigstruct):
    fsize = (14,8)
    plt.figure(figsize=fsize)
    for i in range(4):
        plt.subplot(2, 2, 1+i)
        kcxavg, kcxuperr, kcxlowerr = bigstruct[i].get_kcx_data()
        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, kcxavg, vmin=0, vmax=1)
        plt.title(bigstruct[i].moltype + " : KCX protonation probability")
        plt.colorbar()
    plt.draw()
    plt.savefig("TPSpmfKCXprot.png")

    plt.figure(figsize=fsize)
    for i in range(4):
        plt.subplot(2, 2, 1+i)
        ogavg, oguperr, oglowerr = bigstruct[i].get_og_data()
        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, ogavg, vmin=0, vmax=1)
        plt.title(bigstruct[i].moltype + " : OG protonation probability")
        plt.colorbar()
        
    plt.draw()
    plt.savefig("TPSpmfOGprot.png")
    plt.close('all')

def plot_prot_movement_with_error(bigstruct):
    fsize = 22
    plt.figure(figsize=(fsize,fsize))
    pnum = [[1,3,5],[2,4,6], [7,9,11], [8,10,12]]
    for i in range(4):
        plt.subplot(6, 2, pnum[i][0])
        kcxavg, kcxuperr, kcxlowerr = bigstruct[i].get_kcx_data()
        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, kcxuperr, vmin=0, vmax=1)
        plt.title("UpperError for " + bigstruct[i].moltype + " : KCX protonation probability")
        plt.colorbar()
        plt.subplot(6,2, pnum[i][1])
        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, kcxavg, vmin=0, vmax=1)
        plt.title(bigstruct[i].moltype + " : KCX protonation probability")
        plt.colorbar()
        plt.subplot(6,2, pnum[i][2])
        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, kcxlowerr, vmin=0, vmax=1)
        plt.title("LowerError for " + bigstruct[i].moltype + " : KCX protonation probability")
        plt.colorbar()
    plt.draw()
    plt.savefig("TPSpmfKCXprot.png")

    plt.figure(figsize=(fsize,fsize))
    pnum = [[1,3,5],[2,4,6], [7,9,11], [8,10,12]]
    for i in range(4):
        plt.subplot(6, 2, pnum[i][0])
        ogavg, oguperr, oglowerr = bigstruct[i].get_og_data()
        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, oguperr, vmin=0, vmax=1)
        plt.title("UpperError for " + bigstruct[i].moltype + " : OG protonation probability")
        plt.colorbar()
        plt.subplot(6,2, pnum[i][1])
        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, ogavg, vmin=0, vmax=1)
        plt.title(bigstruct[i].moltype + " : OG protonation probability")
        plt.colorbar()
        plt.subplot(6,2, pnum[i][2])
        plt.pcolor(bigstruct[i].r1points, bigstruct[i].r2points, oglowerr, vmin=0, vmax=1)
        plt.title("LowerError for " + bigstruct[i].moltype + " : OG protonation probability")
        plt.colorbar()
    plt.draw()
    plt.savefig("TPSpmfOGprot.png")
    plt.close('all')

infolist = []
psfpaths = ["/data/sguthrie/imivsdori/dori_sim/sim1/template/dori.psf", "/data/sguthrie/imivsdori/m_imi_sim/from_imi/sim1/template/mimi.psf",
            "/data/sguthrie/imivsdori/sdori_sim/from_dori/sim1/template/sdori.psf", "/data/sguthrie/imivsdori/imi_sim/sim1/template/imi.psf"]
smrefpaths = ["/data/sguthrie/imivsdori/dori_sim/sim2/template/dori.pdb", "/data/sguthrie/imivsdori/m_imi_sim/from_imi/sim1/template/mimi.pdb",
            "/data/sguthrie/imivsdori/sdori_sim/from_dori/sim1/template/sdori.pdb", "/data/sguthrie/imivsdori/imi_sim/sim2/template/imi.pdb"]
rootpaths = ["/data/sguthrie/imivsdori/dori_sim/sim1/tps_getv2", "/data/sguthrie/imivsdori/m_imi_sim/from_imi/sim1/tps_getv",
            "/data/sguthrie/imivsdori/sdori_sim/from_dori/sim1/tps_getv", "/data/sguthrie/imivsdori/imi_sim/sim1/tps_getv"]
isdoris = [True, False, True, False]
moltypes = ['SDR', 'SMI', 'SSD', 'SIM']
picklefiles = ['TPSpmfDori_interactions.pkl', 'TPSpmfMimi_interactions.pkl', 'TPSpmfSDori_interactions.pkl', 'TPSpmfImi_interactions.pkl']
for i in range(4):
    tmp = [psfpaths[i], smrefpaths[i], rootpaths[i], isdoris[i], moltypes[i], picklefiles[i]]
    infolist.append(tmp)
bigrefpath = "/data/sguthrie/imivsdori/OXASetup/RMSD/rmsdcompare.pdb"
#cumul_hbonds structure:
#   overall: dictionary indexed by (r1 value*100, r2 value*100)
#       dictionary hashes to list.
#       0: integer - number of simulations that saved at that (r1,r2) value
#       1: integer - number of simulations that had a proton on KCX
#       2: integer - number of simulations that had a proton on OG
#       3: dictionary - OAI interactions
#       4: dictionary - OAD interactions
#       5: dictionary - KCX interactions
#       6: dictionary - Water interactions
#       7: dictionary - OG interactions
#       8: dictionary - Tail interactions  


# keylists: OAI, OAD, KCX, WAT, OG, Tail

logregex1 = '(-?[0-9]\.[0-9]+) (-?[0-9]\.[0-9]+)'
logregex2 = '(?<=#)(/[-[0-9\.\w]+)+'


try:
    inp = open('TPSpmfInteractions.pkl', 'rb')
    bigstruct = pickle.load(inp)
    all_key_lists = pickle.load(inp)
    inp.close()
except IOError:
    bigstruct = []
    all_key_lists = []
    for foo in range(4):
        psfpath, smrefpath, rootpath, isdori, moltype, pfile = infolist[foo]
        print rootpath
        try:
            inp = open(pfile, 'rb')
            cumul_hbonds = pickle.load(inp)
            keylist = pickle.load(inp)
            inp.close()
        except IOError:
            # keylists: OAI, OAD, KCX, WAT, OG, tail
            keylist = [[] for x in range(6)]
            matrix = []
            for y in range(-60, -5):
                for x in range(-60, 60):
                    matrix.append((x,y))

            cumul_hbonds = {key:[0, 0, 0, {}, {}, {}, {}, {}, {}] for key in matrix}
            for root, dirs, files in os.walk(rootpath):
                if len(dirs) == 0:
                    print root
                    if "tpsv_trajs.dcd" in files:
                        dcdpath = root + "/tpsv_trajs.dcd"
                        traj = analyze_big_traj(dcdpath, psfpath, isdori, moltype)
                        total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds, tailhbonds = traj.analyze_hbonds()
                        hbondlists = [OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds, tailhbonds]
                        for rvals in total:
                            if total[rvals][0] != 0:
                                sm_r1 = round(rvals[0]/10.)
                                sm_r2 = round(rvals[1]/10.)
                                r_rvals = (sm_r1, sm_r2)
                                cumul_hbonds[r_rvals][0] += total[rvals][0]
                                cumul_hbonds[r_rvals][1] += total[rvals][1]
                                cumul_hbonds[r_rvals][2] += total[rvals][2]
                                for i in range(6):
                                    for inter in hbondlists[i][rvals]:
                                        if inter not in keylist[i]:
                                            keylist[i].append(inter)
                                        block_incr_dict(hbondlists[i][rvals][inter], cumul_hbonds[r_rvals][i+3], inter)
                                    keylist[i].sort()

            output = open(pfile, 'wb')
            pickle.dump(cumul_hbonds, output, -1)
            pickle.dump(keylist, output, -1)
            output.close()
            
        all_key_lists.append(keylist)
        data = Trajdata(moltype, cumul_hbonds, keylist, 25)
        bigstruct.append(data)

    output = open('TPSpmfInteractions.pkl', 'wb') 
    pickle.dump(bigstruct, output, -1)
    pickle.dump(all_key_lists, output, -1)
    output.close()

who = ['OAI', 'OAD', 'KCX', 'WAT', 'OG', 'Tail']
keysets = []
for x in range(6):
    foo = set()
    for y in range(4):
        foo |= set(all_key_lists[y][x])
    keysets.append(foo)

plot_prot_movement(bigstruct)
 
for x in range(6):
    for keyname in keysets[x]:
        plot_2Dhist(bigstruct, x, keyname, who[x] + " : " + keyname, ploterr=False)
        print who[x] + " : " + keyname
