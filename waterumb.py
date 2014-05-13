#umbsamp_ensemble_analysis with numpy
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

def print_hbonds(b, inter, r1):
    if debug:
        print b
        print "Is classified to be interacting with " + inter
        print "At r1 = ", r1
        
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

class WaterUmbEnsemble:
    def __init__(self, universe, is_dori, moltype, foldername, r1, r1trace):
        """ if is_dori:
                moltype = SDR or SSD
            else:
                moltype = SIM or SMI"""
        self.universe = universe
        self.is_dori = is_dori
        self.moltype = moltype
        self.name = foldername
        self.r1 = r1
        self.r1trace = r1trace

    def analyze_hbonds(self):
        if self.is_dori:
            examinestr = "(atom A 81 OAI) or (atom A 81 OAD) or (atom A 84 OH1) or (atom A 84 OH1) " + \
                         "or (atom A 84 OH2) or (atom W 277 OH2) or (atom A 81 OG) " + \
                         "or (atom A 81 NAO) or (atom A 81 OAG) or (atom A 81 OAF) or (atom A 81 NAC) or (atom A 81 NAP)"
            #Not sure that these tail atoms can strictly be said to hydrogen bond
            #OG is already a donor and acceptor
            new_donors = ['OAI', 'NAN', 'NAO', 'NAC', 'OH1', 'OH2', 'NAP'] 
            new_acceptors = ['OAD', 'OAI', 'NAN', 'OAH', 'OAE', 'NAC', 'NAP', 'OH2', 'OH1']
            OAIacceptor = self.moltype + '81:OAI'
            OAIdonor = self.moltype + '81:HOI'
            OADacceptor = self.moltype + '81:OAD'
            taildonors = [self.moltype + '81:HN1', self.moltype + '81:HN2', self.moltype + '81:HN3', self.moltype + '81:HN4', self.moltype + '81:HNP']
            tailacceptors = [self.moltype + '81:OAG', self.moltype + '81:OAF']
        else:
            examinestr = "(atom A 81 O62) or (atom A 81 O7) or (atom A 84 OH1) or (atom A 84 OH1) " + \
                         "or (atom A 84 OH2) or (atom W 277 OH2) or (atom A 81 OG) " + \
                         "or (atom A 81 N24) or (atom A 81 N26)"
            new_donors = ['O62', 'N4', 'N24', 'N26', 'OH1', 'OH2']
            new_acceptors = ['O7', 'O62', 'O31', 'O32', 'N4', 'N24', 'N26', 'OH2', 'OH1']
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
                                            acceptors=new_acceptors, angle=150.0)
        hana.run()
        h_bond_results = hana.timeseries
        r1lower = int((self.r1 - 2)*100)
        r1higher = int((self.r1 + 2)*100)

        matrix = range(r1lower, r1higher)

        total = {key:[0, 0, 0] for key in matrix}
        OAIhbonds = {key:{} for key in matrix}
        OADhbonds = {key:{} for key in matrix}
        KCXhbonds = {key:{} for key in matrix}
        WAThbonds = {key:{} for key in matrix}
        OGhbonds = {key:{} for key in matrix}
        tailhbonds = {key:{} for key in matrix}

        r_cov = 1.31
        #Fill dictionaries; 0 if bond is not present; 1 if it is
        for tstep, foo in enumerate(zip(self.r1trace, h_bond_results)):
            r1, frame = foo
            r1 = int(round(r1, 2)*100)
            total[r1][0] += 1

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
                    print_hbonds(hbond, "OAIdonor", r1/100.)
                    incr_dict(OAIhbonds[r1], acceptorstr)
                if hbond[2] == WATdonor1:
                    print_hbonds(hbond, "WATdonor", r1/100.)
                    incr_dict(WAThbonds[r1], acceptorstr)
                if hbond[2] == WATdonor2:
                    print_hbonds(hbond, "WATdonor", r1/100.)
                    incr_dict(WAThbonds[r1], acceptorstr)
                if hbond[2] in taildonors:
                    print_hbonds(hbond, "TailDonor", r1/100.)
                    incr_dict(tailhbonds[r1], acceptorstr)
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
                if hbond[3] in tailacceptors:
                    print_hbonds(hbond, "TailAcceptor", r1/100.)
                    incr_dict(tailhbonds[r1], donorstr)
        return total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds, tailhbonds
    
    def plot(self):
        plt.figure()
        plt.plot(self.r1trace, self.r2trace, 'go')
        plt.draw()



    
def make_key_data(keyname, c_hbonds, inter):
    '''Calculates probability and wilson error bars at each point. Returns np masked arrays '''
    z = 1.96
    zsq = 3.84
    probvals = np.ma.masked_array(np.zeros(len(c_hbonds)), mask=np.zeros(len(c_hbonds)))
    upper_err = np.ma.masked_array(np.zeros(len(c_hbonds)), mask=np.zeros(len(c_hbonds)))
    lower_err = np.ma.masked_array(np.zeros(len(c_hbonds)), mask=np.zeros(len(c_hbonds)))
    rpointlist = c_hbonds.keys()
    rpointlist.sort()
    for i, rpoint in enumerate(rpointlist):
        total = float(c_hbonds[rpoint][0])
        if total == 0:
            probvals[i] = np.ma.masked
            upper_err[i] = np.ma.masked
            lower_err[i] = np.ma.masked
        else:
            if keyname in c_hbonds[rpoint][inter+3]:
                prob = c_hbonds[rpoint][inter+3][keyname]/total
            else:
                prob = 0
            wscorepm = z*math.sqrt(prob*(1-prob)/total + zsq/(4*total**2))
            wmult = 1/(1+zsq/total)
            probvals[i] = prob
            upper_err[i] = wmult*(prob + zsq/(2*total) + wscorepm)
            lower_err[i] = wmult*(prob + zsq/(2*total) - wscorepm)
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

class WaterUmbdata:
    def __init__(self, moltype, cumul_hbonds, keylist):
        z = 1.96
        zsq = 3.84
        self.moltype = moltype
        self.keylist = keylist
        self.r1points = np.zeros(len(cumul_hbonds))
        self.nsims = np.zeros(len(cumul_hbonds))
        self.kcx_attached = np.ma.masked_array(np.zeros(len(cumul_hbonds)), mask=np.zeros(len(cumul_hbonds)))
        self.kcx_upper_err = np.ma.masked_array(np.zeros(len(cumul_hbonds)), mask=np.zeros(len(cumul_hbonds)))
        self.kcx_lower_err = np.ma.masked_array(np.zeros(len(cumul_hbonds)), mask=np.zeros(len(cumul_hbonds)))
        self.og_attached = np.ma.masked_array(np.zeros(len(cumul_hbonds)), mask=np.zeros(len(cumul_hbonds)))
        self.og_upper_err = np.ma.masked_array(np.zeros(len(cumul_hbonds)), mask=np.zeros(len(cumul_hbonds)))
        self.og_lower_err = np.ma.masked_array(np.zeros(len(cumul_hbonds)), mask=np.zeros(len(cumul_hbonds)))
        rpointlist = cumul_hbonds.keys()
        rpointlist.sort()
        for i, rpoint in enumerate(rpointlist):
            self.r1points[i] = rpoint/100.
            total = cumul_hbonds[rpoint][0]
            self.nsims[i] = total
            if cumul_hbonds[rpoint][0] == 0:
                self.kcx_attached[i] = np.ma.masked
                self.og_attached[i] = np.ma.masked
                self.kcx_upper_err[i] = np.ma.masked
                self.kcx_lower_err[i] = np.ma.masked
                self.og_upper_err[i] = np.ma.masked
                self.og_lower_err[i] = np.ma.masked
            else:
                total = float(total)
                probkcx = cumul_hbonds[rpoint][1]/total
                probog = cumul_hbonds[rpoint][2]/total
                kcx_wscorepm = z*math.sqrt(probkcx*(1-probkcx)/total + zsq/(4*total**2))
                og_wscorepm = z*math.sqrt(probog*(1-probog)/total + zsq/(4*total**2))
                wmult = 1/(1+zsq/total)
            
                self.kcx_attached[i] = probkcx
                self.og_attached[i] = probog
                self.kcx_upper_err[i] = wmult*(probkcx + zsq/(2*total) + kcx_wscorepm)
                self.kcx_lower_err[i] = wmult*(probkcx + zsq/(2*total) - kcx_wscorepm)
                self.og_upper_err[i] = wmult*(probog + zsq/(2*total) + og_wscorepm)
                self.og_lower_err[i] = wmult*(probog + zsq/(2*total) - og_wscorepm)

        self.OAI = {}
        self.OAD = {}
        self.KCX = {}
        self.WAT = {}
        self.OG = {}
        self.tail = {}
        self.interdicts = [self.OAI, self.OAD, self.KCX, self.WAT, self.OG, self.tail]

        for i in range(6):
            for keyname in self.keylist[i]:
                probvals, upper_err, lower_err = make_key_data(keyname, cumul_hbonds, i)
                self.interdicts[i][keyname] = [probvals, upper_err, lower_err]

    def get_plot_data(self, inter_i, keyname):
        if keyname in self.interdicts[inter_i]:
            return self.interdicts[inter_i][keyname]
        elif (keyname in ['SDR81', 'SIM81']) and (self.moltype + '81' in self.interdicts[inter_i]):
            return self.interdicts[inter_i][self.moltype + '81']
        else:
            return [None, None, None]
    def get_kcx_data(self):
        return self.kcx_attached, self.kcx_upper_err, self.kcx_lower_err
    def get_og_data(self):
        return self.og_attached, self.og_upper_err, self.og_lower_err
    def get_total(self):
        return self.nsims


def get_dist(atom1, atom2):
    return math.sqrt(sum((atom1[i] - atom2[i])**2 for i in range(3)))

def analyze_umbsamp(dcd_filepath, psf_filepath, is_dori, moltype, r1):
    universe = MDAnalysis.Universe(psf_filepath, dcd_filepath)
    if is_dori:
        aoi = universe.selectAtoms("(atom A 81 CAJ) or (atom W 277 OH2) ")
    else:
        aoi = universe.selectAtoms("(atom A 81 C7) or (atom W 277 OH2)")
    r1trace = []
    rxnatoms = universe.trajectory.timeseries(aoi)
    for t in range(rxnatoms.shape[1]):
        m = get_dist(rxnatoms[0][t], rxnatoms[1][t])
        r1trace.append(-m)

    ensemble = WaterUmbEnsemble(universe, is_dori, moltype, dcd_filepath, r1, r1trace)
    return ensemble

def smooth_avg(avg):
    smoothed = np.ma.copy(avg)
    maskarr = np.ma.getmask(avg)
    for i in range(len(avg)):
        if not maskarr[i]:
            tmp = []
            for across in range(i-2, i+3):
                if not np.ma.getmask(avg[across]):
                    tmp.append(avg[across])
            if len(tmp) == 0:
                smoothed[i] = np.ma.masked
            else:
                smoothed[i] = sum(tmp)/float(len(tmp))
    return smoothed
    
def smooth_all(avg, uperr, lowerr):
    smoothed = np.ma.copy(avg)
    superr = np.ma.copy(uperr)
    slowerr = np.ma.copy(lowerr)
    maskarr = np.ma.getmask(avg)
    for i in range(len(avg)):
        if not maskarr[i]:
            tmp = []
            ltmp = []
            utmp = []
            for across in range(i-2, i+3):
                if not np.ma.getmask(avg[across]):
                    tmp.append(avg[across])
                    utmp.append(uperr[across])
                    ltmp.append(lowerr[across])
            if len(tmp) == 0:
                smoothed[i] = np.ma.masked
                superr[i] = np.ma.masked
                slowerr[i] = np.ma.masked
            else:
                smoothed[i] = sum(tmp)/float(len(tmp))
                superr[i] = sum(utmp)/float(len(utmp))
                slowerr[i] = sum(ltmp)/float(len(ltmp))
    return smoothed, superr, slowerr

def plot_total(bigstruct):
    plt.figure(figsize=(14,8))
    colors = ['b', 'g']
    for i in range(2):
        plt.subplot(2, 1, i+1)
        plt.grid()
        total = bigstruct[i].get_total()
        w = bigstruct[i].r1points[1] - bigstruct[i].r1points[0]
        plt.bar(bigstruct[i].r1points, total, width=w, color=colors[i])
        plt.xlabel('R1')
        plt.ylabel('Number of Simulations')
        plt.axis([-5.5, -1.5, 0, 35])
    plt.subplot(2, 1, 1)
    plt.title('Number of Simulations at Reaction Coordinate Point')
    plt.draw()
    plt.savefig("NumSimsRatchetUmb.png")
    plt.close()

def plot_all_scatter(bigstruct, inter_i, keyname, figname, smooth=True):
    plt.figure(figsize=(14,8))
    colors = ['b', 'g']
    for i in range(2):
        [avg, uperr, lowerr] = bigstruct[i].get_plot_data(inter_i, keyname)
        if avg != None:
            if smooth:
                avg, uperr, lowerr = smooth_all(avg, uperr, lowerr)
            asymmerr = [lowerr, uperr]
            plt.errorbar(bigstruct[i].r1points, avg, yerr=asymmerr, fmt='-o', mfc=colors[i], label=bigstruct[i].moltype)
    plt.legend()
    plt.xlabel('R1')
    plt.ylabel('Probability of Interaction')
    plt.title(figname)
    plt.draw()
    if smooth:
        plt.savefig("SmoothedWaterUmb" + figname + ".png")
    else:
        plt.savefig("WaterUmb" + figname + ".png")
    plt.close()


def plot_all(bigstruct, inter_i, keyname, figname, smooth=True):
    plt.figure(figsize=(14,8))
    colors = ['b', 'g']
    for i in range(2):
        plt.subplot(2, 1, i+1)
        plt.grid()
        w = bigstruct[i].r1points[1] - bigstruct[i].r1points[0]
        [avg, uperr, lowerr] = bigstruct[i].get_plot_data(inter_i, keyname)
        if avg != None:
            if smooth:
                avg, uperr, lowerr = smooth_all(avg, uperr, lowerr)
            asymmerr = [lowerr, uperr]
            plt.bar(bigstruct[i].r1points, avg, width=w, color=colors[i])
        plt.xlabel('R1')
        plt.ylabel('Number of Simulations')
        plt.axis([-5.5, -1.5, 0, 1])
    plt.subplot(2, 1, 1)
    plt.title(figname)
    plt.draw()
    if smooth:
        plt.savefig("SmoothedWaterUmb" + figname + ".png")
    else:
        plt.savefig("WaterUmb" + figname + ".png")
    plt.close()
    

infolist = []
psfpaths = ["/data/sguthrie/imivsdori/dori_sim/explore_wat/2watumb/template/dori.psf",
            "/data/sguthrie/imivsdori/imi_sim/explore_wat/2watumb2/template/imi.psf"]
rootpaths = ["/data/sguthrie/imivsdori/dori_sim/explore_wat/2watumb",
             "/data/sguthrie/imivsdori/imi_sim/explore_wat/2watumb2"]
isdoris = [True, False]
moltypes = ['SDR', 'SIM']
picklefiles = ['WaterUmbDori_interactions.pkl', 'WaterUmbImi_interactions.pkl']
for i in range(2):
    tmp = [psfpaths[i], rootpaths[i], isdoris[i], moltypes[i], picklefiles[i]]
    infolist.append(tmp)

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

logregex1 = '(?<=#)(-?[0-9]\.[0-9]+)'
logregex2 = '(?<=#)(/[-[0-9\.\w]+)+'


try:
    inp = open('WaterUmbInteractions.pkl', 'rb')
    bigstruct = pickle.load(inp)
    all_key_lists = pickle.load(inp)
    inp.close()
except IOError:
    bigstruct = []
    all_key_lists = []
    for foo in range(2):
        psfpath, rootpath, isdori, moltype, pfile = infolist[foo]
        print rootpath
        try:
            inp = open(pfile, 'rb')
            cumul_hbonds = pickle.load(inp)
            keylist = pickle.load(inp)
            inp.close()
        except IOError:
            # keylists: OAI, OAD, KCX, WAT, OG, tail
            keylist = [[] for x in range(6)]

            cumul_hbonds = {key:[0, 0, 0, {}, {}, {}, {}, {}, {}] for key in range(-550, -150)}
            logfile = open(rootpath + "/umbsamp.log", "r")
            i = 0
            for line in logfile:
                m = re.search(logregex1, line)
                if m != None:
                    r1 = m.group(0)
                    x = re.search(logregex2, line)
                    dcdpath = x.group(0) + "/cap_production.dcd"
                    try:
                        #Check if dcd exists
                        print dcdpath
                        f = open(dcdpath, 'r')
                        f.close()
                        ensemble = analyze_umbsamp(dcdpath, psfpath, isdori, moltype, float(r1))
                        total, OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds, tailhbonds = ensemble.analyze_hbonds()
                        hbondlists = [OAIhbonds, OADhbonds, KCXhbonds, WAThbonds, OGhbonds, tailhbonds]
                        for rvals in total:
                            if total[rvals][0] != 0:
                                cumul_hbonds[rvals][0] += total[rvals][0]
                                cumul_hbonds[rvals][1] += total[rvals][1]
                                cumul_hbonds[rvals][2] += total[rvals][2]
                                for i in range(6):
                                    for inter in hbondlists[i][rvals]:
                                        if inter not in keylist[i]:
                                            keylist[i].append(inter)
                                        block_incr_dict(hbondlists[i][rvals][inter], cumul_hbonds[rvals][i+3], inter)
                                    keylist[i].sort()
                    except IOError:
                        pass
            output = open(pfile, 'wb')
            pickle.dump(cumul_hbonds, output, -1)
            pickle.dump(keylist, output, -1)
            output.close()
            
        all_key_lists.append(keylist)
        data = WaterUmbdata(moltype, cumul_hbonds, keylist)
        bigstruct.append(data)

    output = open('WaterUmbInteractions.pkl', 'wb') 
    pickle.dump(bigstruct, output, -1)
    pickle.dump(all_key_lists, output, -1)
    output.close()

who = ['OAI', 'OAD', 'KCX', 'WAT', 'OG', 'Tail']
keysets = []
for x in range(6):
    foo = set()
    for y in range(2):
        foo |= set(all_key_lists[y][x])
    keysets.append(foo)

plot_total(bigstruct)

for x in range(6):
    for keyname in keysets[x]:
        plot_all(bigstruct, x, keyname, who[x] + " : " + keyname, smooth=True)
        plot_all(bigstruct, x, keyname, who[x] + " : " + keyname, smooth=False)
        print who[x] + " : " + keyname





