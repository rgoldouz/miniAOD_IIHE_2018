import os
import ROOT
path = []

f = ROOT.TFile.Open("outfile.root")
#f = ROOT.TFile.Open("/afs/crc.nd.edu/user/r/rgoldouz/ExcitedTopAnalysis/analysis/ANoutput.root")
tree_in = f.Get('IIHEAnalysis')
for event in tree_in:
#   if event.Ch==0:
#        print str(event.Run)+","+str(event.Lumi)+","+str(abs(event.Event))
#    if event.Event==12706422i:
    for i in range(len(event.ph_pt)):
        print str(i)+',' + str(event.ph_pt[i]) + ',' + str(event.ph_eta[i])
        print str(i)+',' + str(event.ph_pt[i]*(event.ph_ecalEnergyPostCorr[i]/event.ph_energy[i])) + ',' + str(event.ph_eta[i])
