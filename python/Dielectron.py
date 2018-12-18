import os, sys
import ROOT
import miniDst
import cuts
from ROOT import ROOTClass


#open root file
Rootfile = TFile('/star/u/tc88qy/AuAu/run17/54GeV/picoDst/miniTree/rootfiles_PicoDst/AF8BE685DC62E95B8D1158145E66157E_0.root')

#retrieve the ntuple of interest
Rootchain = Rootfile.Get('miniDst')
entries = Rootchain.GetEntriesFast()

for jentry in xrange(entires):
	#get the next tree in the chain and verify
	ientry = Rootchain.LoadTree(jentry)
	if ientry < 0 :
		break
	
	#copy next entry into memory and verify
	nb = Rootchain.GetEntry(jentry)
	if nb<=0:
		continue 
	
	#use the values directly from the tree
	nEvent = int(Rootchain.ev)
	if nEvent<0:
		continue
