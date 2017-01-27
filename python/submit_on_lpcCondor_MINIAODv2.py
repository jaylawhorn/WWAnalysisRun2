#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir =  currentDir+"/../";

#inputFolder = "/eos/cms/store/caf/user/lbrianza/WWReducedTree_run2";
inputFolder = "/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/Jan102016";
outputFolder = currentDir+"/output/";
exeName = "produceWWNtuples.exe";
exePathName = currentDir+"/"+exeName;

dryRun = False;
doMC = False;
doData = True;

category = ["mu","el"];
#category = ["el"];
#category = ["mu"];

#data_el_2016_runC_v1  data_mu_2016_runD_v1      Signal_LL       tWch           WJets1200      WJets400        WW_excl
#data_el_2016_runE_v1  data_mu_2016_runE_v1      Signal_LT       tWch_bar       WJets1200ext1  WJets400ext1    WW_excl_amcatnlo
#data_el_2016_runF_v1  data_mu_2016_runF_v1      Signal_TT       tWch_bar_ext1  WJets200       WJets600        WW_excl_ext1
#data_el_2016_runG_v1  data_mu_2016_runG_v1      tch             tWch_ext1      WJets200ext1   WJets600ext1    WZ_excl_amcatnlo
#data_el_2016_runH_v2  data_mu_2016_runH_v3      tch_bar         WJets100       WJets200ext2   WJets800ext1    ZZ_excl_amcatnlo
#data_el_2016_runH_v3  DYJetsToLL_amcatnlo_ext1  TTbar_amcatnlo  WJets100ext1   WJets2500      WJets_amcatnlo
#data_mu_2016_runB_v3  sch                       TTbar_powheg    WJets100ext2   WJets2500ext1  WJets_madgraph


samples = [
    ( 0.25973,		"Signal_LT",		2471400,	0.),
    ( 0.04130,		"Signal_LL",		499800,		0.),
    ( 0.44750,		"Signal_TT",		4988000,	0.),
    ( 49.9970,		"WW_excl",		1999200,	0.),
    ( 49.9970,		"WW_excl_ext1",		6998600,	0.),
    ( 49.9970,		"WW_excl_amcatnlo",	5176114,	0.),
    ( 10.7100,		"WZ_excl_amcatnlo_2",	24221923,	0.),
    ( 10.7100,		"WZ_excl_amcatnlo_1",	24221923,	0.),
    ( 3.22000,		"ZZ_excl_amcatnlo",	15345572,	0.),
    ( 3.36000,		"sch",			1000000,	0.),
    ( 43.8000,		"tch_5",		67240808,	0.),
    ( 43.8000,		"tch_4",		67240808,	0.),
    ( 43.8000,		"tch_3",		67240808,	0.),
    ( 43.8000,		"tch_2",		67240808,	0.),
    ( 43.8000,		"tch_1",		67240808,	0.),
    ( 26.0700,          "tch_bar_3",  		38811017,	0.),
    ( 26.0700,          "tch_bar_2",  		38811017,	0.),
    ( 26.0700,          "tch_bar_1",  		38811017,	0.),
    ( 35.6000,		"tWch_ext1",		6952830,	0.),
    ( 35.6000,		"tWch_bar_ext1",	6933094,	0.),
    ( 831.760,		"TTbar_amcatnlo_4", 	23561608,    	0.),
    ( 831.760,		"TTbar_amcatnlo_3", 	23561608,    	0.),
    ( 831.760,		"TTbar_amcatnlo_2", 	23561608,    	0.),
    ( 831.760,		"TTbar_amcatnlo_1", 	23561608,    	0.),
    ( 831.760,		"TTbar_powheg", 	77229341,    	0.),
    ( 1627.45,       	"WJets100", 		10235198,    	0.),
    ( 1627.45,       	"WJets100ext1", 	29503700,    	0.),
    ( 1627.45,       	"WJets100ext2_2", 	39617787,    	0.),
    ( 1627.45,       	"WJets100ext2_1", 	39617787,    	0.),
    ( 435.24,       	"WJets200", 		4950373,    	0.),
    ( 435.24,       	"WJets200ext1", 	14815928,    	0.),
    ( 435.24,       	"WJets200ext2", 	19914590,    	0.),
    ( 59.18,       	"WJets400", 		1963464,    	0.),
    ( 59.18,       	"WJets400ext1", 	5796237,    	0.),
    ( 14.58,       	"WJets600", 		3779141,    	0.),
    ( 14.58,       	"WJets600ext1", 	14908339,    	0.),
    ( 6.655,       	"WJets800ext1", 	6200954,    	0.),
    ( 1.60809,       	"WJets1200", 		244532,    	0.),
    ( 1.60809,       	"WJets1200ext1", 	6627909,    	0.),
    ( 0.0389136,       	"WJets2500", 		253561,    	0.),
    ( 0.0389136,       	"WJets2500ext1", 	2384260,    	0.)
    ]

nameDataMu = [
"data_mu_2016_runB_v3_1",
"data_mu_2016_runB_v3_10",
"data_mu_2016_runB_v3_11",
"data_mu_2016_runB_v3_12",
"data_mu_2016_runB_v3_13",
"data_mu_2016_runB_v3_14",
"data_mu_2016_runB_v3_15",
"data_mu_2016_runB_v3_16",
"data_mu_2016_runB_v3_17",
"data_mu_2016_runB_v3_18",
"data_mu_2016_runB_v3_19",
"data_mu_2016_runB_v3_2",
"data_mu_2016_runB_v3_3",
"data_mu_2016_runB_v3_4",
"data_mu_2016_runB_v3_5",
"data_mu_2016_runB_v3_6",
"data_mu_2016_runB_v3_7",
"data_mu_2016_runB_v3_8",
"data_mu_2016_runB_v3_9",
"data_mu_2016_runC_v1",
"data_mu_2016_runC_v1_1",
"data_mu_2016_runC_v1_2",
"data_mu_2016_runC_v1_3",
"data_mu_2016_runC_v1_4",
"data_mu_2016_runC_v1_5",
"data_mu_2016_runC_v1_6",
"data_mu_2016_runC_v1_7",
"data_mu_2016_runD_v1_1",
"data_mu_2016_runD_v1_10",
"data_mu_2016_runD_v1_11",
"data_mu_2016_runD_v1_12",
"data_mu_2016_runD_v1_2",
"data_mu_2016_runD_v1_3",
"data_mu_2016_runD_v1_4",
"data_mu_2016_runD_v1_5",
"data_mu_2016_runD_v1_6",
"data_mu_2016_runD_v1_7",
"data_mu_2016_runD_v1_8",
"data_mu_2016_runD_v1_9",
"data_mu_2016_runE_v1_1",
"data_mu_2016_runE_v1_10",
"data_mu_2016_runE_v1_2",
"data_mu_2016_runE_v1_3",
"data_mu_2016_runE_v1_4",
"data_mu_2016_runE_v1_5",
"data_mu_2016_runE_v1_6",
"data_mu_2016_runE_v1_7",
"data_mu_2016_runE_v1_8",
"data_mu_2016_runE_v1_9",
"data_mu_2016_runF_v1_1",
"data_mu_2016_runF_v1_2",
"data_mu_2016_runF_v1_3",
"data_mu_2016_runF_v1_4",
"data_mu_2016_runF_v1_5",
"data_mu_2016_runF_v1_6",
"data_mu_2016_runF_v1_7",
"data_mu_2016_runF_v1_8",
"data_mu_2016_runG_v1_1",
"data_mu_2016_runG_v1_10",
"data_mu_2016_runG_v1_11",
"data_mu_2016_runG_v1_12",
"data_mu_2016_runG_v1_13",
"data_mu_2016_runG_v1_14",
"data_mu_2016_runG_v1_15",
"data_mu_2016_runG_v1_16",
"data_mu_2016_runG_v1_17",
"data_mu_2016_runG_v1_2",
"data_mu_2016_runG_v1_3",
"data_mu_2016_runG_v1_4",
"data_mu_2016_runG_v1_5",
"data_mu_2016_runG_v1_6",
"data_mu_2016_runG_v1_7",
"data_mu_2016_runG_v1_8",
"data_mu_2016_runG_v1_9",
"data_mu_2016_runH_v3"
    ];

nameDataEl = [
"data_el_2016_runC_v1_1",
"data_el_2016_runC_v1_10",
"data_el_2016_runC_v1_11",
"data_el_2016_runC_v1_12",
"data_el_2016_runC_v1_13",
"data_el_2016_runC_v1_14",
"data_el_2016_runC_v1_2",
"data_el_2016_runC_v1_3",
"data_el_2016_runC_v1_4",
"data_el_2016_runC_v1_5",
"data_el_2016_runC_v1_6",
"data_el_2016_runC_v1_7",
"data_el_2016_runC_v1_8",
"data_el_2016_runC_v1_9",
"data_el_2016_runD_v1_1",
"data_el_2016_runD_v1_10",
"data_el_2016_runD_v1_11",
"data_el_2016_runD_v1_12",
"data_el_2016_runD_v1_2",
"data_el_2016_runD_v1_3",
"data_el_2016_runD_v1_4",
"data_el_2016_runD_v1_5",
"data_el_2016_runD_v1_6",
"data_el_2016_runD_v1_7",
"data_el_2016_runD_v1_8",
"data_el_2016_runD_v1_9",
"data_el_2016_runF_v1_1",
"data_el_2016_runF_v1_10",
"data_el_2016_runF_v1_11",
"data_el_2016_runF_v1_12",
"data_el_2016_runF_v1_13",
"data_el_2016_runF_v1_14",
"data_el_2016_runF_v1_15",
"data_el_2016_runF_v1_16",
"data_el_2016_runF_v1_17",
"data_el_2016_runF_v1_2",
"data_el_2016_runF_v1_3",
"data_el_2016_runF_v1_4",
"data_el_2016_runF_v1_5",
"data_el_2016_runF_v1_6",
"data_el_2016_runF_v1_7",
"data_el_2016_runF_v1_8",
"data_el_2016_runF_v1_9",
"data_el_2016_runG_v1_1",
"data_el_2016_runG_v1_10",
"data_el_2016_runG_v1_11",
"data_el_2016_runG_v1_12",
"data_el_2016_runG_v1_13",
"data_el_2016_runG_v1_14",
"data_el_2016_runG_v1_15",
"data_el_2016_runG_v1_16",
"data_el_2016_runG_v1_17",
"data_el_2016_runG_v1_18",
"data_el_2016_runG_v1_19",
"data_el_2016_runG_v1_2",
"data_el_2016_runG_v1_3",
"data_el_2016_runG_v1_4",
"data_el_2016_runG_v1_5",
"data_el_2016_runG_v1_6",
"data_el_2016_runG_v1_7",
"data_el_2016_runG_v1_8",
"data_el_2016_runG_v1_9",
"data_el_2016_runH_v2_1",
"data_el_2016_runH_v2_10",
"data_el_2016_runH_v2_11",
"data_el_2016_runH_v2_12",
"data_el_2016_runH_v2_13",
"data_el_2016_runH_v2_14",
"data_el_2016_runH_v2_15",
"data_el_2016_runH_v2_16",
"data_el_2016_runH_v2_17",
"data_el_2016_runH_v2_18",
"data_el_2016_runH_v2_19",
"data_el_2016_runH_v2_2",
"data_el_2016_runH_v2_3",
"data_el_2016_runH_v2_4",
"data_el_2016_runH_v2_5",
"data_el_2016_runH_v2_6",
"data_el_2016_runH_v2_7",
"data_el_2016_runH_v2_8",
"data_el_2016_runH_v2_9",
"data_el_2016_runH_v3"
];


inputlist="HLT_Ele27_WPLoose_eta2p1.root, MuonID_Z_RunBCD_prompt80X_7p65.root,MuonIso_Z_RunBCD_prompt80X_7p65.root, MuonTrackingEfficiencyratios.root,PU.root, PUxSynch.root, SingleElectron_csc2015.txt,SingleElectron_ecalscn1043093.txt,SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root, SingleMuon_csc2015.txt,SingleMuon_ecalscn1043093.txt, egammaEffi.txt_SF2D_GSF_tracking.root,egammaEffi.txt_SF2D_ID_ISO_RunB_0pt5fb.root,egammaEffi.txt_SF2D_ID_ISO_RunB_5pt0fb.root, ntupleList.txt,pileupDataRun2015D_72mb.root, pileupDataRun2016B_69mb.root,pileupDataRun2016B_71p3mb.root, pileupDataRun2016B_72mb.root, puppiJecCorr.root,python/produceWWNtuples.py,produceWWNtuples.exe, Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

nameData = {"el": nameDataEl, "mu":nameDataMu};

command = "python produceWWNtuples.py -i "+inputFolder+" $*";

outScript = open("runstep2condor.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'cd '+CMSSWDir);
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+"cd -");
outScript.write("\n"+command);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runstep2condor.sh");

outJDL = open("runstep2condor.jdl","w");
outJDL.write("Executable = runstep2condor.sh\n");
outJDL.write("Universe = vanilla\n");
outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("\n");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
outJDL.write("transfer_input_files = "+inputlist+"\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");

for a in range(len(category)):

    #MC
    if( doMC ):
        for i in range(len(samples)):
            outJDL.write("Output = "+str(samples[i][1])+"_"+category[a]+".stdout\n");
            outJDL.write("Error = "+str(samples[i][1])+"_"+category[a]+".stdout\n");
            outJDL.write("Arguments = -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+"_"+category[a]+" -l "+category[a]+" -w "+str(samples[i][0])+" -no "+str(samples[i][2])+" -mass "+str(samples[i][3])+" --ismc 1 -trig 0\n");
            outJDL.write("Queue\n");
    
    #data
    if( doData ):
        for i in range(len(nameData[category[a]])):
            outJDL.write("Output = "+(nameData[category[a]])[i]+".stdout\n");
            outJDL.write("Error = "+(nameData[category[a]])[i]+".stdout\n");
            outJDL.write("Arguments = -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+"_"+category[a]+" -l "+category[a]+" -w 1. -no 1. -mass 0 --ismc 0 -trig 1\n");
            outJDL.write("Queue\n");

outJDL.close();
print "\"condor_submit runstep2condor.jdl\" to submit";
