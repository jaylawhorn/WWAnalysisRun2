# WWAnalysisRun2

The package contains a code to produce ntuple for WW semileptonic final state.
It takes in input ntuples produced from miniAOD with the TreeMaker at this link: https://github.com/lbrianza/RA2_2014


Instructions:

git clone https://github.com/lbrianza/WWAnalysisRun2;

cd WWAnalysisRun2/;

make;

python python/produceWWNtuples.py -n ReducedSelection_RSGraviton4000.root -o RSGraviton4000.root -mc True

# To Do List

1. Add number of jets of AK4 and AK8 after selection.
2. Change/Add VBF & Wjets selections (Try to add 2-3 different type of these selections in same code)
3. Add ZepenFeld Var, Azimuthal
4. Presently electron trigger eff contains the root file for muons so need to change it.
5. Apply muon tracking HIP efficiency
6. Check various efficiency applied in code for ttbar
7. Add some automatic thing to read data from user eos (cernbox) using bjobs
