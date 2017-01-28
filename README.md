# WWAnalysisRun2

The package contains a code to produce ntuple for WW semileptonic final state.
It takes in input ntuples produced from miniAOD with the TreeMaker at this link: https://github.com/ram1123/RA2_2014


Instructions:

	git clone https://github.com/ram1123/WWAnalysisRun2.git;
	cd WWAnalysisRun2/;
	git checkout 80x_devel
	make;
	python python/produceWWNtuples.py -n ReducedSelection_RSGraviton4000.root -o RSGraviton4000.root -mc True

To Do List:

- [ ] Apply latest ID, ISO, & Trigger Efficiency for electron and muons
- [ ] Clean the code
- [ ] make script to extract the number of events from log files
