today=`date +%d%m%Y_%H%M%S`
rm *.dat
for d in ./LSFJOB*/;do 
	#echo $d
	tempVar=`grep "TreeMaker2/PreSelection" $d/STDOUT | awk '{print $5}'`
	tempFile=`grep "TreeMaker2/PreSelection" $d/STDOUT | awk '{print $7}'`
	if [[ "$tempVar" == "el" ]];then
		#echo $tempVar
		#grep "TreeMaker2/PreSelection" $d/STDOUT | awk '{print $7}' > ${tempVar}_${tempFile}.dat
		awk '{ if (lines > 0) {print; --lines; }} /SUMMARY/ {lines = 14}' < $d/STDOUT | awk 'FNR > 5' | awk '{print $(NF-1)}' >> ${tempVar}_${tempFile}.dat
	fi
	if [[ "$tempVar" == "mu" ]];then
		#echo $tempVar
		#grep "TreeMaker2/PreSelection" $d/STDOUT | awk '{print $7}' > ${tempVar}_${tempFile}.dat
		awk '{ if (lines > 0) {print; --lines; }} /SUMMARY/ {lines = 14}' < $d/STDOUT | awk 'FNR > 5' | awk '{print $(NF-1)}' >> ${tempVar}_${tempFile}.dat
	fi
done


##	Sum first column of all the given files
#	awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'  data_el_*.dat
echo  "Adding all data for Electrons...."
echo  "data_el_2016" > el_data_evt.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	el_data_el_2016_run*.dat	>>	el_data_evt.dat
echo  "TTbar_powheg" > el_TTbar_evt.dat 
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	el_TTbar_powheg_*.dat	>>	el_TTbar_evt.dat
echo  "WZ" > el_WZ_evt.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	el_WZ_excl_*.dat	>>	el_WZ_evt.dat
echo  "ZZ" > el_ZZ_evt.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	el_ZZ_excl_*.dat	>>	el_ZZ_evt.dat
echo  "WJets" > el_WJetEvt.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	el_WJets*.dat		>>	el_WJetEvt.dat

(echo WW_excl;	cat el_WW_excl.dat;)	>	el_WW_excl_evt.dat
(echo sch;	cat el_sch.dat;)	>	el_sch_evt.dat
(echo tch_bar;	cat el_tch_bar.dat;)	>	el_tch_bar_evt.dat
(echo tWch;	cat el_tWch.dat;)	>	el_tWch_evt.dat
(echo tWch_bar;	cat el_tWch_bar.dat;)	>	el_tWch_bar_evt.dat
(echo DYJetsToLL; cat el_DYJetsToLL_amcatnlo.dat;)	>	el_DYJetsToLL_amcatnlo_evt.dat

#(echo 0;	cat el_WW_excl.dat;)	>	el_WW_excl_evt.dat
#(echo 0;	cat el_sch.dat;)	>	el_sch_evt.dat
#(echo 0;	cat el_tch_bar.dat;)	>	el_tch_bar_evt.dat
#(echo 0;	cat el_tWch.dat;)	>	el_tWch_evt.dat
#(echo 0;	cat el_tWch_bar.dat;)	>	el_tWch_bar_evt.dat
#paste CutsName.txt	el_data_evt.dat el_TTbar_evt.dat el_WJetEvt.dat el_WZ_evt.dat el_ZZ_evt.dat el_WW_excl_evt.dat el_sch_evt.dat el_tch_bar_evt.dat el_tWch_evt.dat el_tWch_bar_evt.dat > Ele_Event_Selection.dat
paste 	el_data_evt.dat el_TTbar_evt.dat el_WJetEvt.dat el_WZ_evt.dat el_ZZ_evt.dat el_WW_excl_evt.dat el_sch_evt.dat el_tch_bar_evt.dat el_tWch_evt.dat el_tWch_bar_evt.dat el_DYJetsToLL_amcatnlo_evt.dat > Ele_Event_Selection.dat

awk 'BEGIN{print "<html>	\n<head>	\n<style>	\ntable, th, td { \n     border: 1px solid black; \n     border-collapse: collapse; \n} \n</style> \n</head> \n<body>	\n<h1>Cut Flow Table For Electron Channel</h1>	\n<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>\n"}' Ele_Event_Selection.dat > CutFlowTable_${today}.htm

#	TO CREATE TEX FILE
#paste CutsName.txt	el_data_evt.dat el_TTbar_evt.dat el_WJetEvt.dat el_WZ_evt.dat el_ZZ_evt.dat el_WW_excl_evt.dat el_sch_evt.dat el_tch_bar_evt.dat el_tWch_evt.dat el_tWch_bar_evt.dat > el_Event_Selection_Tex.dat
#awk 'BEGIN{OFS="\t&\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' el_Event_Selection_Tex.dat >  CutFlowTable_${today}.tex


echo  "Adding all data for Muons...."
echo  "data_mu_2016" > mu_data_evt.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	mu_data_mu_2016_run*.dat	>>	mu_data_evt.dat
echo  "TTbar_powheg" > mu_TTbar_evt.dat 
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	mu_TTbar_powheg_*.dat	>>	mu_TTbar_evt.dat
echo  "WZ" > mu_WZ_evt.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	mu_WZ_excl_*.dat	>>	mu_WZ_evt.dat
echo  "ZZ" > mu_ZZ_evt.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	mu_ZZ_excl_*.dat	>>	mu_ZZ_evt.dat
echo  "WJets" > mu_WJetEvt.dat
awk '{a[FNR]+=$1;} END{for (i=1; i<=FNR; i++) print a[i]}'	mu_WJets*.dat		>>	mu_WJetEvt.dat


(echo WW_excl;	cat mu_WW_excl.dat;)	>	mu_WW_excl_evt.dat
(echo sch;	cat mu_sch.dat;)	>	mu_sch_evt.dat
(echo tch_bar;	cat mu_tch_bar.dat;)	>	mu_tch_bar_evt.dat
(echo tWch;	cat mu_tWch.dat;)	>	mu_tWch_evt.dat
(echo tWch_bar;	cat mu_tWch_bar.dat;)	>	mu_tWch_bar_evt.dat
(echo DYJetsToLL; cat mu_DYJetsToLL_amcatnlo.dat;)	>	mu_DYJetsToLL_amcatnlo_evt.dat


#paste CutsName.txt	mu_data_evt.dat mu_TTbar_evt.dat mu_WJetEvt.dat mu_WZ_evt.dat mu_ZZ_evt.dat mu_WW_excl_evt.dat mu_sch_evt.dat mu_tch_bar_evt.dat mu_tWch_evt.dat mu_tWch_bar_evt.dat > Mu_Event_Selection.dat
paste mu_data_evt.dat mu_TTbar_evt.dat mu_WJetEvt.dat mu_WZ_evt.dat mu_ZZ_evt.dat mu_WW_excl_evt.dat mu_sch_evt.dat mu_tch_bar_evt.dat mu_tWch_evt.dat mu_tWch_bar_evt.dat mu_DYJetsToLL_amcatnlo_evt.dat > Mu_Event_Selection.dat

awk 'BEGIN{print "\n<h1>Cut Flow Table For Muon Channel</h1>	\n<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>\n</body>\n</html>"}' Mu_Event_Selection.dat >> CutFlowTable_${today}.htm

paste CutsName.txt	mu_data_evt.dat mu_TTbar_evt.dat mu_WJetEvt.dat mu_WZ_evt.dat mu_ZZ_evt.dat mu_WW_excl_evt.dat mu_sch_evt.dat mu_tch_bar_evt.dat mu_tWch_evt.dat mu_tWch_bar_evt.dat > Mu_Event_Selection_Tex.dat
awk 'BEGIN{OFS="\t&\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' Mu_Event_Selection_Tex.dat >>  CutFlowTable_${today}.tex
#awk 'BEGIN{print "<html>	\n<head>	\n<style>	\ntable, th, td { \n     border: 1px solid black; \n     border-collapse: collapse; \n} \n</style> \n</head> \n<body>	\n<h1>Cut Flow Table For Muon Channel</h1>	\n<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" $i"</td>";print "</tr>"} END{print "</table>\n</body>\n</html>"}' Mu_Event_Selection.dat > Mu_Event_Selection.htm
