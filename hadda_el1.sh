hadd -f WWTree_data_el_2016_runC_v1_el.root WWTree_data_el_2016_runC_v1_*_el.root
hadd -f WWTree_data_el_2016_runD_v1_el.root WWTree_data_el_2016_runD_v1_*_el.root
hadd -f WWTree_data_el_2016_runF_v1_el.root WWTree_data_el_2016_runF_v1_*_el.root
hadd -f WWTree_data_el_2016_runG_v1_el.root WWTree_data_el_2016_runG_v1_*_el.root
hadd -f WWTree_data_el_2016_runH_v2_el.root WWTree_data_el_2016_runH_v2_*_el.root
hadd -f WWTree_data_golden_C1D1F1G1H23.root WWTree_data_el_2016_runC_v1_el.root WWTree_data_el_2016_runD_v1_el.root WWTree_data_el_2016_runF_v1_el.root WWTree_data_el_2016_runG_v1_el.root WWTree_data_el_2016_runH_v2_el.root WWTree_data_el_2016_runH_v3_el.root
mv WWTree_data_el_2016_runC_v1_*_el.root others/
mv WWTree_data_el_2016_runD_v1_*_el.root others/
mv WWTree_data_el_2016_runF_v1_*_el.root others/
mv WWTree_data_el_2016_runG_v1_*_el.root others/
mv WWTree_data_el_2016_runH_v2_*_el.root others/
mv WWTree_data_el_2016_runC_v1_el.root WWTree_data_el_2016_runD_v1_el.root WWTree_data_el_2016_runF_v1_el.root WWTree_data_el_2016_runG_v1_el.root WWTree_data_el_2016_runH_v2_el.root WWTree_data_el_2016_runH_v3_el.root others/
hadd -f WWTree_DYJetsToLL_amcatnlo_ext1_el.root WWTree_DYJetsToLL_amcatnlo_ext1_1_el.root WWTree_DYJetsToLL_amcatnlo_ext1_2_el.root
hadd -f WWTree_Signal_LLpLTpTT_el.root	WWTree_Signal_LL_el.root WWTree_Signal_LT_el.root WWTree_Signal_TT_el.root
hadd -f WWTree_tch_el.root	WWTree_tch_1_el.root	WWTree_tch_2_el.root	WWTree_tch_3_el.root	WWTree_tch_4_el.root	WWTree_tch_5_el.root	
hadd -f WWTree_tch_bar_el.root	WWTree_tch_bar_1_el.root	WWTree_tch_bar_2_el.root	WWTree_tch_bar_3_el.root	
hadd -f WWTree_TTbar_amcatnlo_el.root	WWTree_TTbar_amcatnlo_1_el.root	WWTree_TTbar_amcatnlo_2_el.root	WWTree_TTbar_amcatnlo_3_el.root	WWTree_TTbar_amcatnlo_4_el.root	
hadd -f WWTree_TTbar_powheg_el.root	WWTree_TTbar_powheg_1_el.root	WWTree_TTbar_powheg_2_el.root	WWTree_TTbar_powheg_3_el.root	WWTree_TTbar_powheg_4_el.root	
hadd -f WWTree_tWch_bara_el.root	WWTree_tWch_bar_el.root	WWTree_tWch_bar_ext1_el.root	
hadd -f WWTree_tWcha_el.root	WWTree_tWch_el.root	WWTree_tWch_ext1_el.root	
hadd -f WWTree_WJets100a_el.root	WWTree_WJets100_el.root	WWTree_WJets100ext1_el.root	WWTree_WJets100ext2_1_el.root	WWTree_WJets100ext2_2_el.root	
hadd -f WWTree_WJets1200a_el.root	WWTree_WJets1200_el.root	WWTree_WJets1200ext1_el.root	
hadd -f WWTree_WJets200a_el.root	WWTree_WJets200_el.root	WWTree_WJets200ext1_el.root	WWTree_WJets200ext2_el.root	
hadd -f WWTree_WJets2500a_el.root	WWTree_WJets2500_el.root	WWTree_WJets2500ext1_el.root	
hadd -f WWTree_WJets400a_el.root	WWTree_WJets400_el.root	WWTree_WJets400ext1_el.root	
hadd -f WWTree_WJets600a_el.root	WWTree_WJets600_el.root	WWTree_WJets600ext1_el.root	
hadd -f WWTree_WW_excla_el.root	WWTree_WW_excl_el.root	WWTree_WW_excl_ext1_el.root	
hadd -f WWTree_WZ_excl_amcatnlo_el.root	WWTree_WZ_excl_amcatnlo_1_el.root	WWTree_WZ_excl_amcatnlo_2_el.root	

mv WWTree_DYJetsToLL_amcatnlo_ext1_1_el.root WWTree_DYJetsToLL_amcatnlo_ext1_2_el.root others/
mv WWTree_Signal_LL_el.root WWTree_Signal_LT_el.root WWTree_Signal_TT_el.root others/
mv WWTree_tch_1_el.root	WWTree_tch_2_el.root	WWTree_tch_3_el.root	WWTree_tch_4_el.root	WWTree_tch_5_el.root	 others/
mv WWTree_tch_bar_1_el.root	WWTree_tch_bar_2_el.root	WWTree_tch_bar_3_el.root	 others/
mv WWTree_TTbar_amcatnlo_1_el.root	WWTree_TTbar_amcatnlo_2_el.root	WWTree_TTbar_amcatnlo_3_el.root	WWTree_TTbar_amcatnlo_4_el.root	 others/
mv WWTree_TTbar_powheg_1_el.root	WWTree_TTbar_powheg_2_el.root	WWTree_TTbar_powheg_3_el.root	WWTree_TTbar_powheg_4_el.root	 others/
mv WWTree_tWch_bar_el.root	WWTree_tWch_bar_ext1_el.root	 others/
mv WWTree_tWch_el.root	WWTree_tWch_ext1_el.root	 others/
mv WWTree_WJets100_el.root	WWTree_WJets100ext1_el.root	WWTree_WJets100ext2_1_el.root	WWTree_WJets100ext2_2_el.root	 others/
mv WWTree_WJets1200_el.root	WWTree_WJets1200ext1_el.root	 others/
mv WWTree_WJets200_el.root	WWTree_WJets200ext1_el.root	WWTree_WJets200ext2_el.root	 others/
mv WWTree_WJets2500_el.root	WWTree_WJets2500ext1_el.root	 others/
mv WWTree_WJets400_el.root	WWTree_WJets400ext1_el.root	 others/
mv WWTree_WJets600_el.root	WWTree_WJets600ext1_el.root	 others/
mv WWTree_WW_excl_el.root	WWTree_WW_excl_ext1_el.root	 others/
mv WWTree_WZ_excl_amcatnlo_1_el.root	WWTree_WZ_excl_amcatnlo_2_el.root	 others/
