namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes
 /* 
  TString spectra_file = "./TLee_input_fakeset5_EmuCosthetaEhadron/merge.root";
  TString flux_Xs_directory = "./TLee_input_fakeset5_EmuCosthetaEhadron/flux_Xs/";
  TString detector_directory = "./TLee_input_fakeset5_EmuCosthetaEhadron/det/";
  TString mc_directory = "./TLee_input_fakeset5_EmuCosthetaEhadron/mc_stat/";
 */
  /*
  TString spectra_file = "./TLee_input_fakeset5_EmuCosthetaEhadron_check/merge.root";
  TString flux_Xs_directory = "./TLee_input_fakeset5_EmuCosthetaEhadron_check/flux_Xs/";
  TString detector_directory = "./TLee_input_fakeset5_EmuCosthetaEhadron_check/det/";
  TString mc_directory = "./TLee_input_fakeset5_EmuCosthetaEhadron_check/mc_stat/";
  */


  
  


  
  /*
  TString spectra_file = "./TLee_input_EhadronConstraint_FSI_jun17_fix/merge.root";
  TString flux_Xs_directory = "./TLee_input_EhadronConstraint_FSI_jun17_fix/flux_Xs/";
  TString detector_directory = "./TLee_input_EhadronConstraint_FSI_jun17_fix/det/";
  TString mc_directory = "./TLee_input_EhadronConstraint_FSI_jun17_fix/mc_stat/";
  */


  
  
  /* //TString spectra_file = "./TLee_input_farsideband_EmuCosthetaEhadron/merge_Eproton0d8.root";// too small change */
  /* //TString spectra_file = "./TLee_input_farsideband_EmuCosthetaEhadron/merge_Eproton0d5.root";// too large change */
  /* //TString spectra_file = "./TLee_input_farsideband_EmuCosthetaEhadron/merge_Eproton0d7.root"; */
  /* //TString spectra_file = "./TLee_input_farsideband_EmuCosthetaEhadron/merge_Eproton0d67.root"; */

  /*
  TString spectra_file = "./TLee_input_opendata_EmuCosthetaEhadron/merge.root";
  TString flux_Xs_directory = "./TLee_input_opendata_EmuCosthetaEhadron/flux_Xs/";
  TString detector_directory = "./TLee_input_opendata_EmuCosthetaEhadron/det/";
  TString mc_directory = "./TLee_input_opendata_EmuCosthetaEhadron/mc_stat/";
  */

  /*
  TString spectra_file = "./TLee_input_farsideband_EmuCosthetaEhadron/merge.root";
  TString flux_Xs_directory = "./TLee_input_farsideband_EmuCosthetaEhadron/flux_Xs/";
  TString detector_directory = "./TLee_input_farsideband_EmuCosthetaEhadron/det/";
  TString mc_directory = "./TLee_input_farsideband_EmuCosthetaEhadron/mc_stat/";
  */

  /*
  TString spectra_file = "./TLee_input_EhadronConstraint_FSIbug/merge.root";
  TString flux_Xs_directory = "./TLee_input_EhadronConstraint_FSIbug/flux_Xs/";
  TString detector_directory = "./TLee_input_EhadronConstraint_FSIbug/det/";
  TString mc_directory = "./TLee_input_EhadronConstraint_FSIbug/mc_stat/";
  */  

  /*
  TString spectra_file = "./TLee_input_EhadronConstraint_FSIbug/merge.root";
  TString flux_Xs_directory = "./TLee_input_EhadronConstraint_FSIbug/flux_Xs/";
  TString detector_directory = "./TLee_input_EhadronConstraint_FSIbug/det/";
  TString mc_directory = "./TLee_input_EhadronConstraint_FSIbug/mc_stat/";
  */
  
  TString spectra_file = "/data1/wgu/xs_prl/verify_model_reco_dl_v4/merge.root";
  TString flux_Xs_directory = "/data1/wgu/xs_prl/verify_model_reco_dl_v4/hist_rootfiles/XsFlux/";
  TString detector_directory = "/data1/wgu/xs_prl/verify_model_reco_dl_v4/hist_rootfiles/DetVar/";
  TString mc_directory = "/data1/wgu/xs_prl/verify_model_reco_dl_v4/mc_stat/";
  
  int channels_observation = 6;// data channels (=hdata_obsch_# in spectra_file above)
                               // which is equal to the channels after collapse

  int syst_cov_flux_Xs_begin = 1;// files in flux_Xs_directory above
  int syst_cov_flux_Xs_end   = 17;
 
  int syst_cov_mc_stat_begin = 0;// files in mc_directory above
  int syst_cov_mc_stat_end   = 0;
   

  /// some places may need to be changed when use different file-formats
  /// void TLee::Set_Spectra_MatrixCov()
  /// (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  /// (*) map_Lee_ch                     -----> tag LEE channels
  /// (*) map_detectorfile_str           -----> detector files
 
  
  ////////// display graphics flag

  bool flag_display_graphics = 1;
  
  ////////// systematics flag
  
  bool flag_syst_flux_Xs    = 1;
  bool flag_syst_detector   = 0;
  bool flag_syst_additional = 1;
  bool flag_syst_mc_stat    = 1;

  double Lee_strength_for_outputfile_covariance_matrix = 0;
  
  bool flag_plotting_systematics   = 0;
  
  ////////// goodness of fit
  
  double Lee_strength_for_GoF      = 0;

  bool flag_GoF_output2file_default_0 = 0;
  
  bool flag_both_numuCC            = 0;// 1
  bool flag_CCpi0_FC_by_numuCC     = 0;// 2
  bool flag_CCpi0_PC_by_numuCC     = 0;// 3
  bool flag_NCpi0_by_numuCC        = 0;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = 0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = 0;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = 0;// 7
  bool flag_nueCC_FC_by_all        = 0;// 8

  ////////// Lee strength fitting -- data

  bool flag_Lee_strength_data = 0;

  ////////// MicroBooNE suggested

  bool flag_chi2_data_H0 = 0;
  bool flag_dchi2_H0toH1 = 0;
  
  ////////// Advanced tools
  
  ///// void TLee::Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed)
  ///// do the fitting on the spectra and cov_total after constraint ?
  bool flag_Lee_minimization_after_constraint = 0;// hardcoded, only for the standard 7-ch fitting
  
}
