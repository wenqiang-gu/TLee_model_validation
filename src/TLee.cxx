#include "WCPLEEANA/TLee.h"

#include "draw.icc"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace DataBase {
  double x1[101]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                  21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                  31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                  41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                  51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                  61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
                  71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
                  81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                  91, 92, 93, 94, 95, 96, 97, 98, 99, 100};
  
  double yl[101]={0, 0, 0, 0.856632, 1.70317, 2.51005, 3.32075, 4.14046, 4.9693, 5.80646, 6.65117,
                  7.5025, 8.35978, 9.22237, 10.0898, 10.9615, 11.8372, 12.7165, 13.5992, 14.4849, 15.3734,
                  16.2646, 17.1583, 18.0543, 18.9524, 19.8526, 20.7547, 21.6586, 22.5642, 23.4715, 24.3803,
                  25.2906, 26.2023, 27.1153, 28.0297, 28.9452, 29.8619, 30.7797, 31.6987, 32.6187, 33.5396,
                  34.4616, 35.3845, 36.3083, 37.2329, 38.1584, 39.0847, 40.0118, 40.9396, 41.8682, 42.7975,
                  43.7275, 44.6581, 45.5895, 46.5215, 47.454, 48.3873, 49.321, 50.2554, 51.1903, 52.1257,
                  53.0617, 53.9982, 54.9352, 55.8727, 56.8107, 57.7491, 58.6881, 59.6274, 60.5673, 61.5075,
                  62.4482, 63.3892, 64.3307, 65.2725, 66.2148, 67.1575, 68.1005, 69.0438, 69.9876, 70.9317,
                  71.8761, 72.8209, 73.766, 74.7114, 75.6572, 76.6033, 77.5497, 78.4964, 79.4434, 80.3907,
                  81.3383, 82.2862, 83.2342, 84.1827, 85.1314, 86.0804, 87.0296, 87.9791, 88.9288, 89.8788};

  double yh[101]={1.1478, 2.35971, 3.51917, 4.72422, 5.98186, 7.21064, 8.41858, 9.61053, 10.7896, 11.9582, 13.1179,
                  14.27, 15.4155, 16.5552, 17.6898, 18.8197, 19.9454, 21.0673, 22.1858, 23.3011, 24.4133,
                  25.5229, 26.6299, 27.7346, 28.837, 29.9374, 31.0358, 32.1322, 33.2271, 34.3201, 35.4117,
                  36.5017, 37.5904, 38.6776, 39.7635, 40.8483, 41.9318, 43.0141, 44.0955, 45.1757, 46.2549,
                  47.3331, 48.4104, 49.4868, 50.5623, 51.637, 52.7108, 53.7839, 54.8561, 55.9277, 56.9985,
                  58.0686, 59.1381, 60.2068, 61.275, 62.3425, 63.4094, 64.4757, 65.5415, 66.6066, 67.6713,
                  68.7354, 69.7989, 70.862, 71.9246, 72.9866, 74.0483, 75.1094, 76.1701, 77.2304, 78.2902,
                  79.3496, 80.4085, 81.4672, 82.5253, 83.5831, 84.6406, 85.6976, 86.7542, 87.8105, 88.8665,
                  89.9221, 90.9774, 92.0323, 93.0869, 94.1411, 95.1951, 96.2488, 97.3021, 98.3552, 99.4079,
                  100.46, 101.513, 102.564, 103.616, 104.667, 105.718, 106.769, 107.82, 108.87, 109.92};
}


///////////////////////////////////////////////////////// ccc
void TLee::Exe_Fiedman_Cousins_Data(TMatrixD matrix_fakedata, double Lee_true_low, double Lee_true_hgh, double step)
{
  cout<<endl<<" ---> Exe_Feldman_Cousins_Data"<<endl;
  
  Set_fakedata( matrix_fakedata );
  
  //////////////////
  
  Minimization_Lee_strength_FullCov(1, 0);
  
  cout<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
                        minimization_chi2,
                        minimization_Lee_strength_val,
                        minimization_Lee_strength_err
                        )<<endl;
  
  //////////////////

  double Lee_bestFit_data = minimization_Lee_strength_val;
  double Lee_bestFit_err = minimization_Lee_strength_err;
  double chi2_gmin_data = minimization_chi2;
  vector<int>Lee_scan100_data;
  vector<double>chi2_null_scan_data;
  
  TFile *file_data = new TFile("file_data.root", "recreate");
  TTree *tree_data = new TTree("tree_data", "Feldman-Cousins");

  tree_data->Branch( "Lee_bestFit_data", &Lee_bestFit_data, "Lee_bestFit_data/D" );
  tree_data->Branch( "Lee_bestFit_err", &Lee_bestFit_err, "Lee_bestFit_err/D" );
  tree_data->Branch( "chi2_gmin_data", &chi2_gmin_data, "chi2_gmin_data/D");
  tree_data->Branch( "Lee_scan100_data", &Lee_scan100_data );
  tree_data->Branch( "chi2_null_scan_data", &chi2_null_scan_data );

  int num_scan = int(((Lee_true_hgh-Lee_true_low)/step)+0.5) + 1;
  
  for(int idx=1; idx<=num_scan; idx++ ) {
    if( idx%(max(1, num_scan/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./num_scan, idx)<<endl;

    double Lee_strength = Lee_true_low + (idx-1)*step;
    int Lee_strength_scaled100 = (int)(Lee_strength*100 + 0.5);
    
    Minimization_Lee_strength_FullCov(Lee_strength, 1);
    chi2_null_scan_data.push_back( minimization_chi2 );
    Lee_scan100_data.push_back( Lee_strength_scaled100 );
  }
  cout<<endl;
  
  tree_data->Fill();
    
  file_data->cd();
  tree_data->Write();
  file_data->Close();      
}
  
void TLee::Exe_Fledman_Cousins_Asimov(double Lee_true_low, double Lee_true_hgh, double step)
{
  cout<<endl<<" ---> Exe_Feldman_Cousins_Asimov"<<endl<<endl;
  
  ///////////////////////

  int Lee_strength_scaled100  = 0;  
  vector<double>Lee_scan100;
  vector<double>chi2_null_toy;

  TFile *file_Asimov = new TFile("file_Asimov.root", "recreate");
  TTree *tree_Asimov = new TTree("tree_Asimov", "Feldman-Cousins");

  tree_Asimov->Branch( "Lee_strength_scaled100", &Lee_strength_scaled100, "Lee_strength_scaled100/I" );
  tree_Asimov->Branch( "Lee_scan100", &Lee_scan100 );
  tree_Asimov->Branch( "chi2_null_toy", &chi2_null_toy );

  int num_scan = int(((Lee_true_hgh-Lee_true_low)/step)+0.5) + 1;
  
  for(int idx=1; idx<=num_scan; idx++ ) {
    if( idx%(max(1, num_scan/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./num_scan, idx)<<endl;

    double Lee_strength = Lee_true_low + (idx-1)*step;
    Lee_strength_scaled100 = (int)(Lee_strength*100 + 0.5);
    
    Lee_scan100.clear();
    chi2_null_toy.clear();

    scaleF_Lee = Lee_strength;
    Set_Collapse();
    Set_toy_Asimov();

    for(int jdx=1; jdx<=num_scan; jdx++) {
      double val_Lee_scan = Lee_true_low + (jdx-1)*step;
      double val_Lee_scan100 = (int)(val_Lee_scan*100 + 0.5);
      
      Minimization_Lee_strength_FullCov(val_Lee_scan, 1);
      double val_chi2_null = minimization_chi2;

      Lee_scan100.push_back( val_Lee_scan100 );
      chi2_null_toy.push_back( val_chi2_null );
    }

    tree_Asimov->Fill();
    
  }// idx
  
  file_Asimov->cd();
  tree_Asimov->Write();
  file_Asimov->Close();
  
}

void TLee::Exe_Feldman_Cousins(double Lee_true_low, double Lee_true_hgh, double step, int num_toy, int ifile)
{
  cout<<endl<<" ---> Exe_Feldman_Cousins"<<endl<<endl;

  ///////////////////////

  int Lee_strength_scaled100  = 0;
  vector<double>chi2_null_toy;
  vector<double>chi2_gmin_toy;
  vector<double>LeeF_gmin_toy;

  TFile *file_FC = new TFile(Form("file_FC_%06d.root", ifile), "recreate");
  TTree *tree = new TTree("tree", "Feldman-Cousins");

  tree->Branch( "Lee_strength_scaled100", &Lee_strength_scaled100, "Lee_strength_scaled100/I" );
  tree->Branch( "chi2_null_toy", &chi2_null_toy );
  tree->Branch( "chi2_gmin_toy", &chi2_gmin_toy );
  tree->Branch( "LeeF_gmin_toy", &LeeF_gmin_toy );

  int num_idx = int(((Lee_true_hgh-Lee_true_low)/step)+0.5) + 1;
  
  for(int idx=1; idx<=num_idx; idx++ ) {
    if( idx%(max(1, num_idx/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./num_idx, idx)<<endl;

    double Lee_strength = Lee_true_low + (idx-1)*step;

    /////
    scaleF_Lee = Lee_strength;
    Set_Collapse();
    
    Set_Variations( num_toy );

    Lee_strength_scaled100 = (int)(Lee_strength*100 + 0.5);
    
    chi2_null_toy.clear();
    chi2_gmin_toy.clear();    
    LeeF_gmin_toy.clear();
      
    for(int itoy=1; itoy<=num_toy; itoy++) {
      Set_toy_Variation( itoy );

      Minimization_Lee_strength_FullCov(Lee_strength, 1);
      double val_chi2_null = minimization_chi2;
      
      Minimization_Lee_strength_FullCov(Lee_strength, 0);
      double val_chi2_gmin = minimization_chi2;
      if( minimization_status!=0 ) continue;

      chi2_null_toy.push_back( val_chi2_null );
      chi2_gmin_toy.push_back( val_chi2_gmin );
      LeeF_gmin_toy.push_back( minimization_Lee_strength_val );
      
    }// itoy

    tree->Fill();
    
  }// idx
  
  ///////////////////////
    
  file_FC->cd();
  tree->Write();
  file_FC->Close();
  
}

///////////////////////////////////////////////////////// ccc

void TLee::Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed)
{
  TString roostr = "";
  
  ROOT::Minuit2::Minuit2Minimizer min_Lee( ROOT::Minuit2::kMigrad );
  min_Lee.SetPrintLevel(0);
  min_Lee.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
  min_Lee.SetMaxFunctionCalls(500000);
  min_Lee.SetMaxIterations(500000);
  min_Lee.SetTolerance(1e-6); // tolerance*2e-3 = edm precision
  min_Lee.SetPrecision(1e-18); //precision in the target function

  /// set fitting parameters
  ROOT::Math::Functor Chi2Functor_Lee( [&](const double *par) {// FCN
      TString roostr = "";
      double chi2 = 0;
      double Lee_strength = par[0];

      /////////
      int size_map_fake_data = map_fake_data.size();
      TMatrixD matrix_meas(1, size_map_fake_data);
      for(int ibin=0; ibin<size_map_fake_data; ibin++) {
        matrix_meas(0, ibin) = map_fake_data[ibin];     
      }

      /////////
  
      scaleF_Lee = Lee_strength;
      Set_Collapse();
      
      TMatrixD matrix_pred = matrix_pred_newworld;

      /////////
      TMatrixD matrix_cov_syst = matrix_absolute_cov_newworld;
      
      for(int ibin=0; ibin<matrix_cov_syst.GetNrows(); ibin++) {
        double val_stat_cov = 0;        
        double val_meas = matrix_meas(0, ibin);
        double val_pred = matrix_pred(0, ibin);
        
        if( val_meas==0 ) val_stat_cov = val_pred/2;
        else val_stat_cov = 3./( 1./val_meas + 2./val_pred );   
        if( val_meas==0 && val_pred==0 ) val_stat_cov = 1e-6;
        matrix_cov_syst(ibin, ibin) += val_stat_cov;
      }

      /////////
      TMatrixD matrix_cov_total = matrix_cov_syst;
      TMatrixD matrix_cov_total_inv = matrix_cov_total;
      matrix_cov_total_inv.Invert();
      
      ////////
      TMatrixD matrix_delta = matrix_pred - matrix_meas;
      TMatrixD matrix_delta_T( matrix_delta.GetNcols(), matrix_delta.GetNrows() );
      matrix_delta_T.Transpose( matrix_delta );
      
      TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
      chi2 = matrix_chi2(0,0);      
      
      ///////////////////////////////////////////////////////////////////////////

      if( flag_Lee_minimization_after_constraint ) {// do the fitting on the spectra and cov_total after constraint

        if( 0 ) {
          int num_Y = 26+26;
          int num_X = matrix_cov_syst.GetNrows() - num_Y;
          
          matrix_pred.T(); matrix_meas.T();
          TMatrixD matrix_pred_X = matrix_pred.GetSub(num_Y, num_Y+num_X-1, 0, 0);
          TMatrixD matrix_meas_X = matrix_meas.GetSub(num_Y, num_Y+num_X-1, 0, 0);

          TMatrixD matrix_pred_Y = matrix_pred.GetSub(0, num_Y-1, 0, 0);
          TMatrixD matrix_meas_Y = matrix_meas.GetSub(0, num_Y-1, 0, 0);
          matrix_pred.T(); matrix_meas.T();
        
          TMatrixD matrix_XX = matrix_cov_total.GetSub(num_Y, num_Y+num_X-1, num_Y, num_Y+num_X-1);
          TMatrixD matrix_XX_inv = matrix_XX;
          matrix_XX_inv.Invert();
        
          TMatrixD matrix_YY = matrix_cov_total.GetSub(0, num_Y-1, 0, num_Y-1);
        
          TMatrixD matrix_YX = matrix_cov_total.GetSub(0, num_Y-1, num_Y, num_Y+num_X-1);
          TMatrixD matrix_XY(num_X, num_Y); matrix_XY.Transpose(matrix_YX);
        
          TMatrixD matrix_Y_under_X = matrix_pred_Y + matrix_YX * matrix_XX_inv * (matrix_meas_X - matrix_pred_X);
          TMatrixD matrix_YY_under_XX = matrix_YY - matrix_YX * matrix_XX_inv * matrix_XY;
          
          //////
          
          matrix_Y_under_X.T();
          matrix_meas_Y.T();
          
          TMatrixD matrix_wicons_delta = matrix_Y_under_X - matrix_meas_Y;
          TMatrixD matrix_wicons_delta_T = matrix_wicons_delta.T();
          matrix_wicons_delta.T();

          TMatrixD matrix_YY_under_XX_inv = matrix_YY_under_XX;
          matrix_YY_under_XX_inv.Invert();
          
          TMatrixD matrix_wicons_chi2 = matrix_wicons_delta * matrix_YY_under_XX_inv * matrix_wicons_delta_T;
          double val_wicons_chi2 = matrix_wicons_chi2(0,0);
          chi2 = val_wicons_chi2;       
        }
        
        //////

        if( 1 ) {
          
          vector<int> vc_target_detailed_chs;    
          for(int idx=0; idx<8; idx++) vc_target_detailed_chs.push_back( idx );
          for(int idx=26; idx<26+8; idx++) vc_target_detailed_chs.push_back( idx );       
          vector<int> vc_support_detailed_chs;
          for(int idx=8; idx<26; idx++) vc_support_detailed_chs.push_back( idx );
          for(int idx=26+8; idx<26+26; idx++) vc_support_detailed_chs.push_back( idx );
          for(int idx=26+26; idx<137; idx++) vc_support_detailed_chs.push_back( idx );

          // vector<int> vc_target_detailed_chs;         
          // for(int idx=0; idx<26; idx++) vc_target_detailed_chs.push_back( idx );       
          // vector<int> vc_support_detailed_chs;
          // for(int idx=26; idx<137; idx++) vc_support_detailed_chs.push_back( idx );

          //////
          
          int num_Y = vc_target_detailed_chs.size();
          int num_X = vc_support_detailed_chs.size();
          
          TMatrixD matrix_gof_trans( bins_newworld, num_Y+num_X );// oldworld, newworld
          int new_ch = -1;
          
          for(int idx=0; idx<num_Y; idx++) {
            int old_ch = vc_target_detailed_chs.at(idx);
            new_ch++;
            matrix_gof_trans(old_ch, new_ch) = 1;
          }
          
          for(int idx=0; idx<num_X; idx++) {
            int old_ch = vc_support_detailed_chs.at(idx);
            new_ch++;
            matrix_gof_trans(old_ch, new_ch) = 1;
          }
          
          TMatrixD matrix_gof_trans_T = matrix_gof_trans.T();
          matrix_gof_trans.T();

          TMatrixD matrix_pred_userA = matrix_pred * matrix_gof_trans;
          TMatrixD matrix_meas_userA = matrix_meas * matrix_gof_trans;
          TMatrixD matrix_cov_total_userA = matrix_gof_trans_T * matrix_cov_total * matrix_gof_trans;

          //////
          
          matrix_pred_userA.T(); matrix_meas_userA.T();
          TMatrixD matrix_pred_X = matrix_pred_userA.GetSub(num_Y, num_Y+num_X-1, 0, 0);
          TMatrixD matrix_meas_X = matrix_meas_userA.GetSub(num_Y, num_Y+num_X-1, 0, 0);

          TMatrixD matrix_pred_Y = matrix_pred_userA.GetSub(0, num_Y-1, 0, 0);
          TMatrixD matrix_meas_Y = matrix_meas_userA.GetSub(0, num_Y-1, 0, 0);
          matrix_pred_userA.T(); matrix_meas_userA.T();
        
          TMatrixD matrix_XX = matrix_cov_total_userA.GetSub(num_Y, num_Y+num_X-1, num_Y, num_Y+num_X-1);
          TMatrixD matrix_XX_inv = matrix_XX;
          matrix_XX_inv.Invert();
        
          TMatrixD matrix_YY = matrix_cov_total_userA.GetSub(0, num_Y-1, 0, num_Y-1);
        
          TMatrixD matrix_YX = matrix_cov_total_userA.GetSub(0, num_Y-1, num_Y, num_Y+num_X-1);
          TMatrixD matrix_XY(num_X, num_Y); matrix_XY.Transpose(matrix_YX);
        
          TMatrixD matrix_Y_under_X = matrix_pred_Y + matrix_YX * matrix_XX_inv * (matrix_meas_X - matrix_pred_X);
          TMatrixD matrix_YY_under_XX = matrix_YY - matrix_YX * matrix_XX_inv * matrix_XY;
          
          //////
          
          matrix_Y_under_X.T();
          matrix_meas_Y.T();
          
          TMatrixD matrix_wicons_delta = matrix_Y_under_X - matrix_meas_Y;
          TMatrixD matrix_wicons_delta_T = matrix_wicons_delta.T();
          matrix_wicons_delta.T();

          TMatrixD matrix_YY_under_XX_inv = matrix_YY_under_XX;
          matrix_YY_under_XX_inv.Invert();
          
          TMatrixD matrix_wicons_chi2 = matrix_wicons_delta * matrix_YY_under_XX_inv * matrix_wicons_delta_T;
          double val_wicons_chi2 = matrix_wicons_chi2(0,0);
          chi2 = val_wicons_chi2;
          
        }
        
      }// flag_Lee_minimization_after_constraint
      
      ///////////////////////////////////////////////////////////////////////////      
                  
      return chi2;
      
    },// end of FCN
    1 // number of fitting parameters
    );
  
  min_Lee.SetFunction(Chi2Functor_Lee);
  
  min_Lee.SetVariable( 0, "Lee_strength", Lee_initial_value, 1e-2);
  //min_Lee.SetVariableLowerLimit(0, 0);
  min_Lee.SetLowerLimitedVariable(0, "Lee_strength", Lee_initial_value, 1e-2, 0);
  if( flag_fixed ) {
    min_Lee.SetFixedVariable( 0, "Lee_strength", Lee_initial_value );
  }
  
  /// do the minimization
  min_Lee.Minimize();
  int status_Lee = min_Lee.Status();
  const double *par_Lee = min_Lee.X();
  const double *par_Lee_err = min_Lee.Errors();

  if( status_Lee!=0 ) {
    cerr<<" -----------> Lee strength fitting failed "<<endl;
    minimization_status = status_Lee;
  }

  minimization_status = status_Lee;
  minimization_chi2 = min_Lee.MinValue();
  minimization_Lee_strength_val = par_Lee[0];
  minimization_Lee_strength_err = par_Lee_err[0];

  /// MinosError
  // {
  //   min_Lee.SetErrorDef(1);
  //   double minosError_low = 0;
  //   double minosError_hgh = 0;
  //   min_Lee.GetMinosError(0, minosError_low, minosError_hgh);
  //   cout<<TString::Format(" ---> Best fit of Lee strength: %5.2f, (1 sigma) from MinosError: %5.2f %5.2f",
  //                      minimization_Lee_strength_val, minosError_low, minosError_hgh)<<endl;
  // }
  // {
  //   min_Lee.SetErrorDef(4);
  //   double minosError_low = 0;
  //   double minosError_hgh = 0;
  //   min_Lee.GetMinosError(0, minosError_low, minosError_hgh);
  //   cout<<TString::Format(" ---> Best fit of Lee strength: %5.2f, (2 sigma) from MinosError: %5.2f %5.2f",
  //                      minimization_Lee_strength_val, minosError_low, minosError_hgh)<<endl;
  // }
  // {
  //   min_Lee.SetErrorDef(9);
  //   double minosError_low = 0;
  //   double minosError_hgh = 0;
  //   min_Lee.GetMinosError(0, minosError_low, minosError_hgh);
  //   cout<<TString::Format(" ---> Best fit of Lee strength: %5.2f, (3 sigma) from MinosError: %5.2f %5.2f",
  //                      minimization_Lee_strength_val, minosError_low, minosError_hgh)<<endl;
  // }
  
}  

///////////////////////////////////////////////////////// ccc

void TLee::Set_toy_Asimov()
{
  map_fake_data.clear();
  for(int ibin=0; ibin<bins_newworld; ibin++) map_fake_data[ibin] = matrix_pred_newworld(0, ibin);
}

void TLee::Set_toy_Variation(int itoy)
{
  map_fake_data.clear();
  for(int ibin=0; ibin<bins_newworld; ibin++) map_fake_data[ibin] = map_toy_variation[itoy][ibin];
}

void TLee::Set_measured_data()
{
  map_fake_data.clear();
  for(int ibin=0; ibin<bins_newworld; ibin++) map_fake_data[ibin] = matrix_data_newworld(0, ibin);
}
  
void TLee::Set_fakedata(TMatrixD matrix_fakedata)
{
  map_fake_data.clear();
  int cols = matrix_fakedata.GetNcols();
  for(int ibin=0; ibin<cols; ibin++) map_fake_data[ibin] = matrix_fakedata(0, ibin);    
}

///////////////////////////////////////////////////////// ccc

void TLee::Set_Variations(int num_toy)
{
  map_toy_variation.clear();

  //////////////////////////

  TMatrixDSym DSmatrix_cov(bins_newworld);
  for(int ibin=0; ibin<bins_newworld; ibin++) {
    for(int jbin=0; jbin<bins_newworld; jbin++) {
      DSmatrix_cov(ibin, jbin) = matrix_absolute_cov_newworld(ibin, jbin);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TVectorD matrix_eigenvalue = DSmatrix_eigen.GetEigenValues();

  for(int itoy=1; itoy<=num_toy; itoy++) {    
    TMatrixD matrix_element(bins_newworld, 1);    
    for(int ibin=0; ibin<bins_newworld; ibin++) {
      if( matrix_eigenvalue(ibin)>=0 ) {
        matrix_element(ibin,0) = rand->Gaus( 0, sqrt( matrix_eigenvalue(ibin) ) );
      }
      else {
        matrix_element(ibin,0) = 0;
      }      
    }
    TMatrixD matrix_variation = matrix_eigenvector * matrix_element;
    for(int ibin=0; ibin<bins_newworld; ibin++) {
      double val_with_syst = matrix_variation(ibin,0) + map_pred_spectrum_newworld_bin[ibin];// key point
      if( val_with_syst<0 ) val_with_syst = 0;
      map_toy_variation[itoy][ibin] = rand->PoissonD( val_with_syst );
    }
  }
    
}

///////////////////////////////////////////////////////// ccc

// target detailed channels are constrained by supported detailed channels
int TLee::Exe_Goodness_of_fit_detailed(vector<int>vc_target_detailed_chs, vector<int>vc_support_detailed_chs, int index)
{
  int num_Y = vc_target_detailed_chs.size();
  int num_X = vc_support_detailed_chs.size();

  TMatrixD matrix_gof_trans( bins_newworld, num_Y+num_X );// oldworld, newworld
  int new_ch = -1;
  
  for(int idx=0; idx<num_Y; idx++) {
    int old_ch = vc_target_detailed_chs.at(idx);
    new_ch++;
    matrix_gof_trans(old_ch, new_ch) = 1;
  }
  
  for(int idx=0; idx<num_X; idx++) {
    int old_ch = vc_support_detailed_chs.at(idx);
    new_ch++;
    matrix_gof_trans(old_ch, new_ch) = 1;
  }
  
  TMatrixD matrix_gof_trans_T = matrix_gof_trans.T();
  matrix_gof_trans.T();
  
  TMatrixD matrix_pred = matrix_pred_newworld * matrix_gof_trans;
  TMatrixD matrix_data = matrix_data_newworld * matrix_gof_trans;
  TMatrixD matrix_syst = matrix_gof_trans_T * matrix_absolute_cov_newworld * matrix_gof_trans;

  Exe_Goodness_of_fit(num_Y, num_X, matrix_pred, matrix_data, matrix_syst, index);  
  
  return 1;
}
  
///////////////////////////////////////////////////////// ccc

// target constrained by support
int TLee::Exe_Goodness_of_fit(vector<int>vc_target_chs, vector<int>vc_support_chs, int index)
{
  TString roostr = "";

  //cout<<Form(" ---> Goodness of fit, %2d", index)<<endl;

  /////////////////////////////////////////////////////////////////////////////////////////////

  int size_target = vc_target_chs.size();
  int size_support = vc_support_chs.size();

  for(int idx=0; idx<size_target; idx++) {    
    if( size_support==0 ) break;
    for(int jdx=0; jdx<size_support; jdx++) {      
      if( vc_target_chs.at(idx)==vc_support_chs.at(jdx) ) {
        cout<<endl<<" Are you joking? There is same channel for target and support: "<<vc_target_chs.at(idx)<<endl<<endl;
        exit(1);
      }
    }// jdx
  }// idx

  /////////////////////////

  int num_Y = 0;
  int num_X = 0;
  
  for(int idx=0; idx<size_target; idx++) {
    int ch_target = vc_target_chs.at(idx);
    if(  map_data_spectrum_ch_bin.find(ch_target)==map_data_spectrum_ch_bin.end() ) {
      cout<<" There is no such a target channel: "<<ch_target<<endl;
      exit(1);
    }
    int nbins = map_data_spectrum_ch_bin[ch_target].size();
    num_Y += nbins;
  }
  
  for(int idx=0; idx<size_support; idx++) {
    int ch_support = vc_support_chs.at(idx);
    if(  map_data_spectrum_ch_bin.find(ch_support)==map_data_spectrum_ch_bin.end() ) {
      cout<<" There is no such a support channel: "<<ch_support<<endl;
      exit(1);
    }
    int nbins = map_data_spectrum_ch_bin[ch_support].size();
    num_X += nbins;
  }
  
  /////////////////////////
  
  int gline_eff = -1;  
  map<int, map<int, int> >map_local2global;
  for( auto it_ch=map_data_spectrum_ch_bin.begin(); it_ch!=map_data_spectrum_ch_bin.end(); it_ch++ ) {
    int ich = it_ch->first;
    for( auto it_bin=map_data_spectrum_ch_bin[ich].begin(); it_bin!=map_data_spectrum_ch_bin[ich].end(); it_bin++ ) {
      int ibin = it_bin->first;
      gline_eff++;
      map_local2global[ich][ibin] = gline_eff;
    }
  }
  
  /////////////////////////
  
  TMatrixD matrix_gof_trans( bins_newworld, num_Y+num_X );// oldworld, newworld
  int line_newworld = -1;  
  
  for(int idx=0; idx<size_target; idx++) {
    int ch_target = vc_target_chs.at(idx);  
    int nbins = map_data_spectrum_ch_bin[ch_target].size();
    for(int ibin=0; ibin<nbins; ibin++) {
      int gbin = map_local2global[ch_target][ibin];
      line_newworld++;
      matrix_gof_trans(gbin, line_newworld) = 1;
    }
  }

  for(int idx=0; idx<size_support; idx++) {
    int ch_support = vc_support_chs.at(idx);  
    int nbins = map_data_spectrum_ch_bin[ch_support].size();
    for(int ibin=0; ibin<nbins; ibin++) {
      int gbin = map_local2global[ch_support][ibin];
      line_newworld++;
      matrix_gof_trans(gbin, line_newworld) = 1;
    }
  }

  TMatrixD matrix_gof_trans_T = matrix_gof_trans.T();
  matrix_gof_trans.T();
  
  TMatrixD matrix_pred = matrix_pred_newworld * matrix_gof_trans;
  TMatrixD matrix_data = matrix_data_newworld * matrix_gof_trans;
  TMatrixD matrix_syst = matrix_gof_trans_T * matrix_absolute_cov_newworld * matrix_gof_trans;

  Exe_Goodness_of_fit(num_Y, num_X, matrix_pred, matrix_data, matrix_syst, index);
  
  return 1;
}

///////////////////////////////////////////////////////// ccc

// Y constrained by X

int TLee::Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst, int index)
{
  TString roostr = "";

  if( num_Y+num_X!=matrix_syst.GetNrows() ) {cout<<" ERROR"<<endl; exit(10);}
  
  cout<<Form(" ---> Goodness of fit: %2d", index)<<endl;

  int color_no = kRed;
  int color_wi = kBlue;
  int color_data = kBlack;

  //////////////////////////////////// ttt
  
  bool flag_axis_userAA = 0;
  bool flag_axis_userAB = 0;

  bool flag_left_legend = 0;

  TString str_postfix = "";
  
  ///////////
  
  int axis_user_divisions = 508;
  
  double userAA_index_low = 0;
  double userAA_index_hgh = 16;
  double userAA_value_low = -1;
  double userAA_value_hgh = 1;
  double userAA_TickSize       = 0.06;
  double userAA_LabelSize      = 0.078;  
  double userAA_clone_TickSize = 0.05;
  double userAA_wi2no_TickSize = 0.03;;
  
  double userAB_index_low = 0;
  double userAB_index_hgh = 16;
  double userAB_value_low = -1;
  double userAB_value_hgh = 1;
  double userAB_TickSize       = 0.06;
  double userAB_LabelSize      = 0.078;  
  double userAB_clone_TickSize = 0.05;
  double userAB_wi2no_TickSize = 0.03;;
  
  TString title_axis_user = "E_{#nu}^{rec}";

  TLine *line_FC_PC = new TLine(num_Y/2, 0, num_Y/2, 2);
  line_FC_PC->SetLineColor(kGray+1);
  line_FC_PC->SetLineWidth(4);
  //line_FC_PC->SetLineStyle(7);

  if( index==3001 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 1;
    flag_left_legend = 1;

    str_postfix = "_Emuon";
    title_axis_user = "E^{rec}_{#mu} (MeV)";
    axis_user_divisions = 503;
      
    userAA_index_low = 0;
    userAA_index_hgh = 15;
    userAA_value_low = 100;
    userAA_value_hgh = 1600;
    
    userAB_index_low = 15;
    userAB_index_hgh = 30;
    userAB_value_low = 100;
    userAB_value_hgh = 1600;   
  }

  if( index==30011 || index==30012 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;
    flag_left_legend = 0;

    if( index==30011 ) str_postfix = "_Emuon_FC";
    if( index==30012 ) str_postfix = "_Emuon_PC";
    
    title_axis_user = "E^{rec}_{#mu} (MeV)";
    axis_user_divisions = 505;
      
    userAA_index_low = 0;
    userAA_index_hgh = 15;
    userAA_value_low = 100;
    userAA_value_hgh = 1600;
  }

  if( index==3002 || index==5002 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 1;
    flag_left_legend = 1;

    if( index==3002 ) str_postfix = "_theta";
    if( index==5002 ) str_postfix = "_theta_by_Emuon";    
    
    title_axis_user = "cos#theta^{rec}_{#mu}";
    axis_user_divisions = 505;
      
    userAA_index_low = 0;
    userAA_index_hgh = 20;
    userAA_value_low = -1;
    userAA_value_hgh = 0.9999;
    
    userAB_index_low = 20;
    userAB_index_hgh = 40;
    userAB_value_low = -1;
    userAB_value_hgh = 1;    
  }

  if( index==3003 || index==5001 || index==5003 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 1;
    flag_left_legend = 1;

    if( index==3003 ) str_postfix = "_Ehad";
    if( index==5001 ) str_postfix = "_Ehad_by_Emuon";
    if( index==5003 ) str_postfix = "_Ehad_by_EmuonTheta";
	
    title_axis_user = "E^{rec}_{had} (MeV)";
    axis_user_divisions = 503;
      
    userAA_index_low = 0;
    userAA_index_hgh = 16;
    userAA_value_low = 0;
    userAA_value_hgh = 1600;
    
    userAB_index_low = 16;
    userAB_index_hgh = 32;
    userAB_value_low = 0;
    userAB_value_hgh = 1600;    
  }

  if( index==3004 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 1;
    flag_left_legend = 1;

    str_postfix = "_Ev";
    title_axis_user = "E^{rec}_{#nu} (MeV)";
    axis_user_divisions = 503;
      
    userAA_index_low = 0;
    userAA_index_hgh = 26;
    userAA_value_low = 0;
    userAA_value_hgh = 2600;
    
    userAB_index_low = 26;
    userAB_index_hgh = 52;
    userAB_value_low = 0;
    userAB_value_hgh = 2600;    
  }

  if( index==4001 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;
    flag_left_legend = 0;

    str_postfix = "_Emuon_PCbyFC";
    title_axis_user = "E^{rec}_{#mu} (MeV)";
    axis_user_divisions = 505;
      
    userAA_index_low = 0;
    userAA_index_hgh = 15;
    userAA_value_low = 100;
    userAA_value_hgh = 1600;
  }

  if( index==4002 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;
    flag_left_legend = 1;

    str_postfix = "_theta_PCbyFC";
    title_axis_user = "cos#theta^{rec}_{#mu}";
    axis_user_divisions = 505;
      
    userAA_index_low = 0;
    userAA_index_hgh = 20;
    userAA_value_low = -1;
    userAA_value_hgh = 1;  
  }

  if( index==4003 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;
    flag_left_legend = 0;

    str_postfix = "_Ehad_PCbyFC";
    title_axis_user = "E^{rec}_{had} (MeV)";
    axis_user_divisions = 505;
      
    userAA_index_low = 0;
    userAA_index_hgh = 16;
    userAA_value_low = 0;
    userAA_value_hgh = 1600;    
  }

  if( index==4004 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;
    flag_left_legend = 0;

    str_postfix = "_Ev_PCbyFC";
    title_axis_user = "E^{rec}_{#nu} (MeV)";
    axis_user_divisions = 510;
      
    userAA_index_low = 0;
    userAA_index_hgh = 26;
    userAA_value_low = 0;
    userAA_value_hgh = 2600;    
  }

  if( index==5004 || index==5005) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;
    flag_left_legend = 0;

    if(index==5004) str_postfix = "_Ehad_FC_by_EmuonTheta";
    if(index==5005) str_postfix = "_Ehad_PC_by_EmuonTheta";
    
    title_axis_user = "E^{rec}_{had} (MeV)";
    axis_user_divisions = 505;
      
    userAA_index_low = 0;
    userAA_index_hgh = 16;
    userAA_value_low = 0;
    userAA_value_hgh = 1600;    
  }







 
  
  if( index==1 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 1;

    title_axis_user = "Reco neutrino energy (MeV)";
    axis_user_divisions = 503;
      
    userAA_index_low = 0;
    userAA_index_hgh = 26;
    userAA_value_low = 0;
    userAA_value_hgh = 2600;
    
    userAB_index_low = 26;
    userAB_index_hgh = 52;
    userAB_value_low = 0;
    userAB_value_hgh = 2600;    
  }

  if( index==9 ) {    
    flag_axis_userAA = 1;
    flag_axis_userAB = 1;

    title_axis_user = "Reco neutrino energy (MeV)";
    axis_user_divisions = 504;
      
    userAA_index_low = 0;
    userAA_index_hgh = 31;
    userAA_value_low = 0;
    userAA_value_hgh = 3100;
    
    userAB_index_low = 31;
    userAB_index_hgh = 62;
    userAB_value_low = 0;
    userAB_value_hgh = 3100;    
  }

  if( index==2 || index==3 || index==4 ) {
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;

    title_axis_user = "Reco kinetic energy of #pi^{0} (MeV)";
    axis_user_divisions = 508;
      
    userAA_index_low = 0;
    userAA_index_hgh = 11;
    userAA_value_low = 0;
    userAA_value_hgh = 1100;    
  }
  
  if( index==5 ) {
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;

    title_axis_user = "Reco neutrino energy (MeV)";
    axis_user_divisions = 508;
      
    userAA_index_low = 0;
    userAA_index_hgh = 26;
    userAA_value_low = 0;
    userAA_value_hgh = 2600;    
  }
   
  if( index==6 ) {
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;

    title_axis_user = "Reco neutrino energy (MeV)";
    axis_user_divisions = 508;
      
    userAA_index_low = 0;
    userAA_index_hgh = 18;
    userAA_value_low = 800;
    userAA_value_hgh = 2600;    
  }
  
  if( index==7 ) {
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;

    title_axis_user = "Reco neutrino energy (MeV)";
    axis_user_divisions = 508;
      
    userAA_index_low = 0;
    userAA_index_hgh = 8;
    userAA_value_low = 0;
    userAA_value_hgh = 800;    
  }
    
  if( index==8 ) {
    flag_axis_userAA = 1;
    flag_axis_userAB = 0;

    title_axis_user = "Reco neutrino energy (MeV)";
    axis_user_divisions = 508;
      
    userAA_index_low = 0;
    userAA_index_hgh = 26;
    userAA_value_low = 0;
    userAA_value_hgh = 2600;    
  }
  
  
  ///////////
    
  TGaxis *axis_userAA = new TGaxis(userAA_index_low, 0, userAA_index_hgh, 0,   userAA_value_low, userAA_value_hgh, axis_user_divisions, "S");
  axis_userAA->SetName("axis_userAA");
  axis_userAA->SetTickSize(userAA_TickSize);
  axis_userAA->SetLabelSize(userAA_LabelSize);
  axis_userAA->SetLabelFont(42);
  TGaxis *axis_userAA_clone = (TGaxis*)axis_userAA->Clone("axis_userAA_clone");
  axis_userAA_clone->SetLabelSize(0);
  axis_userAA_clone->SetTickSize(userAA_clone_TickSize);
  TGaxis *axis_userAA_wi2no = (TGaxis*)axis_userAA->Clone("axis_userAA_clone");
  axis_userAA_wi2no->SetLabelSize(userAA_clone_TickSize);
  axis_userAA_wi2no->SetTickSize(userAA_wi2no_TickSize);
    
  TGaxis *axis_userAB = new TGaxis(userAB_index_low, 0, userAB_index_hgh, 0,   userAB_value_low, userAB_value_hgh, axis_user_divisions, "S");
  axis_userAB->SetName("axis_userAB");
  axis_userAB->SetTickSize(userAB_TickSize);
  axis_userAB->SetLabelSize(userAB_LabelSize);
  axis_userAB->SetLabelFont(42);
  TGaxis *axis_userAB_clone = (TGaxis*)axis_userAB->Clone("axis_userAB_clone");
  axis_userAB_clone->SetLabelSize(0);
  axis_userAB_clone->SetTickSize(userAB_clone_TickSize);
  TGaxis *axis_userAB_wi2no = (TGaxis*)axis_userAB->Clone("axis_userAB_clone");
  axis_userAB_wi2no->SetLabelSize(userAB_clone_TickSize);
  axis_userAB_wi2no->SetTickSize(userAB_wi2no_TickSize);
    
  ///////////////////////////////////////////////////////////////////////////////////////////// for no-systematics

  for(int ibin=0; ibin<matrix_syst.GetNrows(); ibin++) {
    if( matrix_syst(ibin, ibin)==0 ) matrix_syst(ibin, ibin) = 1e-6;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////// noConstraint
    
  TMatrixD matrix_cov_stat(num_Y+num_X, num_Y+num_X);  
  TMatrixD matrix_cov_total(num_Y+num_X, num_Y+num_X);
  matrix_cov_total = matrix_cov_stat + matrix_syst;
  // for(int idx=1; idx<=num_Y+num_X; idx++) {
  //   if( matrix_cov_total(idx-1, idx-1)==0 ) matrix_cov_total(idx-1, idx-1) = 1e-6;// case inverse
  // }    
  
  matrix_pred.T(); matrix_data.T();
  TMatrixD matrix_pred_Y = matrix_pred.GetSub(0, num_Y-1, 0, 0);
  TMatrixD matrix_data_Y = matrix_data.GetSub(0, num_Y-1, 0, 0);
  matrix_pred.T(); matrix_data.T();

  TMatrixD matrix_YY = matrix_cov_total.GetSub(0, num_Y-1, 0, num_Y-1);

  {
    TFile *userfile = new TFile("file_user_no.root", "recreate");
    TMatrixD matrix_gof_pred = matrix_pred_Y; matrix_gof_pred.T();
    TMatrixD matrix_gof_meas = matrix_data_Y; matrix_gof_meas.T();
    TMatrixD matrix_gof_syst = matrix_YY;
    matrix_gof_pred.Write("matrix_gof_pred");
    matrix_gof_meas.Write("matrix_gof_meas");
    matrix_gof_syst.Write("matrix_gof_syst");
    userfile->Close();
  }
  
  ///////////////////////////// goodness of fit, Pearson's format

  TMatrixD matrix_goodness_cov_total_noConstraint(num_Y, num_Y);
  // for( int i=0; i<num_Y; i++ ) {
  //   double val_pred = matrix_pred_Y(i, 0);
  //   double val_data = matrix_data_Y(i, 0);    
  //   matrix_goodness_cov_total_noConstraint(i,i) = val_pred;    
  //   if( val_data==1 ) {
  //     if( val_pred<0.461 ) {// DocDB-32520, when the prediction is sufficiently low
  //    double numerator = pow(val_pred-val_data, 2);
  //    double denominator = 2*( val_pred - val_data + val_data*log(val_data/val_pred) );
  //    matrix_goodness_cov_total_noConstraint(i,i) = numerator/denominator;
  //     }
  //   }
  //   if( (val_pred==val_data) && (val_pred==0) ) matrix_goodness_cov_total_noConstraint(i,i) = 1e-6;
  // }  
  matrix_goodness_cov_total_noConstraint = matrix_goodness_cov_total_noConstraint + matrix_YY;

  matrix_pred_Y.T(); matrix_data_Y.T();
  TMatrixD matrix_delta_noConstraint = matrix_pred_Y - matrix_data_Y;
  matrix_pred_Y.T(); matrix_data_Y.T();
  TMatrixD matrix_delta_noConstraint_T = matrix_pred_Y - matrix_data_Y;

  TMatrixD matrix_cov_noConstraint_inv = matrix_goodness_cov_total_noConstraint;
  matrix_cov_noConstraint_inv.Invert();
  TMatrixD matrix_chi2_noConstraint = matrix_delta_noConstraint * matrix_cov_noConstraint_inv * matrix_delta_noConstraint_T;
  double val_chi2_noConstraint = matrix_chi2_noConstraint(0,0);
  double p_value_noConstraint = TMath::Prob( val_chi2_noConstraint, num_Y );
  double val_data_noConstraint = 0;
  double val_pred_noConstraint = 0;
  for(int idx=0; idx<num_Y; idx++) {
    val_data_noConstraint += matrix_data_Y(idx,0);
    val_pred_noConstraint += matrix_pred_Y(idx,0);
  }
  cout<<TString::Format(" ---> GOF noConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, p-value %10.8f, meas/pred %4.2f %4.2f",
                        val_chi2_noConstraint, num_Y, val_chi2_noConstraint/num_Y, p_value_noConstraint,
                        val_data_noConstraint, val_pred_noConstraint
                        )<<endl;
  
  val_GOF_noConstrain = val_chi2_noConstraint;
  val_GOF_NDF = num_Y;

  /////////////////////////////

  roostr = TString::Format("h1_pred_Y_noConstraint_%02d", index);
  TH1D *h1_pred_Y_noConstraint = new TH1D(roostr, "", num_Y, 0, num_Y);
  for(int ibin=1; ibin<=num_Y; ibin++) {
    h1_pred_Y_noConstraint->SetBinContent( ibin, matrix_pred_Y(ibin-1, 0) );
    double val_err = sqrt( matrix_YY(ibin-1, ibin-1) );
    h1_pred_Y_noConstraint->SetBinError( ibin, val_err );
  }

  TGraphAsymmErrors *gh_data = new TGraphAsymmErrors();
  TGraphAsymmErrors *gh_ratio_noConstraint = new TGraphAsymmErrors();
  map<int, double>array_val_data_low;
  map<int, double>array_val_data_hgh;

  double max_data_user = 0;
  double max_pred_user = 0;
  
  for(int ibin=1; ibin<=num_Y; ibin++) {
    double val_data = matrix_data_Y(ibin-1, 0);
    double val_pred_noConstraint = matrix_pred_Y(ibin-1, 0);
    
    double val_data_low = 0;
    double val_data_hgh = 0;
    int idx_data = (int)(val_data+0.5);    
    if( idx_data>100 ) {
      val_data_low = val_data - sqrt(val_data);
      val_data_hgh = val_data + sqrt(val_data);
    }
    else {
      val_data_low = DataBase::yl[ idx_data ];
      val_data_hgh = DataBase::yh[ idx_data ];
    }

    val_data_low = val_data - 1e-6;
    val_data_hgh = val_data + 1e-6;

    if( max_data_user < val_data ) max_data_user = val_data;

    double pred_plus_err = val_pred_noConstraint+h1_pred_Y_noConstraint->GetBinError( ibin );
    if( max_pred_user < pred_plus_err  ) max_pred_user = pred_plus_err;
    
    gh_data->SetPoint( ibin-1, ibin-0.5, val_data );
    gh_data->SetPointError( ibin-1, 0.5, 0.5, val_data-val_data_low, val_data_hgh-val_data );

    array_val_data_low[ibin-1] = val_data_low;
    array_val_data_hgh[ibin-1] = val_data_hgh;
    
    ///////

    double val_ratio_no = val_data/val_pred_noConstraint;
    double val_ratio_no_low = val_ratio_no - val_data_low/val_pred_noConstraint;
    double val_ratio_no_hgh = val_data_hgh/val_pred_noConstraint - val_ratio_no;
    if( val_ratio_no!=val_ratio_no || std::isinf(val_ratio_no) ) val_ratio_no = 0;
    gh_ratio_noConstraint->SetPoint( ibin-1, ibin-0.5, val_ratio_no );
    gh_ratio_noConstraint->SetPointError( ibin-1, 0.5, 0.5, val_ratio_no_low, val_ratio_no_hgh );    
  }
  
  ///////
  
  double ymax_pred = 0;
  double ymax_data = 0;
    
  for(int ibin=1; ibin<=num_Y; ibin++) {
    double val_pred = h1_pred_Y_noConstraint->GetBinContent(ibin)
      + h1_pred_Y_noConstraint->GetBinError(ibin);
    if( ymax_pred<val_pred ) ymax_pred = val_pred;
      
    double xx_data(0), yy_data(0);
    gh_data->GetPoint(ibin-1, xx_data, yy_data);
    if( ymax_data<yy_data ) ymax_data = yy_data;      
  }
  
  roostr = TString::Format("canv_spectra_GoF_no_%02d", index);
  TCanvas *canv_spectra_GoF_no = new TCanvas(roostr, roostr, 1000, 950);
  
  ///////
  canv_spectra_GoF_no->cd();
  TPad *pad_top_no = new TPad("pad_top_no", "pad_top_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_no, 0.15, 0.1, 0.1, 0.05);
  pad_top_no->Draw(); pad_top_no->cd();

  TH1D *h1_pred_Y_noConstraint_clone = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_pred_Y_noConstraint_clone");
  
  h1_pred_Y_noConstraint->Draw("e2");
  h1_pred_Y_noConstraint->SetMinimum(0.);
  //h1_pred_Y_noConstraint->SetMaximum(4.6);
  if( ymax_data>ymax_pred*1.05 ) h1_pred_Y_noConstraint->SetMaximum(ymax_data*1.1);
  //h1_pred_Y_noConstraint->SetMaximum(4500);
  h1_pred_Y_noConstraint->SetMarkerSize(0.);
  h1_pred_Y_noConstraint->SetFillColor(color_no); h1_pred_Y_noConstraint->SetFillStyle(3005);
  h1_pred_Y_noConstraint->SetLineColor(color_no);
  func_title_size(h1_pred_Y_noConstraint, 0.065, 0.065, 0.065, 0.065);
  func_xy_title(h1_pred_Y_noConstraint, "", "Entries");
  h1_pred_Y_noConstraint->GetXaxis()->SetLabelOffset(2);
  h1_pred_Y_noConstraint->GetYaxis()->CenterTitle();
  h1_pred_Y_noConstraint->GetYaxis()->SetTitleOffset(1.2);
  h1_pred_Y_noConstraint->GetYaxis()->SetTickLength(0.02);
    
  h1_pred_Y_noConstraint_clone->Draw("same hist");
  h1_pred_Y_noConstraint_clone->SetLineColor(color_no);
  
  gh_data->Draw("same pe");
  gh_data->SetMarkerStyle(20); gh_data->SetMarkerSize(1.12);
  gh_data->SetMarkerColor(color_data); gh_data->SetLineColor(color_data);
  if( num_X==0 ) {
    gh_data->SetMarkerColor(kBlue); gh_data->SetLineColor(kBlue);
  }

  h1_pred_Y_noConstraint->Draw("same axis");

  TLegend *lg_top_no = new TLegend((index==1)?0.2:0.5, 0.60, (index==1)?0.4:0.85, 0.85);
  if( flag_left_legend ) { lg_top_no->SetX1(0.2); lg_top_no->SetX2(0.4); }
  else { lg_top_no->SetX1(0.5); lg_top_no->SetX2(0.85); }
  if( index==1 || index==7 ) { lg_top_no->SetX1(0.2); lg_top_no->SetX2(0.4);}
  lg_top_no->AddEntry(gh_data, "Data", "lep");
  lg_top_no->AddEntry(h1_pred_Y_noConstraint, TString::Format("#color[%d]{Pred no constraint}", color_no), "lf");
  lg_top_no->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_no, val_chi2_noConstraint, num_Y), "");
  lg_top_no->Draw();
  lg_top_no->SetBorderSize(0); lg_top_no->SetFillStyle(0); lg_top_no->SetTextSize(0.065);

  ///////
  canv_spectra_GoF_no->cd();
  TPad *pad_bot_no = new TPad("pad_bot_no", "pad_bot_no", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_no, 0.15, 0.1, 0.05, 0.3);
  pad_bot_no->Draw(); pad_bot_no->cd();

  TH1D *h1_pred_Y_noConstraint_rel_error = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_pred_Y_noConstraint_rel_error");
  h1_pred_Y_noConstraint_rel_error->Reset();
  for(int ibin=1; ibin<=num_Y; ibin++) {    
    double val_cv = h1_pred_Y_noConstraint->GetBinContent(ibin);
    double val_err = h1_pred_Y_noConstraint->GetBinError(ibin);
    double rel_err = val_err/val_cv;
    if( val_cv==0 ) rel_err = 0;
    h1_pred_Y_noConstraint_rel_error->SetBinContent(ibin, 1);
    h1_pred_Y_noConstraint_rel_error->SetBinError(ibin, rel_err);
  }

  h1_pred_Y_noConstraint_rel_error->Draw("e2");
  h1_pred_Y_noConstraint_rel_error->SetMinimum(0); h1_pred_Y_noConstraint_rel_error->SetMaximum(2);
  func_title_size(h1_pred_Y_noConstraint_rel_error, 0.078, 0.078, 0.078, 0.078);
  func_xy_title(h1_pred_Y_noConstraint_rel_error, "Bin index", "Data / Pred");
  h1_pred_Y_noConstraint_rel_error->GetXaxis()->SetTickLength(0.05);
  h1_pred_Y_noConstraint_rel_error->GetXaxis()->SetLabelOffset(0.005);
  h1_pred_Y_noConstraint_rel_error->GetXaxis()->CenterTitle(); h1_pred_Y_noConstraint_rel_error->GetYaxis()->CenterTitle(); 
  h1_pred_Y_noConstraint_rel_error->GetYaxis()->SetTitleOffset(0.99);
  h1_pred_Y_noConstraint_rel_error->GetYaxis()->SetNdivisions(509);

  TF1 *line_no = new TF1("line_no", "1", 0, 1e6); line_no->Draw("same");
  line_no->SetLineColor(kBlack); line_no->SetLineStyle(7);
  
  gh_ratio_noConstraint->Draw("same pe");
  gh_ratio_noConstraint->SetMarkerStyle(20); gh_ratio_noConstraint->SetMarkerSize(1.12);
  gh_ratio_noConstraint->SetMarkerColor(color_no); gh_ratio_noConstraint->SetLineColor(color_no);  

  h1_pred_Y_noConstraint_rel_error->Draw("same axis");

  if( flag_axis_userAA || flag_axis_userAB ) {
    ///////////////////// bot
    h1_pred_Y_noConstraint_rel_error->GetXaxis()->SetTickLength(0);
    h1_pred_Y_noConstraint_rel_error->GetXaxis()->SetLabelSize(0);
    h1_pred_Y_noConstraint_rel_error->SetXTitle( title_axis_user );
    
    if( flag_axis_userAA && flag_axis_userAB ) line_FC_PC->Draw("same");
    
    if( flag_axis_userAA )  axis_userAA->Draw();
    if( flag_axis_userAB )  axis_userAB->Draw();

    ///////////////////// top
    canv_spectra_GoF_no->cd(); pad_top_no->cd();
    h1_pred_Y_noConstraint->GetXaxis()->SetTickLength(0);
    h1_pred_Y_noConstraint->GetXaxis()->SetLabelSize(0);
    if( flag_axis_userAA ) axis_userAA_clone->Draw();
    if( flag_axis_userAB ) axis_userAB_clone->Draw();    
  }  
  
  if( num_X==0 ) {
    gh_ratio_noConstraint->SetMarkerColor(kBlue); gh_ratio_noConstraint->SetLineColor(kBlue);
    
    roostr = TString::Format("canv_spectra_GoF_no_%02d%s.png", index, str_postfix.Data());
    canv_spectra_GoF_no->SaveAs(roostr);        
    return 1;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////// wiConstraint
    
  matrix_pred.T(); matrix_data.T();
  TMatrixD matrix_pred_X = matrix_pred.GetSub(num_Y, num_Y+num_X-1, 0, 0);
  TMatrixD matrix_data_X = matrix_data.GetSub(num_Y, num_Y+num_X-1, 0, 0);
  matrix_pred.T(); matrix_data.T();

  TMatrixD matrix_XX = matrix_cov_total.GetSub(num_Y, num_Y+num_X-1, num_Y, num_Y+num_X-1);
  //for(int ibin=1; ibin<=num_X; ibin++) matrix_XX(ibin-1, ibin-1) += matrix_pred_X(ibin-1, 0);// Pearson's term for statistics
  TMatrixD matrix_XX_inv = matrix_XX;
  matrix_XX_inv.Invert();

  TMatrixD matrix_YX = matrix_cov_total.GetSub(0, num_Y-1, num_Y, num_Y+num_X-1);
  TMatrixD matrix_XY(num_X, num_Y); matrix_XY.Transpose(matrix_YX);

  TMatrixD matrix_Y_under_X = matrix_pred_Y + matrix_YX * matrix_XX_inv * (matrix_data_X - matrix_pred_X);
  TMatrixD matrix_YY_under_XX = matrix_YY - matrix_YX * matrix_XX_inv * matrix_XY;
  // Here, only for systetmaics uncertainty because of no stat in matrix_YY

  {
    TFile *userfile = new TFile("file_user_wi.root", "recreate");
    TMatrixD matrix_gof_pred = matrix_Y_under_X; matrix_gof_pred.T();
    TMatrixD matrix_gof_meas = matrix_data_Y; matrix_gof_meas.T();
    TMatrixD matrix_gof_syst = matrix_YY_under_XX;
    matrix_gof_pred.Write("matrix_gof_pred");
    matrix_gof_meas.Write("matrix_gof_meas");
    matrix_gof_syst.Write("matrix_gof_syst");
    userfile->Close();
  }
  
  // for(int idx=0; idx<num_Y; idx++) {
  //   cout<<" ---> check "<<idx+1<<"\t"<<matrix_Y_under_X (idx, 0)<<"\t"<<matrix_YY_under_XX(idx, idx)<<endl;
  // }

  
  // for(int idx=0; idx<num_Y; idx++) {
  //   for(int jdx=0; jdx<num_Y; jdx++) {
  //     if( idx!=jdx ) matrix_YY_under_XX(idx, jdx) = 0;
  //   }
  // }
  
  // {
  //   TH2D *h2_corr_after = new TH2D("h2_corr_after", "", num_Y, 0, num_Y, num_Y, 0, num_Y);
    
  //   for( int idx=0; idx<num_Y; idx++ ) {
  //     for( int jdx=0; jdx<num_Y; jdx++ ) {
  // 	double val_cov = matrix_YY_under_XX(idx, jdx);
  // 	double sigma_i = sqrt( matrix_YY_under_XX(idx, idx) );
  // 	double sigma_j = sqrt( matrix_YY_under_XX(jdx, jdx) );

  // 	double val_corr = 0;
  // 	if( sigma_i!=0 && sigma_j!=0 ) val_corr = val_cov/sigma_i/sigma_j;
  // 	if( idx==jdx ) val_corr = 1;
  // 	h2_corr_after->SetBinContent( idx+1, jdx+1, val_corr );
  //     }
  //   }

  //   TCanvas *canv_h2_corr_after= new TCanvas("canv_h2_corr_after", "canv_h2_corr_after", 900, 650);
  //   h2_corr_after->Draw("colz");
  //   canv_h2_corr_after->SaveAs("canv_h2_corr_after.png");
  // }
  
  ///////////////////////////// goodness of fit, Pearson's format

  TMatrixD matrix_goodness_cov_total_wiConstraint(num_Y, num_Y);
  // for( int i=0; i<num_Y; i++ ) {
  //   double val_pred = matrix_Y_under_X(i, 0);
  //   double val_data = matrix_data_Y(i, 0);    
  //   matrix_goodness_cov_total_wiConstraint(i,i) = val_pred;    
  //   if( val_data==1 ) {
  //     if( val_pred<0.461 ) {// DocDB-32520, when the prediction is sufficiently low
  // 	double numerator = pow(val_pred-val_data, 2);
  // 	double dewiminator = 2*( val_pred - val_data + val_data*log(val_data/val_pred) );
  // 	matrix_goodness_cov_total_wiConstraint(i,i) = numerator/dewiminator;
  //     }
  //   }
  //   if( (val_pred==val_data) && (val_pred==0) ) matrix_goodness_cov_total_wiConstraint(i,i) = 1e-6;
  // }  
  matrix_goodness_cov_total_wiConstraint = matrix_goodness_cov_total_wiConstraint + matrix_YY_under_XX;

  matrix_Y_under_X.T(); matrix_data_Y.T();
  //matrix_Y_under_X(0, num_Y-1) = matrix_data_Y(0, num_Y-1);
  TMatrixD matrix_delta_wiConstraint = matrix_Y_under_X - matrix_data_Y;
  matrix_Y_under_X.T(); matrix_data_Y.T();
  TMatrixD matrix_delta_wiConstraint_T = matrix_Y_under_X - matrix_data_Y;  

  TMatrixD matrix_cov_wiConstraint_inv = matrix_goodness_cov_total_wiConstraint;
  matrix_cov_wiConstraint_inv.Invert();
  TMatrixD matrix_chi2_wiConstraint = matrix_delta_wiConstraint * matrix_cov_wiConstraint_inv * matrix_delta_wiConstraint_T;
  double val_chi2_wiConstraint = matrix_chi2_wiConstraint(0,0);
  double p_value_wiConstraint = TMath::Prob( val_chi2_wiConstraint, num_Y );
  double val_data_wiConstraint = 0;
  double val_pred_wiConstraint = 0;
  for(int idx=0; idx<num_Y; idx++) {
    val_data_wiConstraint += matrix_data_Y(idx,0);
    val_pred_wiConstraint += matrix_Y_under_X(idx,0);
  }  
  cout<<TString::Format(" ---> GOF wiConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, p-value %10.8f, meas/pred %4.2f %4.2f",
			val_chi2_wiConstraint, num_Y, val_chi2_wiConstraint/num_Y, p_value_wiConstraint,
			val_data_wiConstraint, val_pred_wiConstraint
			)<<endl<<endl;
  
  val_GOF_wiConstrain = val_chi2_wiConstraint;

  /////////////////////////////

  roostr = TString::Format("h1_pred_Y_wiConstraint_%02d", index);
  TH1D *h1_pred_Y_wiConstraint = (TH1D*)h1_pred_Y_noConstraint->Clone(roostr);
  for(int ibin=1; ibin<=num_Y; ibin++) {
    h1_pred_Y_wiConstraint->SetBinContent( ibin, matrix_Y_under_X(ibin-1, 0) );
    double val_err = sqrt( matrix_YY_under_XX(ibin-1, ibin-1) );
    h1_pred_Y_wiConstraint->SetBinError( ibin, val_err );
  }

  TGraphAsymmErrors *gh_ratio_wiConstraint = new TGraphAsymmErrors();
  
  for(int ibin=1; ibin<=num_Y; ibin++) {
    double val_data = matrix_data_Y(ibin-1, 0);
    double val_pred_wiConstraint = matrix_Y_under_X(ibin-1, 0);
    
    double val_data_low = array_val_data_low[ibin-1];
    double val_data_hgh = array_val_data_hgh[ibin-1];  
    gh_data->SetPoint( ibin-1, ibin-0.5, val_data );
    gh_data->SetPointError( ibin-1, 0.5, 0.5, val_data-val_data_low, val_data_hgh-val_data );

    ///////

    double val_ratio_wi = val_data/val_pred_wiConstraint;
    double val_ratio_wi_low = val_ratio_wi - val_data_low/val_pred_wiConstraint;
    double val_ratio_wi_hgh = val_data_hgh/val_pred_wiConstraint - val_ratio_wi;
    if( val_ratio_wi!=val_ratio_wi || std::isinf(val_ratio_wi) ) val_ratio_wi = 0;
    gh_ratio_wiConstraint->SetPoint( ibin-1, ibin-0.5, val_ratio_wi );
    gh_ratio_wiConstraint->SetPointError( ibin-1, 0.5, 0.5, val_ratio_wi_low, val_ratio_wi_hgh );
  }

  roostr = TString::Format("canv_spectra_GoF_wi_%02d", index);
  TCanvas *canv_spectra_GoF_wi = new TCanvas(roostr, roostr, 1000, 950);
  
  ///////
  canv_spectra_GoF_wi->cd();
  TPad *pad_top_wi = new TPad("pad_top_wi", "pad_top_wi", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_wi, 0.15, 0.1, 0.1, 0.05);
  pad_top_wi->Draw(); pad_top_wi->cd();

  TH1D *h1_pred_Y_wiConstraint_clone = (TH1D*)h1_pred_Y_wiConstraint->Clone("h1_pred_Y_wiConstraint_clone");
  
  h1_pred_Y_wiConstraint->Draw("e2");
  h1_pred_Y_wiConstraint->SetFillColor(color_wi); h1_pred_Y_wiConstraint->SetFillStyle(3004);
  h1_pred_Y_wiConstraint->SetLineColor(color_wi);
    
  h1_pred_Y_wiConstraint_clone->Draw("same hist");
  h1_pred_Y_wiConstraint_clone->SetLineColor(color_wi); h1_pred_Y_wiConstraint_clone->SetFillStyle(0);
  
  gh_data->Draw("same pe");
  
  h1_pred_Y_wiConstraint->Draw("same axis");

  TLegend *lg_top_wi = new TLegend(0.5, 0.60, 0.85, 0.85);
  if( flag_left_legend ) { lg_top_wi->SetX1(0.2); lg_top_wi->SetX2(0.4); }
  else { lg_top_wi->SetX1(0.5); lg_top_wi->SetX2(0.85); }
  if( index==1 || index==7 ) { lg_top_wi->SetX1(0.2); lg_top_wi->SetX2(0.4);}
  lg_top_wi->AddEntry(gh_data, "Data", "lep");
  lg_top_wi->AddEntry(h1_pred_Y_wiConstraint, TString::Format("#color[%d]{Pred w/ constraint}", color_wi), "lf");
  lg_top_wi->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_wi, val_chi2_wiConstraint, num_Y), "");
  lg_top_wi->Draw();
  lg_top_wi->SetBorderSize(0); lg_top_wi->SetFillStyle(0); lg_top_wi->SetTextSize(0.065);

  ///////
  canv_spectra_GoF_wi->cd();
  TPad *pad_bot_wi = new TPad("pad_bot_wi", "pad_bot_wi", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_wi, 0.15, 0.1, 0.05, 0.3);
  pad_bot_wi->Draw(); pad_bot_wi->cd();

  TH1D *h1_pred_Y_wiConstraint_rel_error = (TH1D*)h1_pred_Y_noConstraint_rel_error->Clone("h1_pred_Y_wiConstraint_rel_error");
  h1_pred_Y_wiConstraint_rel_error->Reset();
  for(int ibin=1; ibin<=num_Y; ibin++) {    
    double val_cv = h1_pred_Y_noConstraint->GetBinContent(ibin);
    double val_err = h1_pred_Y_wiConstraint->GetBinError(ibin);
    double rel_err = val_err/val_cv;
    if( val_cv==0 ) rel_err = 0;
    h1_pred_Y_wiConstraint_rel_error->SetBinContent(ibin, 1);
    h1_pred_Y_wiConstraint_rel_error->SetBinError(ibin, rel_err);
  }

  h1_pred_Y_wiConstraint_rel_error->Draw("e2");
  h1_pred_Y_wiConstraint_rel_error->SetFillColor(color_wi); h1_pred_Y_wiConstraint_rel_error->SetFillStyle(3004);
    
  TF1 *line_wi = new TF1("line_wi", "1", 0, 1e6); line_wi->Draw("same");
  line_wi->SetLineColor(kBlack); line_wi->SetLineStyle(7);
  
  gh_ratio_wiConstraint->Draw("same pe");
  gh_ratio_wiConstraint->SetMarkerStyle(20); gh_ratio_wiConstraint->SetMarkerSize(1.12);
  gh_ratio_wiConstraint->SetMarkerColor(color_wi); gh_ratio_wiConstraint->SetLineColor(color_wi);

  h1_pred_Y_wiConstraint_rel_error->Draw("same axis");
   
  if( flag_axis_userAA || flag_axis_userAB ) {
    ///////////////////// bot
    
    if( flag_axis_userAA && flag_axis_userAB ) line_FC_PC->Draw("same");
    
    if( flag_axis_userAA )  axis_userAA->Draw();
    if( flag_axis_userAB )  axis_userAB->Draw();
    
    ///////////////////// top
    canv_spectra_GoF_wi->cd(); pad_top_wi->cd();
    if( flag_axis_userAA ) axis_userAA_clone->Draw();
    if( flag_axis_userAB ) axis_userAB_clone->Draw();    
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////

  canv_spectra_GoF_no->cd(); pad_top_no->cd(); pad_top_no->Update(); double y2_no = gPad->GetUymax();
  canv_spectra_GoF_wi->cd(); pad_top_wi->cd(); pad_top_wi->Update(); double y2_wi = gPad->GetUymax();

  if( y2_no > y2_wi ) {
    canv_spectra_GoF_wi->cd(); pad_top_wi->cd(); h1_pred_Y_wiConstraint->SetMaximum(y2_no);
    h1_pred_Y_wiConstraint->Draw("same e2"); h1_pred_Y_wiConstraint_clone->Draw("same hist");
    gh_data->Draw("same pe"); h1_pred_Y_wiConstraint->Draw("same axis");  
  }
  else {
    canv_spectra_GoF_no->cd(); pad_top_no->cd(); h1_pred_Y_noConstraint->SetMaximum(y2_wi);
    h1_pred_Y_noConstraint->Draw("same e2"); h1_pred_Y_noConstraint_clone->Draw("same hist");
    gh_data->Draw("same pe"); h1_pred_Y_noConstraint->Draw("same axis");  
  }

  roostr = TString::Format("canv_spectra_GoF_no_%02d%s.png", index, str_postfix.Data()); canv_spectra_GoF_no->SaveAs(roostr);
  roostr = TString::Format("canv_spectra_GoF_wi_%02d%s.png", index, str_postfix.Data()); canv_spectra_GoF_wi->SaveAs(roostr);
  
  ////////////////

  roostr = TString::Format("h1_spectra_wi2no_%02d", index);
  TString roostr_wi2no = roostr;
  //TH1D *h1_spectra_wi2no = (TH1D*)h1_pred_Y_noConstraint->Clone(roostr);
  TH1D *h1_spectra_wi2no = new TH1D(roostr, "", num_Y, 0, num_Y);
  h1_spectra_wi2no->SetFillStyle(0);
  for(int ibin=1; ibin<=num_Y; ibin++) {
    double val_wiConstraint = h1_pred_Y_wiConstraint->GetBinContent(ibin);
    double val_noConstraint = h1_pred_Y_noConstraint->GetBinContent(ibin);
    double val_wi2no = val_wiConstraint/val_noConstraint;
    if( val_noConstraint==0 ) val_wi2no = 0;
    h1_spectra_wi2no->SetBinContent(ibin, val_wi2no);
  }
  
  roostr = TString::Format("canv_spectra_wi2no_%02d", index);
  TCanvas *canv_spectra_wi2no = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_spectra_wi2no, 0.15, 0.1, 0.1, 0.15);
  canv_spectra_wi2no->SetGridy();
      
  h1_spectra_wi2no->Draw("hist");
  h1_spectra_wi2no->SetLineColor(color_wi);
  
  h1_spectra_wi2no->SetMinimum(0);
  h1_spectra_wi2no->SetMaximum(2);
  func_title_size(h1_spectra_wi2no, 0.05, 0.05, 0.05, 0.05);
  
  func_xy_title(h1_spectra_wi2no, "Bin index", "Prediction wi/no constraint");
  h1_spectra_wi2no->GetXaxis()->CenterTitle(); h1_spectra_wi2no->GetYaxis()->CenterTitle(); 
  h1_spectra_wi2no->GetYaxis()->SetTitleOffset(1.18);
  h1_spectra_wi2no->GetYaxis()->SetNdivisions(509);
  h1_spectra_wi2no->GetYaxis()->SetTickLength(0.03);    
  if( flag_axis_userAA || flag_axis_userAB ) {/// ttt
    func_xy_title(h1_spectra_wi2no, title_axis_user, "Prediction wi/no constraint");
    h1_spectra_wi2no->GetXaxis()->SetLabelSize(0);
    if( flag_axis_userAA )  axis_userAA_wi2no->Draw();
    if( flag_axis_userAB )  axis_userAB_wi2no->Draw(); 
  }
  
  roostr = TString::Format("canv_spectra_wi2no_%02d%s.png", index, str_postfix.Data()); canv_spectra_wi2no->SaveAs(roostr);

  //h1_spectra_wi2no->SaveAs(roostr_wi2no+".root");

  ////////////////
  
  roostr = TString::Format("canv_spectra_GoF_total_%02d", index);
  TCanvas *canv_spectra_GoF_total = new TCanvas(roostr, roostr, 1000, 950);
  
  ///////
  canv_spectra_GoF_total->cd();
  TPad *pad_top_total = new TPad("pad_top_total", "pad_top_total", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_total, 0.15, 0.1, 0.1, 0.05);
  pad_top_total->Draw(); pad_top_total->cd();

  double ymax_no_total = 0;
  double ymax_data_total = 0;
  for(int ibin=1; ibin<=num_Y; ibin++) {
    double val_no = h1_pred_Y_noConstraint->GetBinContent(ibin)
      + h1_pred_Y_noConstraint->GetBinError(ibin);
    if( ymax_no_total<val_no ) ymax_no_total = val_no;

    double xx_data(0), yy_data(0);
    gh_data->GetPoint(ibin-1, xx_data, yy_data);
    if( ymax_data_total<yy_data ) ymax_data_total = yy_data;
  }
  
  h1_pred_Y_wiConstraint->Draw("e2");

  if( ymax_data_total > ymax_no_total ) {
    h1_pred_Y_wiConstraint->SetMaximum( ymax_data_total*1.1 );
  }
  h1_pred_Y_wiConstraint->SetXTitle("Energy (#times 100 MeV)");
  //if( index==7 ) h1_pred_Y_wiConstraint->SetMaximum(25);
  //if( index==9 ) h1_pred_Y_wiConstraint->SetMaximum(50);
  h1_pred_Y_noConstraint->Draw("same e2");  
  h1_pred_Y_wiConstraint_clone->Draw("same hist");
  h1_pred_Y_noConstraint_clone->Draw("same hist");  
  gh_data->Draw("same pe");
  h1_pred_Y_wiConstraint->Draw("same axis");

  //TFile *file_temp = new TFile( TString::Format("file_temp_%d.root", index), "recreate" );
  //gh_data->SetName("gh_data");
  //gh_data->Write();
  //h1_pred_Y_noConstraint_clone->Write();
  //h1_pred_Y_wiConstraint_clone->Write();  
  //file_temp->Close();
    
  TLegend *lg_top_total = new TLegend(0.5, 0.45, 0.85, 0.85);
  if( flag_left_legend ) { lg_top_total->SetX1(0.2); lg_top_total->SetX2(0.4); }
  else { lg_top_total->SetX1(0.5); lg_top_total->SetX2(0.85); }
  if( index==1 || index==7 ) { lg_top_total->SetX1(0.2); lg_top_total->SetX2(0.4);}
  lg_top_total->AddEntry(gh_data, "Data", "lep");
  lg_top_total->AddEntry(h1_pred_Y_noConstraint, TString::Format("#color[%d]{Pred no constraint}", color_no), "lf");
  lg_top_total->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_no, val_chi2_noConstraint, num_Y), "");
  lg_top_total->AddEntry(h1_pred_Y_wiConstraint, TString::Format("#color[%d]{Pred w/ constraint}", color_wi), "lf");
  lg_top_total->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_wi, val_chi2_wiConstraint, num_Y), "");
  lg_top_total->Draw();
  lg_top_total->SetBorderSize(0); lg_top_total->SetFillStyle(0); lg_top_total->SetTextSize(0.065);

  ///////
  canv_spectra_GoF_total->cd();
  TPad *pad_bot_total = new TPad("pad_bot_total", "pad_bot_total", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_total, 0.15, 0.1, 0.05, 0.3);
  pad_bot_total->Draw(); pad_bot_total->cd();
  //pad_bot_total->SetTicky();

  h1_pred_Y_noConstraint_rel_error->Draw("e2");
  //h1_pred_Y_noConstraint_rel_error->SetYTitle("Uncertainties");
  h1_pred_Y_wiConstraint_rel_error->Draw("same e2");
    
  TF1 *line_total = new TF1("line_total", "1", 0, 1e6); line_total->Draw("same");
  line_total->SetLineColor(kBlack); line_total->SetLineStyle(7);
  
  gh_ratio_noConstraint->Draw("same pe");
  gh_ratio_wiConstraint->Draw("same pe");
      
  if( flag_axis_userAA || flag_axis_userAB ) {
    ///////////////////// bot
    
    if( flag_axis_userAA && flag_axis_userAB ) line_FC_PC->Draw("same");
    
    if( flag_axis_userAA )  axis_userAA->Draw();
    if( flag_axis_userAB )  axis_userAB->Draw();

    ///////////////////// top
    canv_spectra_GoF_total->cd(); pad_top_total->cd();
    if( flag_axis_userAA ) axis_userAA_clone->Draw();
    if( flag_axis_userAB ) axis_userAB_clone->Draw();    
  }
  
  // h1_spectra_wi2no->Draw("same");
  // h1_spectra_wi2no->SetLineColor(kGreen+1);  
  // TLegend *lg_wi2no = new TLegend(0.92, 0.15, 0.94, 0.60);
  // lg_wi2no->SetHeader( TString::Format("#color[%d]{Prediction wi/wo}", kGreen+1) );
  // lg_wi2no->Draw("same"); lg_wi2no->SetTextSize(0.078); lg_wi2no->SetTextAngle(90);
  // lg_wi2no->SetBorderSize(0);
  
  // TLatex *latex = new TLatex(0.5, 0.5, TString::Format("#color[%d]{Predictioin wi/wo}", kGreen+1));
  // latex->Draw("same"); latex->SetTextSize(0.078); //latex->SetTextAngle(90);

  // h1_pred_Y_noConstraint_rel_error->Draw("same axis");  
  // TLegend *lg_bot_total = new TLegend(0.5, 0.85, 0.85, 0.93);
  // lg_bot_total->AddEntry(h1_pred_Y_noConstraint_rel_error, TString::Format("#color[%d]{Prediction wi/no}", kGreen+1), "l");
  // lg_bot_total->Draw();
  // lg_bot_total->SetBorderSize(0); lg_bot_total->SetTextSize(0.078);
  // lg_bot_total->SetFillColor(10); 

  roostr = TString::Format("canv_spectra_GoF_total_%02d%s.png", index, str_postfix.Data()); canv_spectra_GoF_total->SaveAs(roostr);
  //roostr = TString::Format("canv_spectra_GoF_total_%02d.pdf", index); canv_spectra_GoF_total->SaveAs(roostr);
      
  ////////////////

  roostr = TString::Format("h1_spectra_relerr_%02d", index);
  TH1D *h1_spectra_relerr = new TH1D(roostr, "", num_Y, 0, num_Y);

  roostr = TString::Format("h1_spectra_relerr_wi_%02d", index);
  TH1D *h1_spectra_relerr_wi = new TH1D(roostr, "", num_Y, 0, num_Y);

  for(int ibin=1; ibin<=num_Y; ibin++) {    
    double val_noConstraint = h1_pred_Y_noConstraint_rel_error->GetBinError(ibin);
    double val_wiConstraint = h1_pred_Y_wiConstraint_rel_error->GetBinError(ibin);
    h1_spectra_relerr->SetBinContent(ibin, val_noConstraint);
    h1_spectra_relerr_wi->SetBinContent(ibin, val_wiConstraint);    
  }
  
  roostr = TString::Format("canv_spectra_relerr_%02d", index);
  TCanvas *canv_spectra_relerr = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_spectra_relerr, 0.15, 0.1, 0.1, 0.15);
  canv_spectra_relerr->SetGridy();
      
  h1_spectra_relerr->Draw("hist");
  h1_spectra_relerr->SetLineColor(color_no);
       
  h1_spectra_relerr_wi->Draw("same hist");
  h1_spectra_relerr_wi->SetLineColor(color_wi);
  
  h1_spectra_relerr->SetMinimum(0);
  h1_spectra_relerr->SetMaximum(0.5);
  func_title_size(h1_spectra_relerr, 0.05, 0.05, 0.05, 0.05);
  
  func_xy_title(h1_spectra_relerr, "Bin index", "Rel.Err to Pred no constraint");
  h1_spectra_relerr->GetXaxis()->CenterTitle(); h1_spectra_relerr->GetYaxis()->CenterTitle(); 
  h1_spectra_relerr->GetYaxis()->SetTitleOffset(1.18);
  h1_spectra_relerr->GetYaxis()->SetNdivisions(509);
  h1_spectra_relerr->GetYaxis()->SetTickLength(0.03);    
  if( flag_axis_userAA || flag_axis_userAB ) {/// ttt
    func_xy_title(h1_spectra_relerr, title_axis_user,"Rel.Err to Pred no constraint");
    
    if( flag_axis_userAA && flag_axis_userAB ) {
      line_FC_PC->Draw("same");
      line_FC_PC->SetY2(0.5);
    }
    
    h1_spectra_relerr->GetXaxis()->SetLabelSize(0);
    h1_spectra_relerr->GetXaxis()->SetTickSize(0);
    
    if( flag_axis_userAA )  {
      axis_userAA->Draw();
      axis_userAA->SetTickSize(0.06);
      axis_userAA->SetLabelSize(0.05);
    }
    if( flag_axis_userAB )  {
      axis_userAB->Draw();
      axis_userAB->SetTickSize(0.06);
      axis_userAB->SetLabelSize(0.05);
    }
  } 
  
  roostr = TString::Format("canv_spectra_relerr_%02d%s.png", index, str_postfix.Data()); canv_spectra_relerr->SaveAs(roostr);
  //roostr = TString::Format("canv_h1_spectra_relerr_wi_%02d.root", index); h1_spectra_relerr_wi->SaveAs(roostr);
  
  return 1;
}
  
///////////////////////////////////////////////////////// ccc

void TLee::Plotting_systematics()
{
  cout<<" ---> Plotting_systematics"<<endl<<endl;

  int color_flux       = kRed;
  int color_Xs         = kBlue;
  int color_detector   = kMagenta;
  int color_additional = kOrange-3;
  int color_mc_stat    = kGreen+1;
  int color_total      = kBlack;
    
  int rows = bins_newworld;
  int num_ch = map_data_spectrum_ch_bin.size();
  
  if( color_flux+color_Xs+color_detector+color_additional+color_mc_stat+color_total+rows+num_ch==0 ) cout<<" test "<<endl;
  
  ///////////////////////////////////////////

  map<int, double>line_xy;
  map<int, TLine*>line_root_xx;
  map<int, TLine*>line_root_yy;
  
  for(int ich=1; ich<num_ch; ich++) {
    for(int jch=1; jch<=ich; jch++) {
      line_xy[ich] += (int)(map_data_spectrum_ch_bin[jch].size());
    }
    //cout<<Form(" ---> line xy %2d: %4.0f", ich, line_xy[ich])<<endl;    
    line_root_xx[ich] = new TLine( line_xy[ich], 0, line_xy[ich], rows );
    line_root_xx[ich]->SetLineWidth(1); line_root_xx[ich]->SetLineColor(kBlack); line_root_xx[ich]->SetLineStyle(7);
    line_root_yy[ich] = new TLine( 0, line_xy[ich], rows, line_xy[ich]);
    line_root_yy[ich]->SetLineWidth(1); line_root_yy[ich]->SetLineColor(kBlack); line_root_yy[ich]->SetLineStyle(7);
  }

  ///////////////////////////////////////////
  TH2D *h2_covariance_total = new TH2D("h2_covariance_total", "", rows, 0, rows, rows, 0, rows);
  TH2D *h2_correlation_total = new TH2D("h2_correlation_total", "", rows, 0, rows, rows, 0, rows);
  
  TH1D *h1_total_relerr = new TH1D("h1_total_relerr", "", rows, 0, rows);
  TH1D *h1_flux_relerr = new TH1D("h1_flux_relerr", "", rows, 0, rows);
  TH1D *h1_Xs_relerr = new TH1D("h1_Xs_relerr", "", rows, 0, rows);
  TH1D *h1_detector_relerr = new TH1D("h1_detector_relerr", "", rows, 0, rows);
  TH1D *h1_mc_stat_relerr = new TH1D("h1_mc_stat_relerr", "", rows, 0, rows);
  TH1D *h1_additional_relerr = new TH1D("h1_additional_relerr", "", rows, 0, rows);
  
  TH1D *h1_flux_fraction = new TH1D("h1_flux_fraction", "", rows, 0, rows);
  TH1D *h1_Xs_fraction = new TH1D("h1_Xs_fraction", "", rows, 0, rows);
  TH1D *h1_detector_fraction = new TH1D("h1_detector_fraction", "", rows, 0, rows);
  TH1D *h1_mc_stat_fraction = new TH1D("h1_mc_stat_fraction", "", rows, 0, rows);
  TH1D *h1_additional_fraction = new TH1D("h1_additional_fraction", "", rows, 0, rows);
  
  TH1D *h1_pred_totalsyst = new TH1D("h1_pred_totalsyst", "", rows, 0, rows);
  TH1D *h1_meas = new TH1D("h1_meas", "", rows, 0, rows);
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = matrix_absolute_cov_newworld(ibin-1,jbin-1);      
      double cov_i  = matrix_absolute_cov_newworld(ibin-1,ibin-1);
      double cov_j  = matrix_absolute_cov_newworld(jbin-1,jbin-1);      
      
      double val_correlation = cov_ij/sqrt(cov_i*cov_j);
      if( cov_i==0 || cov_j==0 ) val_correlation = 0;
   
      h2_covariance_total->SetBinContent(ibin, jbin, cov_ij);
      h2_correlation_total->SetBinContent(ibin, jbin, val_correlation);

      if( ibin==jbin ) {
	double val_cv = matrix_pred_newworld(0, ibin-1);

	double cov_total      = matrix_absolute_cov_newworld(ibin-1, ibin-1);
	double cov_flux       = matrix_absolute_flux_cov_newworld(ibin-1, ibin-1);
	double cov_Xs         = matrix_absolute_Xs_cov_newworld(ibin-1, ibin-1);
	double cov_detector   = matrix_absolute_detector_cov_newworld(ibin-1, ibin-1);
	double cov_mc_stat    = matrix_absolute_mc_stat_cov_newworld(ibin-1, ibin-1);
	double cov_additional = matrix_absolute_additional_cov_newworld(ibin-1, ibin-1);

	// if( val_cv==0 || val_cv<1e-3) {
	//   cout<<" CV==0 at bin "<<ibin<<" "<<val_cv<<endl;
	// }
	
	if(val_cv!=0) {
	  h1_total_relerr->SetBinContent( ibin, sqrt( cov_total )/val_cv );
	  h1_flux_relerr->SetBinContent( ibin, sqrt( cov_flux )/val_cv );
	  h1_Xs_relerr->SetBinContent( ibin, sqrt(cov_Xs  )/val_cv );
	  h1_detector_relerr->SetBinContent( ibin, sqrt( cov_detector )/val_cv );
	  h1_mc_stat_relerr->SetBinContent( ibin, sqrt( cov_mc_stat )/val_cv );
	  h1_additional_relerr->SetBinContent( ibin, sqrt( cov_additional )/val_cv );
	}

	if( cov_total!=0 ) {
	  h1_flux_fraction->SetBinContent(ibin, cov_flux*100./cov_total );
	  h1_Xs_fraction->SetBinContent(ibin, cov_Xs*100./cov_total );	  
	  h1_detector_fraction->SetBinContent(ibin, cov_detector*100./cov_total );
	  h1_mc_stat_fraction->SetBinContent(ibin, cov_mc_stat*100./cov_total );
	  h1_additional_fraction->SetBinContent(ibin, cov_additional*100./cov_total );
	}
	
	h1_pred_totalsyst->SetBinContent( ibin, val_cv ); h1_pred_totalsyst->SetBinError( ibin, sqrt(cov_total) );
	h1_meas->SetBinContent( ibin, matrix_data_newworld(0, ibin-1) );	  
	
      }// ibin==jbin      
    }// jbin
  }// ibin

  ///////////////////////
  
  TCanvas *canv_h2_correlation_total = new TCanvas("canv_h2_correlation_total", "canv_h2_correlation_total", 800, 700);
  func_canv_margin(canv_h2_correlation_total, 0.15, 0.15, 0.1, 0.15);
  h2_correlation_total->Draw("colz");
  func_title_size(h2_correlation_total, 0.05, 0.05, 0.05, 0.05);
  h2_correlation_total->GetZaxis()->SetLabelSize(0.05);
  h2_correlation_total->GetZaxis()->SetRangeUser(-1, 1);
  func_xy_title(h2_correlation_total, "Bin index", "Bin index");
  h2_correlation_total->GetXaxis()->CenterTitle(); h2_correlation_total->GetYaxis()->CenterTitle();
  h2_correlation_total->GetXaxis()->SetTitleOffset(1.2); h2_correlation_total->GetYaxis()->SetTitleOffset(1.2);
  
  for(int idx=1; idx<num_ch; idx++) {
    line_root_xx[idx]->Draw("same");
    line_root_yy[idx]->Draw("same");
  }
  
  canv_h2_correlation_total->SaveAs("canv_h2_correlation_total.png");
  
  ///////////////////////

  TH2D *h2_relerr_total = new TH2D("h2_relerr_total", "", rows, 0, rows, 100, 0, 1);

  TCanvas *canv_h2_relerr_total = new TCanvas("canv_h2_relerr_total", "canv_h2_relerr_total", 1300, 700);
  func_canv_margin(canv_h2_relerr_total, 0.15, 0.2, 0.11, 0.15);
  h2_relerr_total->Draw();
  func_title_size(h2_relerr_total, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_relerr_total, "Bin index", "Relative error");
  h2_relerr_total->GetXaxis()->CenterTitle(); h2_relerr_total->GetYaxis()->CenterTitle();
  h2_relerr_total->GetXaxis()->SetTitleOffset(1.2); h2_relerr_total->GetYaxis()->SetTitleOffset(1.);
   
  h1_total_relerr->Draw("same hist"); h1_total_relerr->SetLineColor(color_total); h1_total_relerr->SetLineWidth(4);  
  h1_additional_relerr->Draw("same hist"); h1_additional_relerr->SetLineColor(color_additional);  
  h1_mc_stat_relerr->Draw("same hist"); h1_mc_stat_relerr->SetLineColor(color_mc_stat);  
  h1_flux_relerr->Draw("same hist"); h1_flux_relerr->SetLineColor(color_flux);  
  h1_Xs_relerr->Draw("same hist"); h1_Xs_relerr->SetLineColor(color_Xs);  
  h1_detector_relerr->Draw("same hist"); h1_detector_relerr->SetLineColor(color_detector);
  
  for(int idx=1; idx<num_ch; idx++) {
    line_root_xx[idx]->Draw(); line_root_xx[idx]->SetLineStyle(7); line_root_xx[idx]->SetY2(1);
  }

  TLegend *lg_relerr_total = new TLegend(0.82, 0.5, 0.95, 0.89);
  lg_relerr_total->AddEntry(h1_total_relerr, "Total", "l");
  lg_relerr_total->AddEntry(h1_flux_relerr, "Flux", "l");
  lg_relerr_total->AddEntry(h1_Xs_relerr, "Xs", "l");
  lg_relerr_total->AddEntry(h1_detector_relerr, "Detector", "l");
  lg_relerr_total->AddEntry(h1_mc_stat_relerr, "Stat", "l");
  //lg_relerr_total->AddEntry(h1_additional_relerr, "Dirt", "l");
  lg_relerr_total->Draw();
  lg_relerr_total->SetTextSize(0.05);
    
  h2_relerr_total->Draw("same axis");
  canv_h2_relerr_total->SaveAs("canv_h2_relerr_total.png");

  /////////////////////////
  
  THStack *h1_stack_fraction = new THStack("h1_stack_fraction", "");
  h1_stack_fraction->Add(h1_flux_fraction);
  h1_flux_fraction->SetFillColor(color_flux); h1_flux_fraction->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_Xs_fraction);
  h1_Xs_fraction->SetFillColor(color_Xs); h1_Xs_fraction->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_detector_fraction);
  h1_detector_fraction->SetFillColor(color_detector); h1_detector_fraction->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_mc_stat_fraction);
  h1_mc_stat_fraction->SetFillColor(color_mc_stat); h1_mc_stat_fraction->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_additional_fraction);
  h1_additional_fraction->SetFillColor(color_additional); h1_additional_fraction->SetLineColor(kBlack);
    
  TH2D *h2_basic_fraction = new TH2D("h2_basic_fraction", "", rows, 0, rows, 110, 0, 110);

  TCanvas *canv_h2_basic_fraction = new TCanvas("canv_h2_basic_fraction", "canv_h2_basic_fraction", 1300, 700);
  func_canv_margin(canv_h2_basic_fraction, 0.15, 0.2, 0.11, 0.15);
  h2_basic_fraction->Draw();
  func_title_size(h2_basic_fraction, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h2_basic_fraction, "Bin index", "Syst. percentage");
  h2_basic_fraction->GetXaxis()->CenterTitle(); h2_basic_fraction->GetYaxis()->CenterTitle();
  h2_basic_fraction->GetXaxis()->SetTitleOffset(1.2); h2_basic_fraction->GetYaxis()->SetTitleOffset(1.05);
  
  h1_stack_fraction->Draw("same");
  
  for(int idx=1; idx<num_ch; idx++) {
    line_root_xx[idx]->Draw(); line_root_xx[idx]->SetLineStyle(7); line_root_xx[idx]->SetY2(110);
  }

  TLegend *lg_fraction_total = new TLegend(0.82, 0.55, 0.95, 0.89);
  lg_fraction_total->AddEntry(h1_flux_fraction, "Flux", "f");
  lg_fraction_total->AddEntry(h1_Xs_fraction, "Xs", "f");  
  lg_fraction_total->AddEntry(h1_detector_fraction, "Detector", "f");
  lg_fraction_total->AddEntry(h1_mc_stat_fraction, "Stat", "f");
  //lg_fraction_total->AddEntry(h1_additional_fraction, "Dirt", "f");
  lg_fraction_total->Draw();
  lg_fraction_total->SetTextSize(0.05);
      
  h2_basic_fraction->Draw("same axis");
  canv_h2_basic_fraction->SaveAs("canv_h2_basic_fraction.png");
  
  /////////////////////////

  TCanvas *canv_h1_pred_totalsyst = new TCanvas("canv_h1_pred_totalsyst", "canv_h1_pred_totalsyst", 1300, 700);
  func_canv_margin(canv_h1_pred_totalsyst, 0.15, 0.1, 0.11, 0.15);
  canv_h1_pred_totalsyst->SetLogy();
  
  TH1D *h1_pred_totalsyst_clone = (TH1D*)h1_pred_totalsyst->Clone("h1_pred_totalsyst_clone");
  
  h1_pred_totalsyst->Draw("e2"); h1_pred_totalsyst->SetFillColor(kRed); h1_pred_totalsyst->SetMarkerSize(0);
  h1_pred_totalsyst->SetMinimum(1e-3);
  func_title_size(h1_pred_totalsyst, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_pred_totalsyst, "Bin index", "Entries");
  h1_pred_totalsyst->GetXaxis()->CenterTitle(); h1_pred_totalsyst->GetYaxis()->CenterTitle();
  h1_pred_totalsyst->GetXaxis()->SetTitleOffset(1.2); h1_pred_totalsyst->GetYaxis()->SetTitleOffset(1.05);  
  
  h1_pred_totalsyst_clone->Draw("same hist"); h1_pred_totalsyst_clone->SetLineColor(kBlack);

  canv_h1_pred_totalsyst->cd(); canv_h1_pred_totalsyst->Update(); double ymax_canv_h1_pred_totalsyst = gPad->GetUymax();
  for(int idx=1; idx<num_ch; idx++) {
    line_root_xx[idx]->Draw(); line_root_xx[idx]->SetLineStyle(7);
    line_root_xx[idx]->SetY2( pow(10, ymax_canv_h1_pred_totalsyst) );
    line_root_xx[idx]->SetY1(1e-3);
  }

  canv_h1_pred_totalsyst->SaveAs("canv_h1_pred_totalsyst.png");
  
  
    
}

///////////////////////////////////////////////////////// ccc

void TLee::Set_Collapse()
{
  //////////////////////////////////////// pred

  TMatrixD matrix_transform_Lee = matrix_transform;
  for(int ibin=0; ibin<matrix_transform_Lee.GetNrows(); ibin++) {
    for(int jbin=0; jbin<matrix_transform_Lee.GetNcols(); jbin++) {
      if( map_Lee_oldworld.find(ibin)!=map_Lee_oldworld.end() ) matrix_transform_Lee(ibin, jbin) *= scaleF_Lee;
    }
  }
  
  map_pred_spectrum_newworld_bin.clear();
  TMatrixD matrix_pred_oldworld(1, bins_oldworld);
  for(int ibin=0; ibin<bins_oldworld; ibin++) matrix_pred_oldworld(0, ibin) = map_input_spectrum_oldworld_bin[ibin];

  matrix_pred_newworld.Clear();
  matrix_pred_newworld.ResizeTo(1, bins_newworld);
  matrix_pred_newworld = matrix_pred_oldworld * matrix_transform_Lee;

  //matrix_pred_newworld(0, bins_newworld-1) = matrix_data_newworld(0, bins_newworld-1);
  
  if( bins_newworld!=matrix_pred_newworld.GetNcols() ) { cerr<<"bins_newworld!=matrix_pred_newworld.GetNcols()"<<endl; exit(1); }
  for(int ibin=0; ibin<bins_newworld; ibin++) map_pred_spectrum_newworld_bin[ibin] = matrix_pred_newworld(0, ibin);

  
  ////////////////////////////////////////
  
  matrix_absolute_cov_oldworld.Clear();
  matrix_absolute_cov_oldworld.ResizeTo( bins_oldworld, bins_oldworld );

  TMatrixD matrix_transform_Lee_T( bins_newworld, bins_oldworld );
  matrix_transform_Lee_T.Transpose( matrix_transform_Lee );

  matrix_absolute_cov_newworld.Clear();
  matrix_absolute_cov_newworld.ResizeTo(bins_newworld, bins_newworld);

  TMatrixD matrix_cov_flux_Xs_newworld = matrix_transform_Lee_T * matrix_input_cov_flux_Xs * matrix_transform_Lee;
  TMatrixD matrix_cov_detector_newworld = matrix_transform_Lee_T * matrix_input_cov_detector * matrix_transform_Lee;
  TMatrixD matrix_cov_additional_newworld = matrix_transform_Lee_T * matrix_input_cov_additional * matrix_transform_Lee;

  TMatrixD matrix_absolute_data_stat_cov(bins_newworld, bins_newworld);
  TMatrixD matrix_absolute_pred_stat_cov(bins_newworld, bins_newworld);
  
  if( flag_syst_mc_stat ) {
    //for(int ibin=0; ibin<bins_newworld; ibin++) {
    //double val_mc_stat_cov = gh_mc_stat_bin[ibin]->Eval( scaleF_Lee );
    //matrix_absolute_cov_newworld(ibin, ibin) += val_mc_stat_cov;
    //matrix_absolute_cov_newworld(ibin, ibin) += val_mc_stat_cov/4.;      
    //}

    for(int idx=0; idx<bins_newworld; idx++) {
      for(int jdx=0; jdx<bins_newworld; jdx++) {

	/// Data_stat
	double data_sigma_i = sqrt( matrix_pred_newworld(0, idx) );
	double data_sigma_j = sqrt( matrix_pred_newworld(0, jdx) );
	double data_correlation = matrix_data_MCstat_correlation(idx, jdx);
	double data_stat_cov = data_correlation * data_sigma_i * data_sigma_j;

	/// MC_stat
	double pred_sigma_i = sqrt( gh_mc_stat_bin[idx]->Eval( scaleF_Lee ) );
	double pred_sigma_j = sqrt( gh_mc_stat_bin[jdx]->Eval( scaleF_Lee ) );		
	double pred_correlation = matrix_pred_MCstat_correlation(idx, jdx);
	double pred_stat_cov = pred_correlation * pred_sigma_i * pred_sigma_j;

	///
	//matrix_absolute_cov_newworld(idx, jdx) += data_stat_cov;
	//matrix_absolute_cov_newworld(idx, jdx) += pred_stat_cov;
	matrix_absolute_data_stat_cov(idx, jdx) = data_stat_cov;
	matrix_absolute_pred_stat_cov(idx, jdx) = pred_stat_cov;
      }
    }



    
    if( 0 ) {
      //matrix_absolute_data_stat_cov.SaveAs("matrix_absolute_data_stat_cov.root");
      //matrix_absolute_pred_stat_cov.SaveAs("matrix_absolute_pred_stat_cov.root");

      cout<<" -----------------> test "<<endl;
      cout<<endl;

      TMatrixD matrix_trans_user(bins_newworld, 15*2 + 20*2 + 16*2);
      for(int idx=1; idx<=30; idx++) matrix_trans_user(idx-1, idx-1) = 1;
      for(int idx=1; idx<=20; idx++) matrix_trans_user(15*2+idx-1, 30+idx-1) = 1;
      for(int idx=1; idx<=20; idx++) matrix_trans_user(15*2+21+idx-1, 30+20+idx-1) = 1;
      for(int idx=1; idx<=32; idx++) matrix_trans_user(15*2+21*2+idx-1, 30+20*2+idx-1) = 1;

      TMatrixD matrix_trans_user_T = matrix_trans_user.T(); matrix_trans_user.T();

      TMatrixD matrix_user_stat_cov = matrix_trans_user_T * matrix_absolute_pred_stat_cov * matrix_trans_user;

      int rows_user = matrix_user_stat_cov.GetNrows();

      TH2D *h2_user_stat_cov = new TH2D("h2_user_stat_cov", "", rows_user, 0, rows_user, rows_user, 0, rows_user);
      for(int idx=0; idx<rows_user; idx++) {
	for(int jdx=0; jdx<rows_user; jdx++) {
	  double cov = matrix_user_stat_cov(idx, jdx);
	  h2_user_stat_cov->SetBinContent(idx+1, jdx+1, cov);
	}
      }
      
      TCanvas *canv_h2_user_stat_cov = new TCanvas("canv_h2_user_stat_cov", "canv_h2_user_stat_cov", 800, 700);
      func_canv_margin(canv_h2_user_stat_cov, 0.15, 0.15, 0.1, 0.15);
      canv_h2_user_stat_cov->SetLogz();
      h2_user_stat_cov->Draw("colz");
      func_title_size(h2_user_stat_cov, 0.05, 0.05, 0.05, 0.05);
      h2_user_stat_cov->GetZaxis()->SetLabelSize(0.05);
      //h2_user_stat_cov->GetZaxis()->SetRangeUser(-1, 1);
      func_xy_title(h2_user_stat_cov, "Bin index", "Bin index");
      h2_user_stat_cov->GetXaxis()->CenterTitle(); h2_user_stat_cov->GetYaxis()->CenterTitle();
      h2_user_stat_cov->GetXaxis()->SetTitleOffset(1.2); h2_user_stat_cov->GetYaxis()->SetTitleOffset(1.2);
    }


    
  }
  
  /*
  for(int idx=1; idx<=16; idx++) {
      int ibin = 15*2 + 21*2 + idx;
      int jbin = 15*2 + 21*2 + 16 + idx;
      
      for(int ii=1; ii<=bins_newworld; ii++) {
	
	matrix_cov_flux_Xs_newworld( ibin-1, ii-1 ) = 0; matrix_cov_flux_Xs_newworld( ii-1, ibin-1 ) = 0;
	matrix_cov_flux_Xs_newworld( jbin-1, ii-1 ) = 0; matrix_cov_flux_Xs_newworld( ii-1, jbin-1 ) = 0;
	
	matrix_cov_detector_newworld( ibin-1, ii-1 ) = 0; matrix_cov_detector_newworld( ii-1, ibin-1 ) = 0;
	matrix_cov_detector_newworld( jbin-1, ii-1 ) = 0; matrix_cov_detector_newworld( ii-1, jbin-1 ) = 0;
	
	matrix_cov_additional_newworld( ibin-1, ii-1 ) = 0; matrix_cov_additional_newworld( ii-1, ibin-1 ) = 0;
	matrix_cov_additional_newworld( jbin-1, ii-1 ) = 0; matrix_cov_additional_newworld( ii-1, jbin-1 ) = 0;

	matrix_absolute_data_stat_cov( ibin-1, ii-1 ) = 0; matrix_absolute_data_stat_cov( ii-1, ibin-1 ) = 0;
	matrix_absolute_data_stat_cov( jbin-1, ii-1 ) = 0; matrix_absolute_data_stat_cov( ii-1, jbin-1 ) = 0;
	matrix_absolute_pred_stat_cov( ibin-1, ii-1 ) = 0; matrix_absolute_pred_stat_cov( ii-1, ibin-1 ) = 0;
	matrix_absolute_pred_stat_cov( jbin-1, ii-1 ) = 0; matrix_absolute_pred_stat_cov( ii-1, jbin-1 ) = 0;
	
      }
  }
  */
  
  if( flag_syst_flux_Xs ) matrix_absolute_cov_newworld += matrix_cov_flux_Xs_newworld;
  if( flag_syst_detector ) matrix_absolute_cov_newworld += matrix_cov_detector_newworld;
  if( flag_syst_additional ) matrix_absolute_cov_newworld += matrix_cov_additional_newworld;
  if( flag_syst_mc_stat ) {
    matrix_absolute_cov_newworld += matrix_absolute_data_stat_cov;
    matrix_absolute_cov_newworld += matrix_absolute_pred_stat_cov;    
  }

  //matrix_data_MCstat_correlation
  
  ////////////////////////////////////////

  if( flag_individual_cov_newworld ) {
    cout<<" ---> Producing the systematics for plotting (should appear only one time)"<<endl;
    cout<<" ---> The LEE strength used for the producing is corresponding to the one in the Configure_LEE.h"<<endl<<endl;
    
    flag_individual_cov_newworld = false;

    matrix_absolute_flux_cov_newworld.Clear();
    matrix_absolute_Xs_cov_newworld.Clear();
    matrix_absolute_detector_cov_newworld.Clear();
    matrix_absolute_mc_stat_cov_newworld.Clear();
    matrix_absolute_additional_cov_newworld.Clear();
    
    matrix_absolute_flux_cov_newworld.ResizeTo( bins_newworld, bins_newworld );
    matrix_absolute_Xs_cov_newworld.ResizeTo( bins_newworld, bins_newworld );
    matrix_absolute_detector_cov_newworld.ResizeTo( bins_newworld, bins_newworld );
    matrix_absolute_mc_stat_cov_newworld.ResizeTo( bins_newworld, bins_newworld );
    matrix_absolute_additional_cov_newworld.ResizeTo( bins_newworld, bins_newworld );
                
    for(auto it=matrix_input_cov_detector_sub.begin(); it!=matrix_input_cov_detector_sub.end(); it++) {
      int idx = it->first;
      matrix_absolute_detector_sub_cov_newworld[idx].Clear();
      matrix_absolute_detector_sub_cov_newworld[idx].ResizeTo( bins_newworld, bins_newworld );
      matrix_absolute_detector_sub_cov_newworld[idx] = matrix_transform_Lee_T * matrix_input_cov_detector_sub[idx] * matrix_transform_Lee;
    }
      
    matrix_absolute_flux_cov_newworld = matrix_transform_Lee_T * matrix_input_cov_flux * matrix_transform_Lee;
    matrix_absolute_Xs_cov_newworld = matrix_transform_Lee_T * matrix_input_cov_Xs * matrix_transform_Lee;
    matrix_absolute_detector_cov_newworld = matrix_transform_Lee_T * matrix_input_cov_detector * matrix_transform_Lee;
    matrix_absolute_additional_cov_newworld = matrix_transform_Lee_T * matrix_input_cov_additional * matrix_transform_Lee;

    matrix_absolute_mc_stat_cov_newworld += (matrix_absolute_data_stat_cov + matrix_absolute_pred_stat_cov);
    
    // for(int ibin=0; ibin<bins_newworld; ibin++) {
    //   double val_mc_stat_cov = gh_mc_stat_bin[ibin]->Eval( scaleF_Lee );
    //   matrix_absolute_mc_stat_cov_newworld(ibin, ibin) = val_mc_stat_cov;
    // }// ibin    
  }
  
}

///////////////////////////////////////////////////////// ccc

void TLee::Set_TransformMatrix()
{
   cout<<endl<<" ---> Set_TransformMatrix"<<endl<<endl;

   ////////////////////////////// correponding to "Set_Spectra_MatrixCov"

}

///////////////////////////////////////////////////////// ccc

void TLee::Set_POT_implement()
{
  cout<<endl<<" ---> Set_POT_implement"<<endl;
  
  ////////////////////////////// pred

  int line_pred = -1;
  for( auto it_ch=map_input_spectrum_ch_bin.begin(); it_ch!=map_input_spectrum_ch_bin.end(); it_ch++ ) {
    int ich = it_ch->first;
    for(int ibin=0; ibin<(int)map_input_spectrum_ch_bin[ich].size(); ibin++) {
      line_pred++;
      map_input_spectrum_ch_bin[ich][ibin] *= scaleF_POT;
      map_input_spectrum_oldworld_bin[line_pred] *= scaleF_POT;
    }// ibin
  }// ich
  
  ////////////////////////////// data

  int line_data = -1;
  for( auto it_ch=map_data_spectrum_ch_bin.begin(); it_ch!=map_data_spectrum_ch_bin.end(); it_ch++ ) {
    int ich = it_ch->first;
    for( int ibin=0; ibin<(int)map_data_spectrum_ch_bin[ich].size(); ibin++ ) {
      line_data++;
      map_data_spectrum_ch_bin[ich][ibin] *= scaleF_POT;
      map_data_spectrum_newworld_bin[line_data] *= scaleF_POT;
    }// ibin
  }// ich

  matrix_data_newworld.Clear();
  matrix_data_newworld.ResizeTo(1, bins_newworld);
  for(int ibin=0; ibin<bins_newworld; ibin++) {
    matrix_data_newworld(0, ibin) = map_data_spectrum_newworld_bin[ibin];
  }
    
  ////////////////////////////// flux_Xs, detector, additional, mc_stat

  double scaleF_POT2 = scaleF_POT * scaleF_POT;
  
  for(int ibin=0; ibin<bins_oldworld; ibin++) {
    for(int jbin=0; jbin<bins_oldworld; jbin++) {      
      matrix_input_cov_flux_Xs(ibin, jbin) *= scaleF_POT2;
      matrix_input_cov_flux(ibin, jbin) *= scaleF_POT2;
      matrix_input_cov_Xs(ibin, jbin) *= scaleF_POT2;      
      matrix_input_cov_detector(ibin, jbin) *= scaleF_POT2;
      matrix_input_cov_additional(ibin, jbin) *= scaleF_POT2;
            
      for(auto it=matrix_input_cov_detector_sub.begin(); it!=matrix_input_cov_detector_sub.end(); it++) {
	int idx = it->first;
	matrix_input_cov_detector_sub[idx](ibin, jbin) *= scaleF_POT2;
      }
  
    }// jbin
  }// ibin

  for(auto it=gh_mc_stat_bin.begin(); it!=gh_mc_stat_bin.end(); it++) {
    int ibin = it->first; //cout<<Form(" ---> check %3d, %3d", ibin, gh_mc_stat_bin[ibin]->GetN())<<endl;
    for(int idx=0; idx<gh_mc_stat_bin[ibin]->GetN(); idx++) {
      double x(0), y(0);
      gh_mc_stat_bin[ibin]->GetPoint(idx, x, y);
      gh_mc_stat_bin[ibin]->SetPoint(idx, x, y*scaleF_POT2);
    }// ipoint
  }// ibin
  
}

///////////////////////////////////////////////////////// ccc

void TLee::Set_config_file_directory(TString spectra_file_, TString flux_Xs_directory_, TString detector_directory_, TString mc_directory_)
{
  cout<<endl<<" ---> Set_config_file_directory"<<endl<<endl;

  spectra_file       = spectra_file_;
  flux_Xs_directory  = flux_Xs_directory_;
  detector_directory = detector_directory_;
  mc_directory       = mc_directory_;

  cout<<Form(" spectra_file       %-10s", spectra_file.Data() )<<endl;
  cout<<Form(" flux_Xs_directory  %-10s", flux_Xs_directory.Data() )<<endl;
  cout<<Form(" detector_directory %-10s", detector_directory.Data() )<<endl;
  cout<<Form(" mc_directory       %-10s", mc_directory.Data() )<<endl;  
}


void TLee::Set_Spectra_MatrixCov()
{
  /// spectra should be consist with matrix-cov order
  
  cout<<endl<<" ---> Set_Spectra_MatrixCov"<<endl<<endl;
  TString roostr = "";
    
  //////////////////////////////////////
  /*
  map_input_spectrum_ch_str[1] = "FC_Emuon";
  map_input_spectrum_ch_str[2] = "PC_Emuon";
  map_input_spectrum_ch_str[3] = "FC_theta";
  map_input_spectrum_ch_str[4] = "PC_theta"; 
  map_input_spectrum_ch_str[5] = "FC_Ehadron"; 
  map_input_spectrum_ch_str[6] = "PC_Ehadron"; 
  */
  
  //for(int idx=1; idx<=16; idx++) map_input_spectrum_ch_str[idx] = TString::Format("pred_%02d", idx);
  //for(int idx=1; idx<=6; idx++) map_input_spectrum_ch_str[idx] = TString::Format("pred_%02d", idx);
  
  for(int idx=1; idx<=12; idx++) map_input_spectrum_ch_str[idx] = TString::Format("pred_%02d", idx);
  // for(int idx=1; idx<=32; idx++) map_input_spectrum_ch_str[idx] = TString::Format("pred_%02d", idx);
  // map_Lee_ch[14] = 1;
  // map_Lee_ch[15] = 1;
  // map_Lee_ch[16] = 1;
  // map_Lee_ch[17] = 1;
  // map_Lee_ch[18] = 1;
  // map_Lee_ch[19] = 1;
 
  ////////////////////////////////////// pred
  /*
  map_input_spectrum_ch_str[1] = "nueCC_FC_norm";
  map_input_spectrum_ch_str[2] = "nueCC_PC_norm";
  map_input_spectrum_ch_str[3] = "numuCC_FC_norm";
  map_input_spectrum_ch_str[4] = "numuCC_PC_norm";
  map_input_spectrum_ch_str[5] = "CCpi0_FC_norm";
  map_input_spectrum_ch_str[6] = "CCpi0_PC_norm";
  map_input_spectrum_ch_str[7] = "NCpi0_norm";
  map_input_spectrum_ch_str[8] = "Lee_FC";
  map_input_spectrum_ch_str[9] = "Lee_PC"; 
  map_input_spectrum_ch_str[10]= "nueCC_FC_ext";
  map_input_spectrum_ch_str[11]= "nueCC_PC_ext";
  map_input_spectrum_ch_str[12]= "numuCC_FC_ext";
  map_input_spectrum_ch_str[13]= "numuCC_PC_ext";
  map_input_spectrum_ch_str[14]= "CCpi0_FC_ext";
  map_input_spectrum_ch_str[15]= "CCpi0_PC_ext";
  map_input_spectrum_ch_str[16]= "NCpi0_ext";

  /// flag for LEE channels corresponding to the cov_input.txt
  map_Lee_ch[8] = 1;
  map_Lee_ch[9] = 1;
  */
  /////////////////////////////////////// case: separate nueCC signal and bkg
  /*
  for(int idx=1; idx<=18; idx++) map_input_spectrum_ch_str[idx] = TString::Format("pred_%02d", idx);
  map_Lee_ch[8] = 1;
  map_Lee_ch[9] = 1;
  */
  /////////////////////////////////////// case: fake data
  /*
  map_input_spectrum_ch_str[1] = "nueCC_FC_norm";
  map_input_spectrum_ch_str[2] = "nueCC_PC_norm";
  map_input_spectrum_ch_str[3] = "numuCC_FC_norm";
  map_input_spectrum_ch_str[4] = "numuCC_PC_norm";
  map_input_spectrum_ch_str[5] = "CCpi0_FC_norm";
  map_input_spectrum_ch_str[6] = "CCpi0_PC_norm";
  map_input_spectrum_ch_str[7] = "NCpi0_norm";
  map_input_spectrum_ch_str[8] = "Lee_FC";
  map_input_spectrum_ch_str[9] = "Lee_PC";
  
  /// flag for LEE channels corresponding to the cov_input.txt
  map_Lee_ch[8] = 1;
  map_Lee_ch[9] = 1;
  */
  /////////////////////////////////////// case: 1u0p and 1uNp
  /* 
  map_input_spectrum_ch_str[1] = "nueCC_FC_norm";
  map_input_spectrum_ch_str[2] = "nueCC_PC_norm";
  map_input_spectrum_ch_str[3] = "numuCC_FC_1u0p_norm";
  map_input_spectrum_ch_str[4] = "numuCC_PC_1u0p_norm";
  map_input_spectrum_ch_str[5] = "numuCC_FC_1uNp_norm";
  map_input_spectrum_ch_str[6] = "numuCC_PC_1uNp_norm";
  map_input_spectrum_ch_str[7] = "CCpi0_FC_norm";
  map_input_spectrum_ch_str[8] = "CCpi0_PC_norm";
  map_input_spectrum_ch_str[9] = "NCpi0_norm";
  map_input_spectrum_ch_str[10] = "Lee_FC";
  map_input_spectrum_ch_str[11] = "Lee_PC";

  /// flag for LEE channels corresponding to the cov_input.txt
  map_Lee_ch[10] = 1;
  map_Lee_ch[11] = 1;
  */ 
  //////////////////
  //////////////////
  
  roostr = spectra_file;
  TFile *file_spectra = new TFile(roostr, "read");

  ///
  TMatrixD *mat_collapse = (TMatrixD*)file_spectra->Get("mat_collapse");
  matrix_transform.Clear();
  matrix_transform.ResizeTo( mat_collapse->GetNrows(), mat_collapse->GetNcols() );
  matrix_transform = (*mat_collapse);

  ///
  cout<<" Predictions"<<endl;
  for(int ich=1; ich<=(int)map_input_spectrum_ch_str.size(); ich++) {
    roostr = TString::Format("histo_%d", ich);
    TH1F *h1_spectrum = (TH1F*)file_spectra->Get(roostr);
    int bins = h1_spectrum->GetNbinsX() + 1;    
    cout<<Form(" %2d  %-20s   bin-num %2d", ich, map_input_spectrum_ch_str[ich].Data(), bins)<<endl;
    
    for(int ibin=1; ibin<=bins; ibin++) {
      double content = h1_spectrum->GetBinContent(ibin);
      map_input_spectrum_ch_bin[ich][ibin-1] = content;
    }
    
  }
  cout<<endl;
  
  ////////////////////
  ////////////////////

  bins_oldworld = 0;
  for(auto it_ch=map_input_spectrum_ch_bin.begin(); it_ch!=map_input_spectrum_ch_bin.end(); it_ch++) {
    int ich = it_ch->first;
      for(int ibin=0; ibin<(int)map_input_spectrum_ch_bin[ich].size(); ibin++) {
	bins_oldworld++;
	int index_oldworld = bins_oldworld - 1;	
	map_input_spectrum_oldworld_bin[ index_oldworld ] = map_input_spectrum_ch_bin[ich][ibin];
	if( map_Lee_ch.find(ich)!=map_Lee_ch.end() ) map_Lee_oldworld[index_oldworld] = 1;
    }// ibin
  }// ich

  ////////////////////////////////////// data
  
  cout<<" Observations"<<endl;
  
  int line_data = -1;
  bins_newworld = 0;
  for(int ich=1; ich<=channels_observation; ich++) {
    roostr = TString::Format("hdata_obsch_%d", ich);
    TH1F *h1_spectrum = (TH1F*)file_spectra->Get(roostr);
    cout<<Form(" %2d  %-20s   bin-num %2d", ich, roostr.Data(), h1_spectrum->GetNbinsX()+1)<<endl;
        
    for(int ibin=1; ibin<=h1_spectrum->GetNbinsX()+1; ibin++) {
      map_data_spectrum_ch_bin[ich][ibin-1] = h1_spectrum->GetBinContent(ibin);

      line_data++;
      bins_newworld++;
      map_data_spectrum_newworld_bin[line_data] = map_data_spectrum_ch_bin[ich][ibin-1];

    }// ibin
  }// ich
  cout<<endl;
  
  ////////////////////////////////////////// flux_Xs
  
  cout<<" Flux and Xs systematics"<<endl;
    
  //https://www.phy.bnl.gov/xqian/talks/wire-cell/LEEana/configurations/cov_input.txt  
  map<int, TFile*>map_file_flux_Xs_frac;  
  map<int, TMatrixD*>map_matrix_flux_Xs_frac;
  
  TMatrixD matrix_flux_Xs_frac(bins_oldworld, bins_oldworld);
  TMatrixD matrix_flux_frac(bins_oldworld, bins_oldworld);
  TMatrixD matrix_Xs_frac(bins_oldworld, bins_oldworld);
  
  for(int idx=syst_cov_flux_Xs_begin; idx<=syst_cov_flux_Xs_end; idx++) {
    roostr = TString::Format(flux_Xs_directory+"cov_%d.root", idx);
    map_file_flux_Xs_frac[idx] = new TFile(roostr, "read");
    map_matrix_flux_Xs_frac[idx] = (TMatrixD*)map_file_flux_Xs_frac[idx]->Get(TString::Format("frac_cov_xf_mat_%d", idx));
    cout<<TString::Format(" %2d %s", idx, roostr.Data())<<endl;

    //if( idx!=17 )
    matrix_flux_Xs_frac += (*map_matrix_flux_Xs_frac[idx]);    
    
    if( idx<=16 ) {// flux
      matrix_flux_frac += (*map_matrix_flux_Xs_frac[idx]);
    }
    else {// interaction
      matrix_Xs_frac += (*map_matrix_flux_Xs_frac[idx]);
    }    
  }
  cout<<endl;  
  
  ////////////////////////////////////////// detector

  cout<<" Detector systematics"<<endl;
  cout << " (FIXME !!!!)" << endl;
    
  map<int, TString>map_detectorfile_str;
  
  map_detectorfile_str[1] = detector_directory+"cov_LYDown.root";
  map_detectorfile_str[2] = detector_directory+"cov_LYRayleigh.root";
  map_detectorfile_str[3] = detector_directory+"cov_Recomb2.root";
  map_detectorfile_str[4] = detector_directory+"cov_SCE.root";
  //map_detectorfile_str[5] = detector_directory+"cov_WMdEdx.root";
  map_detectorfile_str[6] = detector_directory+"cov_WMThetaXZ.root";
  map_detectorfile_str[7] = detector_directory+"cov_WMThetaYZ.root";
  map_detectorfile_str[8] = detector_directory+"cov_WMX.root";
  map_detectorfile_str[9] = detector_directory+"cov_WMYZ.root";
  map_detectorfile_str[10]= detector_directory+"cov_LYatt.root";
  
  
  map<int, TFile*>map_file_detector_frac;
  map<int, TMatrixD*>map_matrix_detector_frac;
  TMatrixD matrix_detector_frac(bins_oldworld, bins_oldworld);
  map<int, TMatrixD>matrix_detector_sub_frac;
  
  for( auto it=map_detectorfile_str.begin(); it!=map_detectorfile_str.end(); it++ ) {
    int idx = it->first;
    if(idx==5)  continue;    
    roostr = map_detectorfile_str[idx];
    cout<<TString::Format(" %2d %s", idx, roostr.Data())<<endl;
    
    map_file_detector_frac[idx] = new TFile(roostr, "read");
    map_matrix_detector_frac[idx] = (TMatrixD*)map_file_detector_frac[idx]->Get(TString::Format("frac_cov_det_mat_%d", idx));

    matrix_detector_frac += (*map_matrix_detector_frac[idx]);

    matrix_detector_sub_frac[idx].Clear();
    matrix_detector_sub_frac[idx].ResizeTo(bins_oldworld, bins_oldworld);
    matrix_detector_sub_frac[idx] = (*map_matrix_detector_frac[idx]);
  }
  cout<<endl;
        
  ////////////////////////////////////////// additional

  TMatrixD *matrix_additional_abs_point = (TMatrixD*)file_spectra->Get("cov_mat_add");
  TMatrixD matrix_additional_abs = (*matrix_additional_abs_point);
    
  //////////////////////////////////////////

  matrix_input_cov_flux_Xs.Clear();
  matrix_input_cov_flux.Clear();
  matrix_input_cov_Xs.Clear();
  matrix_input_cov_detector.Clear();
  matrix_input_cov_additional.Clear();
  
  matrix_input_cov_flux_Xs.ResizeTo( bins_oldworld, bins_oldworld );
  matrix_input_cov_flux.ResizeTo( bins_oldworld, bins_oldworld );
  matrix_input_cov_Xs.ResizeTo( bins_oldworld, bins_oldworld );  
  matrix_input_cov_detector.ResizeTo( bins_oldworld, bins_oldworld );
  matrix_input_cov_additional.ResizeTo( bins_oldworld, bins_oldworld );

  for(auto it=matrix_detector_sub_frac.begin(); it!=matrix_detector_sub_frac.end(); it++) {
    int idx = it->first;
    matrix_input_cov_detector_sub[idx].Clear();
    matrix_input_cov_detector_sub[idx].ResizeTo( bins_oldworld, bins_oldworld );
  }
  
  for(int ibin=0; ibin<bins_oldworld; ibin++) {
    for(int jbin=0; jbin<bins_oldworld; jbin++) {
      double val_i = map_input_spectrum_oldworld_bin[ibin];
      double val_j = map_input_spectrum_oldworld_bin[jbin];
      double val_cov = 0;
      
      val_cov = matrix_flux_Xs_frac(ibin, jbin);
      matrix_input_cov_flux_Xs(ibin, jbin) = val_cov * val_i * val_j;

      val_cov = matrix_flux_frac(ibin, jbin);
      matrix_input_cov_flux(ibin, jbin) = val_cov * val_i * val_j;

      val_cov = matrix_Xs_frac(ibin, jbin);
      matrix_input_cov_Xs(ibin, jbin) = val_cov * val_i * val_j;      
      
      val_cov = matrix_detector_frac(ibin, jbin);
      matrix_input_cov_detector(ibin, jbin) = val_cov * val_i * val_j;
      
      for(auto it=matrix_input_cov_detector_sub.begin(); it!=matrix_input_cov_detector_sub.end(); it++) {
	int idx = it->first;
	val_cov = matrix_detector_sub_frac[idx](ibin, jbin);
	matrix_input_cov_detector_sub[idx](ibin, jbin) = val_cov * val_i * val_j;
      }
  
    }
  }

  matrix_input_cov_additional = matrix_additional_abs;
  
  ////////////////////////////////////////// MC statistics

  int mc_file_begin = syst_cov_mc_stat_begin;
  int mc_file_end = syst_cov_mc_stat_end;
  
  cout<<TString::Format(" MC statistics. Files:  %d.log - %d.log", mc_file_begin, mc_file_end)<<endl;
  
  map<int, map<int, double> >map_mc_stat_file_bin_Lee;
  map<int, map<int, double> >map_mc_stat_file_bin_mcStat;
  int gbins_mc_stat = 0;
    
  for(int ifile=mc_file_begin; ifile<=mc_file_end; ifile++) {
    roostr = TString::Format(mc_directory+"%d.log", ifile);
    
    ifstream InputFile_aa(roostr, ios::in);
    if(!InputFile_aa) { cerr<<" No input-list"<<endl; exit(1); }

    /////////////////////// check

    int count_check = 0;
    string line_check;    
    ifstream file_check(roostr);
    while( getline(file_check, line_check) ) count_check++;
    //cout<<endl<<" Numbers of lines in the mc_stat file: "<<count_check<<endl<<endl;
    
    gbins_mc_stat = count_check -1;
    if( gbins_mc_stat!=bins_newworld ) {
      cout<<" Error gbins_mc_stat!=bins_newworld: "<<roostr<<endl;
      cerr<<" Error gbins_mc_stat!=bins_newworld: "<<roostr<<endl;
      exit(1);
    }
    
    ///////////////////////
    
    int line = 0;    
    double Lee = 1; double run = 1;
    
    for(int idx=1; idx<=gbins_mc_stat+1; idx++) {            
      int gbin = 0; int lbin = 0; double val_pred = 0; double mc_stat = 0; double nn_stat = 0;
      if(idx==1) { InputFile_aa>>Lee>>run; }
      else {
	InputFile_aa>>gbin>>lbin>>val_pred>>mc_stat>>nn_stat;
	line++;
	map_mc_stat_file_bin_Lee[ifile][line-1] = Lee;
	map_mc_stat_file_bin_mcStat[ifile][line-1] = mc_stat;
      }
    }
  }
  
  /// gh_mc_stat_bin
  for(int ibin=0; ibin<gbins_mc_stat; ibin++) {
    gh_mc_stat_bin[ibin] = new TGraph(); gh_mc_stat_bin[ibin]->SetName(TString::Format("gh_mc_stat_bin_%03d", ibin));
    
    for(auto it=map_mc_stat_file_bin_Lee.begin(); it!=map_mc_stat_file_bin_Lee.end(); it++) {
      int ifile = it->first;
      double Lee = map_mc_stat_file_bin_Lee[ifile][ibin];
      double mc_stat = map_mc_stat_file_bin_mcStat[ifile][ibin];
      gh_mc_stat_bin[ibin]->SetPoint( gh_mc_stat_bin[ibin]->GetN(), Lee, mc_stat );
    }
    
    double x,y;
    gh_mc_stat_bin[ibin]->GetPoint( gh_mc_stat_bin[ibin]->GetN()-1, x, y);
    gh_mc_stat_bin[ibin]->SetPoint( gh_mc_stat_bin[ibin]->GetN(), x+1, y);
  }  

  ///////////////////

  matrix_data_MCstat_correlation.Clear();
  matrix_data_MCstat_correlation.ResizeTo(bins_newworld, bins_newworld);
  
  matrix_pred_MCstat_correlation.Clear();
  matrix_pred_MCstat_correlation.ResizeTo(bins_newworld, bins_newworld);
  
  TFile *roofile_data = new TFile(mc_directory+"../"+"run_pred_stat.root", "read");
  TMatrixD *matrix_datastat = (TMatrixD*)roofile_data->Get("cov_mat_0");
  // TMatrixD *matrix_datastat = new TMatrixD(bins_newworld, bins_newworld); // empty matrix

  cout<<endl<<" check matrix_datastat "<<matrix_datastat->GetNrows()<<endl<<endl;
  
  for(int idx=0; idx<matrix_datastat->GetNrows(); idx++ ) {
    for(int jdx=0; jdx<matrix_datastat->GetNrows(); jdx++ ) {
      double val_cov = (*matrix_datastat)(idx, jdx);      
      double ii = sqrt((*matrix_datastat)(idx,idx));
      double jj = sqrt((*matrix_datastat)(jdx,jdx));
      if( ii!=0 && jj!=0 ) matrix_data_MCstat_correlation(idx, jdx) = val_cov/ii/jj;
      if( idx==jdx) matrix_data_MCstat_correlation(idx, jdx) = 1;
    }
  }

  matrix_pred_MCstat_correlation = matrix_data_MCstat_correlation;

  //matrix_pred_MCstat_correlation.Draw("colz");
    
  /////////////////////////////////////////////////////////////////
    
  cout<<endl;
  cout<<" ---> Complete the initialization"<<endl;
  cout<<" ---> Complete the initialization"<<endl;  
  cout<<endl;
  
}
