  #include <admodel.h>
  #include <time.h>
  void solveBaranov( const int& t, const int& iter);
  // Flag to control whether header written to mcmc output file.
  int mcmcHeader = 0;
  
  // Flag to control which parameters are done in MCMC phase.
  int mcmcFlag = 1;
  //ofstream mcout("mcout.dat");
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <sableOpMod.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  likeYear.allocate("likeYear");
 ad_comm::change_datafile_name("opModCatch.dat");
  nT.allocate("nT");
  nFisheries.allocate("nFisheries");
  firstYearCatch.allocate(1,nFisheries-2,"firstYearCatch");
  landCatchMatrix.allocate(1,nT,1,nFisheries+2,"landCatchMatrix");
  releaseMatrix.allocate(1,nT,1,nFisheries+2,"releaseMatrix");
  totLandCatch.allocate(1,nT);
  propCatchFish.allocate(1,nFisheries);
 ad_comm::change_datafile_name("opModIndex.dat");
  nIndexSeries.allocate("nIndexSeries");
  idxIndex.allocate(1,nIndexSeries,"idxIndex");
  idxLikeWeight.allocate(1,nIndexSeries,"idxLikeWeight");
  idxFirstYear.allocate(1,nIndexSeries,"idxFirstYear");
  idxLastYear.allocate(1,nIndexSeries,"idxLastYear");
  fracYearSurvey.allocate(1,nIndexSeries,"fracYearSurvey");
  idxSeries.allocate(1,nIndexSeries,1,nT,"idxSeries");
  nobs.allocate(1,nIndexSeries);
 nobs=idxLastYear-idxFirstYear+1;
  validObs.allocate(1,nIndexSeries);
 ad_comm::change_datafile_name("opModTVq.dat");
  firstqdev_g.allocate(1,nIndexSeries,"firstqdev_g");
  lastqdev_g.allocate(1,nIndexSeries,"lastqdev_g");
  ph_log_qdev_g.allocate(1,nIndexSeries,"ph_log_qdev_g");
  priorSD_qdev.allocate(1,nIndexSeries,"priorSD_qdev");
  checkTVq.allocate("checkTVq");
 if(checkTVq != 999) cout << "Bad Check in TVq!" << endl;
 if(checkTVq != 999) exit(1);
 ad_comm::change_datafile_name("opModAges.dat");
  nAgeSeries.allocate("nAgeSeries");
 cout << "nAgeSeries == " << nAgeSeries << endl;
  plusGroupAge.allocate("plusGroupAge");
  ages.allocate(1,plusGroupAge);
 ages.fill_seqadd(1,1);
  minAge.allocate(1,nAgeSeries,"minAge");
  maxAge.allocate(1,nAgeSeries,"maxAge");
  fisheryAgeID.allocate(1,nAgeSeries,"fisheryAgeID");
  ageLikeWeight_m.allocate(1,nAgeSeries,"ageLikeWeight_m");
  ageLikeWeight_f.allocate(1,nAgeSeries,"ageLikeWeight_f");
  fracYearAges.allocate(1,nAgeSeries,"fracYearAges");
  useAgeError_m.allocate(1,nAgeSeries,"useAgeError_m");
  useAgeError_f.allocate(1,nAgeSeries,"useAgeError_f");
  firstAge_m.allocate(1,nAgeSeries,1,nT,"firstAge_m");
  firstAge_f.allocate(1,nAgeSeries,1,nT,"firstAge_f");
  lastAge_m.allocate(1,nAgeSeries,1,nT,"lastAge_m");
  lastAge_f.allocate(1,nAgeSeries,1,nT,"lastAge_f");
  nObsAge_m.allocate(1,nAgeSeries,1,nT,"nObsAge_m");
  nObsAge_f.allocate(1,nAgeSeries,1,nT,"nObsAge_f");
  ageObsProp_m.allocate(1,nAgeSeries,1,nT,minAge,maxAge,"ageObsProp_m");
  ageObsProp_f.allocate(1,nAgeSeries,1,nT,minAge,maxAge,"ageObsProp_f");
  prop_threshold = 0.005;
  accBin_m.allocate(1,nAgeSeries,1,nT);
  accBin_f.allocate(1,nAgeSeries,1,nT);
  eof_Ages.allocate("eof_Ages");
 cout << "opModAges/eof_Ages = " << eof_Ages << endl;
 if(eof_Ages != 123) cout << "Bad Check in Ages!" << endl;
 if(eof_Ages != 123) exit(1);
 ad_comm::change_datafile_name("ageErrorMatrix.dat");
  Q.allocate(2,35,2,plusGroupAge,"Q");
 ad_comm::change_datafile_name("opModSelPriors.dat");
  firstSelPrior_g.allocate(1,nFisheries-2);
  lastSelPrior_g.allocate(1,nFisheries-2);
  firstSelPrior_1.allocate("firstSelPrior_1");
 firstSelPrior_g(1) = firstSelPrior_1;
  firstSelPrior_2.allocate("firstSelPrior_2");
 firstSelPrior_g(2) = firstSelPrior_2;
  firstSelPrior_3.allocate("firstSelPrior_3");
 firstSelPrior_g(3) = firstSelPrior_3;
  lastSelPrior_1.allocate("lastSelPrior_1");
 lastSelPrior_g(1) = lastSelPrior_1;
  lastSelPrior_2.allocate("lastSelPrior_2");
 lastSelPrior_g(2) = lastSelPrior_2;
  lastSelPrior_3.allocate("lastSelPrior_3");
 lastSelPrior_g(3) = lastSelPrior_3;
  ph_log_alpha_dev_g.allocate(1,nFisheries-2,"ph_log_alpha_dev_g");
  ph_log_beta_dev_g.allocate(1,nFisheries-2,"ph_log_beta_dev_g");
  priorSD_alpha_dev.allocate("priorSD_alpha_dev");
  priorSD_beta_dev.allocate("priorSD_beta_dev");
  ph_alpha_sexdev_g.allocate(1,nFisheries,"ph_alpha_sexdev_g");
  ph_beta_sexdev_g.allocate(1,nFisheries,"ph_beta_sexdev_g");
  priorSD_alpha_sexdev.allocate("priorSD_alpha_sexdev");
  priorSD_beta_sexdev.allocate("priorSD_beta_sexdev");
  p_alpha_1t.allocate(firstSelPrior_1,lastSelPrior_1,"p_alpha_1t");
  p_alpha_2t.allocate(firstSelPrior_2,lastSelPrior_2,"p_alpha_2t");
  p_alpha_3t.allocate(firstSelPrior_3,lastSelPrior_3,"p_alpha_3t");
  p_beta_1t.allocate(firstSelPrior_1,lastSelPrior_1,"p_beta_1t");
  p_beta_2t.allocate(firstSelPrior_2,lastSelPrior_2,"p_beta_2t");
  p_beta_3t.allocate(firstSelPrior_3,lastSelPrior_3,"p_beta_3t");
  sd_alpha_1t.allocate(firstSelPrior_1,lastSelPrior_1,"sd_alpha_1t");
  sd_alpha_2t.allocate(firstSelPrior_2,lastSelPrior_2,"sd_alpha_2t");
  sd_alpha_3t.allocate(firstSelPrior_3,lastSelPrior_3,"sd_alpha_3t");
  sd_beta_1t.allocate(firstSelPrior_1,lastSelPrior_1,"sd_beta_1t");
  sd_beta_2t.allocate(firstSelPrior_2,lastSelPrior_2,"sd_beta_2t");
  sd_beta_3t.allocate(firstSelPrior_3,lastSelPrior_3,"sd_beta_3t");
  checkSelPriors.allocate("checkSelPriors");
 if(checkSelPriors != 999) cout << "Bad Check in selPriors!" << endl;
 if(checkSelPriors != 999) exit(1);
 ad_comm::change_datafile_name("opModSelBlocks.dat");
  nSelBlocks_g1.allocate("nSelBlocks_g1");
  nSelBlocks_g2.allocate("nSelBlocks_g2");
  nSelBlocks_g3.allocate("nSelBlocks_g3");
  nSelBlocks.allocate(1,3);
 nSelBlocks(1) = nSelBlocks_g1;
 nSelBlocks(2) = nSelBlocks_g2;
 nSelBlocks(3) = nSelBlocks_g3;
  ph_alpha_blockDev.allocate(1,3,"ph_alpha_blockDev");
  ph_beta_blockDev.allocate(1,3,"ph_beta_blockDev");
  blockDevPriorSd.allocate(1,3,"blockDevPriorSd");
  selBlockStart_g1.allocate(1,nSelBlocks_g1,"selBlockStart_g1");
  selBlockEnd_g1.allocate(1,nSelBlocks_g1,"selBlockEnd_g1");
  selBlockStart_g2.allocate(1,nSelBlocks_g2,"selBlockStart_g2");
  selBlockEnd_g2.allocate(1,nSelBlocks_g2,"selBlockEnd_g2");
  selBlockStart_g3.allocate(1,nSelBlocks_g3,"selBlockStart_g3");
  selBlockEnd_g3.allocate(1,nSelBlocks_g3,"selBlockEnd_g3");
  checkSelBlock.allocate("checkSelBlock");
 if(checkSelBlock != 999) cout << "Bad Check in selBlocks!" << endl;
 if(checkSelBlock != 999) exit(1);
 ad_comm::change_datafile_name("opModPropDirected.dat");
  pDirOn.allocate(1,nFisheries-2,"pDirOn");
  pDir_t1_in.allocate(1,nFisheries-2,"pDir_t1_in");
  pDir_t1.allocate(1,nFisheries-2);
 pDir_t1 = pDir_t1_in;
  pDir_t2_in.allocate(1,nFisheries-2,"pDir_t2_in");
  pDir_t2.allocate(1,nFisheries-2);
 pDir_t2= pDir_t2_in;
  pDirScale.allocate("pDirScale");
  selAlphaDirDev_base.allocate(1,nFisheries-2,"selAlphaDirDev_base");
  ph_logit_pDir_in.allocate(1,nFisheries-2,"ph_logit_pDir_in");
  ph_logit_pDir.allocate(1,nFisheries-2);
 ph_logit_pDir = ph_logit_pDir_in;
  ph_log_selAlphaDir_in.allocate(1,nFisheries-2,"ph_log_selAlphaDir_in");
  ph_log_selAlphaDirDev.allocate(1,nFisheries-2);
 ph_log_selAlphaDirDev = ph_log_selAlphaDir_in;
  ph_log_selBetaDir_in.allocate(1,nFisheries-2,"ph_log_selBetaDir_in");
  ph_log_selBetaDirDev.allocate(1,nFisheries-2);
 ph_log_selBetaDirDev = ph_log_selBetaDir_in;
  priorSD_logit_pDir.allocate("priorSD_logit_pDir");
  priorSD_seldev_pDir.allocate("priorSD_seldev_pDir");
  penaliseUndir.allocate("penaliseUndir");
  undirK.allocate("undirK");
  checkPDir.allocate("checkPDir");
 if(checkPDir != 999) cout << "Bad Check in pDir!" << endl;
 if(checkPDir != 999) exit(1);
 ad_comm::change_datafile_name("opModReportRate.dat");
  lenRep95.allocate(1,nFisheries-2,"lenRep95");
  lenRep50.allocate(1,nFisheries-2,"lenRep50");
  phzLenRep50_g.allocate(1,nFisheries-2,"phzLenRep50_g");
  phzLenRep95_g.allocate(1,nFisheries-2,"phzLenRep95_g");
  checkRepRate.allocate("checkRepRate");
 if(checkPDir != 999) cout << "Bad Check in report rates!" << endl;
 if(checkPDir != 999) exit(1);
 ad_comm::change_datafile_name("opModControl.ctl");
  verbose.allocate("verbose");
  useBaranov.allocate("useBaranov");
 ph_logF = 1;
 if( useBaranov ) ph_logF = -1;
  baranovIter.allocate("baranovIter");
  baranovSteps.allocate(1,nFisheries,"baranovSteps");
  objFunc_scale.allocate("objFunc_scale");
  objFunc_scale_last.allocate("objFunc_scale_last");
  objFunc_scale_sd.allocate("objFunc_scale_sd");
  sizeLim.allocate("sizeLim");
  sizeLimYear.allocate("sizeLimYear");
  dM.allocate(1,nFisheries,"dM");
  lInf_m.allocate("lInf_m");
  lInf_f.allocate("lInf_f");
  vonK_m.allocate("vonK_m");
  vonK_f.allocate("vonK_f");
  sigmaL_m.allocate("sigmaL_m");
  sigmaL_f.allocate("sigmaL_f");
  L1_m.allocate("L1_m");
  L1_f.allocate("L1_f");
  wt_a.allocate(1,2,"wt_a");
  wt_b.allocate(1,2,"wt_b");
  aMat50.allocate("aMat50");
  aMat95.allocate("aMat95");
  ph_log_avgR.allocate("ph_log_avgR");
  ph_logit_h.allocate("ph_logit_h");
  lb_logit_h.allocate("lb_logit_h");
  ub_logit_h.allocate("ub_logit_h");
  prior_h.allocate(1,2,"prior_h");
  ph_log_SSB0.allocate("ph_log_SSB0");
  lb_log_SSB0.allocate("lb_log_SSB0");
  ub_log_SSB0.allocate("ub_log_SSB0");
  ph_log_M.allocate("ph_log_M");
  ph_log_q.allocate("ph_log_q");
  phRecDevs.allocate("phRecDevs");
  firstRecDev.allocate("firstRecDev");
  lastRecDev.allocate("lastRecDev");
  selType.allocate(1,nFisheries,"selType");
 cout << "selType = " << selType << endl;
  ph_log_alpha_g.allocate(1,nFisheries,"ph_log_alpha_g");
  lb_log_alpha_g.allocate(1,nFisheries,"lb_log_alpha_g");
  ub_log_alpha_g.allocate(1,nFisheries,"ub_log_alpha_g");
  ph_log_beta_g.allocate(1,nFisheries,"ph_log_beta_g");
  lb_log_beta_g.allocate(1,nFisheries,"lb_log_beta_g");
  ub_log_beta_g.allocate(1,nFisheries,"ub_log_beta_g");
  useHighgrading.allocate(1,nFisheries-2,"useHighgrading");
  ph_log_hg50.allocate("ph_log_hg50");
  ageLikeType.allocate("ageLikeType");
  sigma_R.allocate("sigma_R");
  sigma_alpha.allocate(1,nFisheries,"sigma_alpha");
  sigma_beta.allocate(1,nFisheries,"sigma_beta");
  tauRelLambda.allocate(1,nFisheries-2,"tauRelLambda");
  tauIndexLambda.allocate(1,nIndexSeries,"tauIndexLambda");
  tauAgeLambda.allocate(1,nFisheries-1,"tauAgeLambda");
 cout << "tauAgeLambda = " << tauAgeLambda << endl;
  prior_mean_M.allocate(1,2,"prior_mean_M");
  prior_sd_M.allocate(1,2,"prior_sd_M");
  priorSD_F.allocate("priorSD_F");
  catError.allocate("catError");
  nFixedRecDev.allocate("nFixedRecDev");
  tFixedRecDev.allocate(1,nFixedRecDev,"tFixedRecDev");
  logFixedRecDev.allocate(1,nFixedRecDev,"logFixedRecDev");
  fracYearFleet.allocate(1,nFisheries,"fracYearFleet");
 cout << "fracYearFleet = " << fracYearFleet << endl;
  checkCtl.allocate("checkCtl");
 cout << "checkCtl = " << checkCtl << endl;
 if(checkCtl != 999) cout << "Bad Check in opModControl!" << endl;
 if(checkCtl != 999) exit(1);
 ad_comm::change_datafile_name("opModRelDevs.dat");
  relDevInitYr.allocate(1,nFisheries-2,"relDevInitYr");
  relDevLastYr.allocate(1,nFisheries-2,"relDevLastYr");
  relDevPenalty.allocate("relDevPenalty");
  checkRelDevs.allocate("checkRelDevs");
 if(checkRelDevs != 999) cout << "Bad Check in relDevs!" << endl;
 if(checkRelDevs != 999) exit(1);
  matAge_f.allocate(1,plusGroupAge);
  lenAge_m.allocate(1,plusGroupAge);
  lenAge_f.allocate(1,plusGroupAge);
  wtAge_m.allocate(1,plusGroupAge);
  wtAge_f.allocate(1,plusGroupAge);
  Ctg.allocate(1,nT,1,nFisheries);
  obs_relCtg.allocate(1,nT,1,nFisheries);
  sublegalDtg.allocate(1,nFisheries,1,nT);
  pRet_m.allocate(1,nT,1,plusGroupAge);
  pRet_f.allocate(1,nT,1,plusGroupAge);
  validObs_idx.allocate(1,nIndexSeries);
  validObs_rel.allocate(1,nFisheries-2);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  objFunc.allocate("objFunc");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  log_SSB0.allocate(lb_log_SSB0,ub_log_SSB0,ph_log_SSB0,"log_SSB0");
 cout << "log_SSB0 = " << log_SSB0 << endl;
  logit_h.allocate(lb_logit_h,ub_logit_h,ph_logit_h,"logit_h");
 cout << "logit_h = " << logit_h << endl;
  log_avgR.allocate(0.7,2.,ph_log_avgR,"log_avgR");
 cout << "log_avgR = " << log_avgR << endl;
  log_M.allocate(1,2,-5,-1,ph_log_M,"log_M");
 cout << "log_M = " << log_M << endl;
  log_Fg1.allocate(1,nFisheries-2,-10.0,-0.7,ph_logF,"log_Fg1");
 cout << "log_Fg1 = " << log_Fg1 << endl;
  log_Fdevs_1t.allocate(10,nT,-5.,5.,ph_logF,"log_Fdevs_1t");
  log_Fdevs_2t.allocate(2,nT,-5.,5.,ph_logF,"log_Fdevs_2t");
  log_Fdevs_3t.allocate(2,nT,-5.,5.,ph_logF,"log_Fdevs_3t");
  logit_rho.allocate(1,nFisheries-2,-10.,10.,-1,"logit_rho");
  recDevs.allocate(firstRecDev,lastRecDev,-3.,3.,phRecDevs,"recDevs");
 cout << "recDevs = " << recDevs << endl;
  log_alpha_g.allocate(1,nFisheries,lb_log_alpha_g,ub_log_alpha_g,ph_log_alpha_g,"log_alpha_g");
 cout << "log_alpha_g = " << log_alpha_g << endl;
  log_beta_g.allocate(1,nFisheries,lb_log_beta_g,ub_log_beta_g,ph_log_beta_g,"log_beta_g");
 cout << "log_beta_g = " << log_beta_g << endl;
  log_alpha_dev_gt.allocate(1,nFisheries-2,firstSelPrior_g,lastSelPrior_g,-2,2,ph_log_alpha_dev_g,"log_alpha_dev_gt");
  log_beta_dev_gt.allocate(1,nFisheries-2,firstSelPrior_g,lastSelPrior_g,-1.,1.,ph_log_beta_dev_g,"log_beta_dev_gt");
 cout << "ph_log_alpha_dev_g = " << ph_log_alpha_dev_g << endl;
 cout << "ph_log_beta_dev_g = " << ph_log_beta_dev_g << endl;
  lnq.allocate(1,nIndexSeries,-2,2,ph_log_q,"lnq");
 cout << "lnq = " << lnq << endl;
  log_qdev_gt.allocate(1,nIndexSeries,firstqdev_g,lastqdev_g,-1.,1.,ph_log_qdev_g,"log_qdev_gt");
 cout << "log_qdev_1 = " << log_qdev_gt(1) << endl;
 cout << "log_qdev_2 = " << log_qdev_gt(2) << endl;
 cout << "log_qdev_3 = " << log_qdev_gt(3) << endl;
  log_hg50.allocate(1,nFisheries-2,4.0,4.1,ph_log_hg50,"log_hg50");
 cout << "log_hg50 = " << log_hg50 << endl;
  log_alpha_blockDev_g.allocate(1,nFisheries-2,1,nSelBlocks,-2,2,ph_alpha_blockDev,"log_alpha_blockDev_g");
 cout << "log_alpha_blockDev_1 = " << log_alpha_blockDev_g(1) << endl;
 cout << "log_alpha_blockDev_2 = " << log_alpha_blockDev_g(2) << endl;
 cout << "log_alpha_blockDev_3 = " << log_alpha_blockDev_g(3) << endl;
  log_beta_blockDev_g.allocate(1,nFisheries-2,1,nSelBlocks,-2,2,ph_beta_blockDev,"log_beta_blockDev_g");
  log_alpha_sexdev_g.allocate(1,nFisheries,-2,2,ph_alpha_sexdev_g,"log_alpha_sexdev_g");
  log_beta_sexdev_g.allocate(1,nFisheries,-2,2,ph_beta_sexdev_g,"log_beta_sexdev_g");
 cout << "log_alpha_sexdev_1 = " << log_alpha_sexdev_g(1) << endl;
 cout << "log_alpha_sexdev_2 = " << log_alpha_sexdev_g(2) << endl;
 cout << "log_alpha_sexdev_3 = " << log_alpha_sexdev_g(3) << endl;
  logit_pDir_g.allocate(1,nFisheries-2,-10.,10.,ph_logit_pDir,"logit_pDir_g");
 cout << "logit_pDir_1 = " << logit_pDir_g(1) << endl;
 cout << "logit_pDir_2 = " << logit_pDir_g(2) << endl;
 cout << "logit_pDir_3 = " << logit_pDir_g(3) << endl;
  log_selAlphaDirDev_g.allocate(1,nFisheries-2,-5,5,ph_log_selAlphaDirDev,"log_selAlphaDirDev_g");
 cout << "log_selAlphaDirDev_1 = " << log_selAlphaDirDev_g(1) << endl;
 cout << "log_selAlphaDirDev_2 = " << log_selAlphaDirDev_g(2) << endl;
 cout << "log_selAlphaDirDev_3 = " << log_selAlphaDirDev_g(3) << endl;
  log_selBetaDirDev_g.allocate(1,nFisheries-2,-5,5,ph_log_selBetaDirDev,"log_selBetaDirDev_g");
 cout << "logit_pDir_g = " << logit_pDir_g << endl;
 cout << "log_selAlphaDirDev_g = " << log_selAlphaDirDev_g << endl;
 cout << "log_selBetaDirDev_g = " << log_selBetaDirDev_g << endl;
  log_epsLenRep50_g.allocate(1,nFisheries-2,-10.,10.,phzLenRep50_g,"log_epsLenRep50_g");
  log_epsLenRep95_g.allocate(1,nFisheries-2,-10.,10.,phzLenRep95_g,"log_epsLenRep95_g");
  reportRate_ga_m.allocate(1,nFisheries,1,plusGroupAge,"reportRate_ga_m");
  #ifndef NO_AD_INITIALIZE
    reportRate_ga_m.initialize();
  #endif
  reportRate_ga_f.allocate(1,nFisheries,1,plusGroupAge,"reportRate_ga_f");
  #ifndef NO_AD_INITIALIZE
    reportRate_ga_f.initialize();
  #endif
  B0.allocate("B0");
  #ifndef NO_AD_INITIALIZE
  B0.initialize();
  #endif
  SSB0.allocate("SSB0");
  R0.allocate("R0");
  D2015.allocate("D2015");
  B2015.allocate("B2015");
  h.allocate("h");
  M_m.allocate("M_m");
  M_f.allocate("M_f");
  LHR2015.allocate("LHR2015");
  SLHR2015.allocate("SLHR2015");
  legalBt.allocate(1,nT+1,"legalBt");
  subLegalBt.allocate(1,nT+1,"subLegalBt");
  landedCtg.allocate(1,nT,1,nFisheries,"landedCtg");
  #ifndef NO_AD_INITIALIZE
    landedCtg.initialize();
  #endif
  undirLandCtg.allocate(1,nT,1,nFisheries,"undirLandCtg");
  #ifndef NO_AD_INITIALIZE
    undirLandCtg.initialize();
  #endif
  dirLandCtg.allocate(1,nT,1,nFisheries,"dirLandCtg");
  #ifndef NO_AD_INITIALIZE
    dirLandCtg.initialize();
  #endif
  undirRelCtg.allocate(1,nT,1,nFisheries,"undirRelCtg");
  #ifndef NO_AD_INITIALIZE
    undirRelCtg.initialize();
  #endif
  dirRelCtg.allocate(1,nT,1,nFisheries,"dirRelCtg");
  #ifndef NO_AD_INITIALIZE
    dirRelCtg.initialize();
  #endif
  subLegalHR.allocate(1,nT,"subLegalHR");
  #ifndef NO_AD_INITIALIZE
    subLegalHR.initialize();
  #endif
  legalHR.allocate(1,nT,"legalHR");
  #ifndef NO_AD_INITIALIZE
    legalHR.initialize();
  #endif
  log_Fdevs_gt.allocate(1,nFisheries-2,1,nT,"log_Fdevs_gt");
  #ifndef NO_AD_INITIALIZE
    log_Fdevs_gt.initialize();
  #endif
  rho.allocate(1,nFisheries-2,"rho");
  #ifndef NO_AD_INITIALIZE
    rho.initialize();
  #endif
  avgR.allocate("avgR");
  #ifndef NO_AD_INITIALIZE
  avgR.initialize();
  #endif
  phiSSB.allocate("phiSSB");
  #ifndef NO_AD_INITIALIZE
  phiSSB.initialize();
  #endif
  bh_a.allocate("bh_a");
  #ifndef NO_AD_INITIALIZE
  bh_a.initialize();
  #endif
  bh_b.allocate("bh_b");
  #ifndef NO_AD_INITIALIZE
  bh_b.initialize();
  #endif
  ch_sigmaF.allocate(1,3,1,3,"ch_sigmaF");
  #ifndef NO_AD_INITIALIZE
    ch_sigmaF.initialize();
  #endif
  sigmaF.allocate(1,3,1,3,"sigmaF");
  #ifndef NO_AD_INITIALIZE
    sigmaF.initialize();
  #endif
  row_means_F.allocate(1,3,"row_means_F");
  #ifndef NO_AD_INITIALIZE
    row_means_F.initialize();
  #endif
  alpha_gt_m.allocate(1,nFisheries,1,nT,"alpha_gt_m");
  #ifndef NO_AD_INITIALIZE
    alpha_gt_m.initialize();
  #endif
  beta_gt_m.allocate(1,nFisheries,1,nT,"beta_gt_m");
  #ifndef NO_AD_INITIALIZE
    beta_gt_m.initialize();
  #endif
  alpha_gt_f.allocate(1,nFisheries,1,nT,"alpha_gt_f");
  #ifndef NO_AD_INITIALIZE
    alpha_gt_f.initialize();
  #endif
  beta_gt_f.allocate(1,nFisheries,1,nT,"beta_gt_f");
  #ifndef NO_AD_INITIALIZE
    beta_gt_f.initialize();
  #endif
  alphaUndir_gt_m.allocate(1,nFisheries,1,nT,"alphaUndir_gt_m");
  #ifndef NO_AD_INITIALIZE
    alphaUndir_gt_m.initialize();
  #endif
  alphaUndir_gt_f.allocate(1,nFisheries,1,nT,"alphaUndir_gt_f");
  #ifndef NO_AD_INITIALIZE
    alphaUndir_gt_f.initialize();
  #endif
  betaUndir_gt_m.allocate(1,nFisheries,1,nT,"betaUndir_gt_m");
  #ifndef NO_AD_INITIALIZE
    betaUndir_gt_m.initialize();
  #endif
  betaUndir_gt_f.allocate(1,nFisheries,1,nT,"betaUndir_gt_f");
  #ifndef NO_AD_INITIALIZE
    betaUndir_gt_f.initialize();
  #endif
  alpha_g1.allocate(1,nFisheries,"alpha_g1");
  #ifndef NO_AD_INITIALIZE
    alpha_g1.initialize();
  #endif
  beta_g1.allocate(1,nFisheries,"beta_g1");
  #ifndef NO_AD_INITIALIZE
    beta_g1.initialize();
  #endif
  sel_gta_m.allocate(1,nFisheries,1,nT,1,plusGroupAge,"sel_gta_m");
  #ifndef NO_AD_INITIALIZE
    sel_gta_m.initialize();
  #endif
  sel_gta_f.allocate(1,nFisheries,1,nT,1,plusGroupAge,"sel_gta_f");
  #ifndef NO_AD_INITIALIZE
    sel_gta_f.initialize();
  #endif
  selUndir_gta_m.allocate(1,nFisheries,1,nT,1,plusGroupAge,"selUndir_gta_m");
  #ifndef NO_AD_INITIALIZE
    selUndir_gta_m.initialize();
  #endif
  selUndir_gta_f.allocate(1,nFisheries,1,nT,1,plusGroupAge,"selUndir_gta_f");
  #ifndef NO_AD_INITIALIZE
    selUndir_gta_f.initialize();
  #endif
  lnq_it.allocate(1,nIndexSeries,1,nT,"lnq_it");
  #ifndef NO_AD_INITIALIZE
    lnq_it.initialize();
  #endif
  SSBt.allocate(1,nT+1,"SSBt");
  Rt.allocate(1,nT+1,"Rt");
  #ifndef NO_AD_INITIALIZE
    Rt.initialize();
  #endif
  omega_t.allocate(1,nT,"omega_t");
  #ifndef NO_AD_INITIALIZE
    omega_t.initialize();
  #endif
  Nta_m.allocate(1,nT+1,1,plusGroupAge,"Nta_m");
  #ifndef NO_AD_INITIALIZE
    Nta_m.initialize();
  #endif
  Nta_f.allocate(1,nT+1,1,plusGroupAge,"Nta_f");
  #ifndef NO_AD_INITIALIZE
    Nta_f.initialize();
  #endif
  Bta_m.allocate(1,nT+1,1,plusGroupAge,"Bta_m");
  #ifndef NO_AD_INITIALIZE
    Bta_m.initialize();
  #endif
  Bta_f.allocate(1,nT+1,1,plusGroupAge,"Bta_f");
  #ifndef NO_AD_INITIALIZE
    Bta_f.initialize();
  #endif
  Zta_m.allocate(1,nT+1,1,plusGroupAge,"Zta_m");
  #ifndef NO_AD_INITIALIZE
    Zta_m.initialize();
  #endif
  Zta_f.allocate(1,nT+1,1,plusGroupAge,"Zta_f");
  #ifndef NO_AD_INITIALIZE
    Zta_f.initialize();
  #endif
  pCtg.allocate(1,nT,1,nFisheries-2,"pCtg");
  #ifndef NO_AD_INITIALIZE
    pCtg.initialize();
  #endif
  log_Ftg.allocate(1,nT,1,nFisheries-2,"log_Ftg");
  #ifndef NO_AD_INITIALIZE
    log_Ftg.initialize();
  #endif
  Ftg.allocate(1,nT,1,nFisheries-2,"Ftg");
  #ifndef NO_AD_INITIALIZE
    Ftg.initialize();
  #endif
  uCgta_m.allocate(1,nAgeSeries,1,nT,minAge,plusGroupAge,"uCgta_m");
  #ifndef NO_AD_INITIALIZE
    uCgta_m.initialize();
  #endif
  uCgta_f.allocate(1,nAgeSeries,1,nT,minAge,plusGroupAge,"uCgta_f");
  #ifndef NO_AD_INITIALIZE
    uCgta_f.initialize();
  #endif
  predObsProp_m.allocate(1,nAgeSeries,1,nT,minAge,maxAge,"predObsProp_m");
  #ifndef NO_AD_INITIALIZE
    predObsProp_m.initialize();
  #endif
  predObsProp_f.allocate(1,nAgeSeries,1,nT,minAge,maxAge,"predObsProp_f");
  #ifndef NO_AD_INITIALIZE
    predObsProp_f.initialize();
  #endif
  relCtg.allocate(1,nT,1,nFisheries,"relCtg");
  #ifndef NO_AD_INITIALIZE
    relCtg.initialize();
  #endif
  pR_gta_m.allocate(1,nFisheries,1,nT,1,plusGroupAge,"pR_gta_m");
  #ifndef NO_AD_INITIALIZE
    pR_gta_m.initialize();
  #endif
  pR_gta_f.allocate(1,nFisheries,1,nT,1,plusGroupAge,"pR_gta_f");
  #ifndef NO_AD_INITIALIZE
    pR_gta_f.initialize();
  #endif
  pDir_tg.allocate(1,nT,1,nFisheries-2,"pDir_tg");
  #ifndef NO_AD_INITIALIZE
    pDir_tg.initialize();
  #endif
  expBit.allocate(1,nFisheries,1,nT,"expBit");
  #ifndef NO_AD_INITIALIZE
    expBit.initialize();
  #endif
  indexLikelihood.allocate(1,nIndexSeries,"indexLikelihood");
  #ifndef NO_AD_INITIALIZE
    indexLikelihood.initialize();
  #endif
  ssQ.allocate("ssQ");
  #ifndef NO_AD_INITIALIZE
  ssQ.initialize();
  #endif
  ssR.allocate("ssR");
  #ifndef NO_AD_INITIALIZE
  ssR.initialize();
  #endif
  ss.allocate(1,nIndexSeries,"ss");
  #ifndef NO_AD_INITIALIZE
    ss.initialize();
  #endif
  tauSquareIndex.allocate(1,nIndexSeries,"tauSquareIndex");
  #ifndef NO_AD_INITIALIZE
    tauSquareIndex.initialize();
  #endif
  ageLikelihood_m.allocate(1,nAgeSeries,"ageLikelihood_m");
  #ifndef NO_AD_INITIALIZE
    ageLikelihood_m.initialize();
  #endif
  tauSquareAges_m.allocate(1,nAgeSeries,"tauSquareAges_m");
  #ifndef NO_AD_INITIALIZE
    tauSquareAges_m.initialize();
  #endif
  ageLikelihood_f.allocate(1,nAgeSeries,"ageLikelihood_f");
  #ifndef NO_AD_INITIALIZE
    ageLikelihood_f.initialize();
  #endif
  tauSquareAges_f.allocate(1,nAgeSeries,"tauSquareAges_f");
  #ifndef NO_AD_INITIALIZE
    tauSquareAges_f.initialize();
  #endif
  releaseLikelihood.allocate(1,nFisheries,"releaseLikelihood");
  #ifndef NO_AD_INITIALIZE
    releaseLikelihood.initialize();
  #endif
  tauSquareReleases.allocate(1,nFisheries,"tauSquareReleases");
  #ifndef NO_AD_INITIALIZE
    tauSquareReleases.initialize();
  #endif
  mnLL_m.allocate("mnLL_m");
  #ifndef NO_AD_INITIALIZE
  mnLL_m.initialize();
  #endif
  mnLL_f.allocate("mnLL_f");
  #ifndef NO_AD_INITIALIZE
  mnLL_f.initialize();
  #endif
  dataLikelihood.allocate("dataLikelihood");
  #ifndef NO_AD_INITIALIZE
  dataLikelihood.initialize();
  #endif
  mPrior.allocate("mPrior");
  #ifndef NO_AD_INITIALIZE
  mPrior.initialize();
  #endif
  ssbPrior.allocate("ssbPrior");
  #ifndef NO_AD_INITIALIZE
  ssbPrior.initialize();
  #endif
  selPrior.allocate("selPrior");
  #ifndef NO_AD_INITIALIZE
  selPrior.initialize();
  #endif
  recPrior.allocate("recPrior");
  #ifndef NO_AD_INITIALIZE
  recPrior.initialize();
  #endif
  fLike.allocate("fLike");
  #ifndef NO_AD_INITIALIZE
  fLike.initialize();
  #endif
  catLike.allocate("catLike");
  #ifndef NO_AD_INITIALIZE
  catLike.initialize();
  #endif
  dirCatPrior.allocate("dirCatPrior");
  #ifndef NO_AD_INITIALIZE
  dirCatPrior.initialize();
  #endif
  undirCatPrior.allocate("undirCatPrior");
  #ifndef NO_AD_INITIALIZE
  undirCatPrior.initialize();
  #endif
  hPrior.allocate("hPrior");
  #ifndef NO_AD_INITIALIZE
  hPrior.initialize();
  #endif
  pDirPrior.allocate("pDirPrior");
  #ifndef NO_AD_INITIALIZE
  pDirPrior.initialize();
  #endif
  repRatePrior.allocate("repRatePrior");
  #ifndef NO_AD_INITIALIZE
  repRatePrior.initialize();
  #endif
  qdevPrior.allocate(1,nIndexSeries,"qdevPrior");
  #ifndef NO_AD_INITIALIZE
    qdevPrior.initialize();
  #endif
  idxLam.allocate("idxLam");
  #ifndef NO_AD_INITIALIZE
  idxLam.initialize();
  #endif
  x.allocate("x");
  #ifndef NO_AD_INITIALIZE
  x.initialize();
  #endif
  lim1.allocate("lim1");
  #ifndef NO_AD_INITIALIZE
  lim1.initialize();
  #endif
  lim2.allocate("lim2");
  #ifndef NO_AD_INITIALIZE
  lim2.initialize();
  #endif
 int firstSelPrior = min(firstSelPrior_g);
 int lastSelPrior = max(lastSelPrior_g);
  p_alpha_gt.allocate(1,3,firstSelPrior,lastSelPrior,"p_alpha_gt");
  #ifndef NO_AD_INITIALIZE
    p_alpha_gt.initialize();
  #endif
  p_beta_gt.allocate(1,3,firstSelPrior,lastSelPrior,"p_beta_gt");
  #ifndef NO_AD_INITIALIZE
    p_beta_gt.initialize();
  #endif
  sd_alpha_gt.allocate(1,3,firstSelPrior,lastSelPrior,"sd_alpha_gt");
  #ifndef NO_AD_INITIALIZE
    sd_alpha_gt.initialize();
  #endif
  sd_beta_gt.allocate(1,3,firstSelPrior,lastSelPrior,"sd_beta_gt");
  #ifndef NO_AD_INITIALIZE
    sd_beta_gt.initialize();
  #endif
  iz_m.allocate(1,nAgeSeries,1,nT,minAge,maxAge,"iz_m");
  #ifndef NO_AD_INITIALIZE
    iz_m.initialize();
  #endif
  iz_f.allocate(1,nAgeSeries,1,nT,minAge,maxAge,"iz_f");
  #ifndef NO_AD_INITIALIZE
    iz_f.initialize();
  #endif
  pred_accBin_m.allocate(1,nAgeSeries,1,nT,"pred_accBin_m");
  #ifndef NO_AD_INITIALIZE
    pred_accBin_m.initialize();
  #endif
  pred_accBin_f.allocate(1,nAgeSeries,1,nT,"pred_accBin_f");
  #ifndef NO_AD_INITIALIZE
    pred_accBin_f.initialize();
  #endif
 cout << "Parameter Section Complete" << endl;
}

void model_parameters::set_runtime(void)
{
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  for ( int t=1;t<=landCatchMatrix.rowmax();t++ )
  {
    // katch: Rows are different fisheries.  Columns are tSteps.
    for ( int i=3;i<=landCatchMatrix.colmax();i++ )
    {
        Ctg(t,i-2) = landCatchMatrix(landCatchMatrix(t,2),i)/1000.;
      // converted catch units (t) to model units (1,000 t)
        obs_relCtg(t,i-2) = releaseMatrix(releaseMatrix(t,2),i)/1000.;      
    }
  }
  // Length-at-age
  lenAge_m = lInf_m+(L1_m-lInf_m)*exp( -vonK_m*(ages-1.));
  lenAge_f = lInf_f+(L1_m-lInf_f)*exp( -vonK_f*(ages-1.));
  // Weight-at-age
  wtAge_m  = wt_a(1)*pow(lenAge_m,wt_b(1));
  wtAge_f  = wt_a(2)*pow(lenAge_f,wt_b(2));
  
  // Maturity-at-age females
  matAge_f = 1./(1.+ exp( -log(19.)*(ages-aMat50)/(aMat95-aMat50)));
  // Proportion released-at-age function
  pR_gta_m.initialize(); pR_gta_f.initialize();
  pRet_m.initialize(); pRet_f.initialize();
  pRet_m = 1.; pRet_f = 1.;
  nSteps = 5;
  // Implement 55 cm comm. size limit in sizeLimYear
  for( int t=sizeLimYear; t<=nT; t++ )
  {
    for( int a=1; a<=plusGroupAge; a++ )
    {
      // Proportion retained-at-age by integrating the normal
      // pdf of length for each age and taking the proportion
      // above the size limit
      funnel_dvariable pL;
      funnel_dvariable pU;
      ad_begin_funnel();
      x      = lenAge_m(a);
      sigmaL = sigmaL_m;
      lim1 = 0.;
      lim2 = sizeLim;
      pL    = adromb( &model_parameters::fz, lim1, lim2, nSteps);
      lim1 = sizeLim;
      lim2 = lInf_m*(1. + 4.*sigmaL_m);
      pU    = adromb( &model_parameters::fz, lim1, lim2, nSteps);
      pRet_m(t,a) = value( pU/(pL+pU) );
      x    = lenAge_f(a);
      sigmaL = sigmaL_f;
      lim1 = 0.;
      lim2 = sizeLim;
      pL    = adromb( &model_parameters::fz, lim1, lim2, nSteps);
      lim1 = sizeLim;
      lim2 = lInf_f*(1. + 4.*sigmaL_f);
      pU    = adromb( &model_parameters::fz, lim1, lim2, nSteps);
      pRet_f(t,a) = value( pU/(pL+pU) );
    }
  }
  //cout << row(pRet_m,sizeLimYear+1) << endl;
  //cout << row(pR_gta_m(1),sizeLimYear+1) << endl;
  //cout << row(pR_gta_f(1),sizeLimYear+1) << endl;
  //exit(1);
  iz_m.initialize();
  iz_f.initialize();
  for( int i=1; i<=nAgeSeries; i++ )
  {
    iz_f(i) = 1.; iz_m(i) = 1.;
    for( int t=1; t<=nT; t++ )
    {
      // MALES--------------------------
      // check whether this year has data
      accBin_m(i,t) = 0.;
      if( firstAge_m(i,t) > 0 )
      {
        // vector of all age props for this year
        dvector obsPat_m  = row(ageObsProp_m(i),t)(firstAge_m(i,t),lastAge_m(i,t));
        for( int a=firstAge_m(i,t); a<=lastAge_m(i,t); a++ )
        {
          if( obsPat_m(a) <= prop_threshold )
          {
            // add obsPat to accBin
            accBin_m(i,t) += obsPat_m(a);
            // set indicator to 0 so (1) this bin is added to accumulator in predicted
            // and (2) multiplying obs by iz will eliminate these values in the 
            // likelihood and residual count
            iz_m(i,t,a) = 0.;
          }
        } // end a-loop
      }
      // FEMALES--------------------------
      // check whether this year has data
      accBin_f(i,t) = 0.;
      if( firstAge_f(i,t) > 0 )
      {
        // vector of all age props for this year
        dvector obsPat_f  = row(ageObsProp_f(i),t)(firstAge_f(i,t),lastAge_f(i,t));
        for( int a=firstAge_f(i,t); a<=lastAge_f(i,t); a++ )
        {
          if( obsPat_f(a) < prop_threshold )
          {
            // add obsPat to accBin
            accBin_f(i,t) += obsPat_f(a);
            // set indicator to 0 so (1) this bin is added to accumulator in predicted
            // and (2) multiplying obs by iz will eliminate these values in the 
            // likelihood and residual count
            iz_f(i,t,a) = 0;
          }
        } // end a-loop
      }
    }     // end t-loop
  }       // end i-loop
  
  neval = 0;
}

void model_parameters::userfunction(void)
{
  objFunc =0.0;
  // make sel prior matrix
  makeSelPriorMtx();
  if (verbose) cout << "makeSelPriorMtx() done..." << endl;
  // Calculate proportion directed
  calc_pDir_tg();
  if (verbose) cout << "calc_pDir_tg() done..." << endl;
  // Reporting rate functions
  calcRepRate();
  calcRepRatePrior();
  initModelParameters();
  if (verbose) cout << "initModelParameters done..." << endl;
  popInit_SSB0();
  if (verbose) cout << "popInit done..." << endl;
  popDynamics_MF();
  if (verbose) cout << "popDynamics done..." << endl;  
  calc_age_likelihood();
  if (verbose) cout << "age_likelihood = " <<ageLikelihood_m<<" "<<ageLikelihood_f<< endl;
  calc_index_likelihood();
  if (verbose) cout << "index_likelihood = " <<indexLikelihood<< endl;
  calc_release_likelihood();
  if (verbose) cout << "release_likelihood = " <<releaseLikelihood<<endl;
  if( active(recDevs) )
    calc_rec_prior();
  if (verbose) cout << "rec_prior = " <<recPrior<< endl;
  calc_sel_prior();
  if (verbose) cout << "sel_prior = " <<selPrior<< endl;
  if( !useBaranov )
  {
    calc_F_likelihood();
    calc_catch_likelihood();
    if (verbose) cout << "catLike = " <<catLike<< endl;
  }
  if (verbose) cout << "fLike = " <<fLike<< endl;
  if( active(log_M) )
    calc_M_prior();
  calc_pDirPrior();
  calc_SSB0_prior();
  for(int i = 1; i<= nIndexSeries; i++)
  {
    qdevPrior(i) = 0;
    if(active(log_qdev_gt(i)))
      for( int t = firstqdev_g(i); t <= lastqdev_g(i); t++)
        if( log_qdev_gt(i,t) != 0 )
          qdevPrior(i) += 0.5*pow( log_qdev_gt(i,t)/priorSD_qdev(i), 2 );
  } // end i loop
  objFunc  = 0.;  
  objFunc  = sum( elem_prod(idxLikeWeight,indexLikelihood ) );
  objFunc += sum( elem_prod(ageLikeWeight_m,ageLikelihood_m) ) 
          +  sum( elem_prod(ageLikeWeight_f,ageLikelihood_f) );
  objFunc += sum(releaseLikelihood);
  // Save data-only likelihood
  dataLikelihood = objFunc;
  // add priors to objFunc
  objFunc += fLike;
  objFunc += catLike;
  objFunc += mPrior + recPrior + selPrior;
  objFunc += 100.*ssbPrior;
  objFunc += hPrior;
  objFunc += sum(qdevPrior);
  objFunc += pDirPrior;
  objFunc += repRatePrior;
  if( last_phase() )
    objFunc *= objFunc_scale_last;
  else
    objFunc *= objFunc_scale;
  //cout << objFunc << endl;
  if (verbose) cout << "obj_func done..." << endl;
  // mcmc output when in mceval phase  
  if ( mceval_phase() )
  {
    ofstream mcout("mcout.dat", ios::app );
    // On first call, build header for MCMC file output 
    if ( mcmcHeader == 0 )
    {
      int i = 1;
      mcout << "SSB0 h M_m M_f";
      mcout << " projSSB projLegalB projSublegalB projDep legalHR.T sublegalHR.T";
      for ( i=1; i<=nT+1; i++ )
        mcout << " spawnB" << i;
      for ( i=firstRecDev; i<=lastRecDev; i++ )
        mcout << " recDev" << i;
      for ( i=1; i<=nIndexSeries; i++ )
        mcout << " tauIndex" << i;
      for ( i=1; i<=nAgeSeries; i++ )
        mcout << " tauAges_m" << i;
      for ( i=1; i<=nAgeSeries; i++ )
        mcout << " tauAges_f" << i;
      for ( i=1; i<=nFisheries; i++ )
        mcout << " tauReleases" << i;
      for ( i=1; i<=nIndexSeries; i++ )
        mcout << " qg" << i;      
      for( int g = 1; g <= nFisheries - 2; g++)
        for ( i=1; i<=nT; i++ )
          mcout << " F_g" << g <<"_t" << i;      
      for( int g = 1; g <= nFisheries; g++ )
        mcout << " alpha_g" << g;
      for( int g = 1; g <= nFisheries; g++ )
        mcout << " beta_g" << g;        
      mcout << endl;      
      mcmcHeader = 1;
    }
    mcout << SSB0 << " " << h << " " << M_m << " " << M_f << " ";
    mcout << SSBt(nT+1) << " " << legalBt(nT+1) << " " << subLegalBt(nT+1) << " "
          << SSBt(nT+1)/SSBt(1) << " " << legalHR(nT) << " " << subLegalHR(nT) << " " ;
    mcout << SSBt;
    mcout << recDevs;
    mcout << sqrt(tauSquareIndex);
    mcout << sqrt(tauSquareAges_m);
    mcout << sqrt(tauSquareAges_f);
    mcout << sqrt(tauSquareReleases);
    mcout << exp(lnq); 
    for( int g = 1; g <= nFisheries - 2; g++ )  
      mcout << column(Ftg,g) << " ";
    mcout << alpha_g1;
    mcout << beta_g1;
    mcout << endl; 
  }
}

void model_parameters::between_phases_calculations(void)
{
}

dvariable model_parameters::fz(const dvariable& z)
{
  double tmp1; double tmp;
  double mx;
  tmp1 = value( 1./sqrt(2.*3.141593*pow(x*sigmaL,2.)) );
  mx   = value(z);
  tmp  = value( mfexp( -0.5*pow(x-mx,2)/pow(x*sigmaL,2)) );
  //tmp  *= tmp1;
  return tmp;
}

void model_parameters::calc_SSB0_prior(void)
{
  ssbPrior = -log(1./SSB0);
}

void model_parameters::calc_M_prior(void)
{
  mPrior = 0.;
  mPrior += pow(M_m - prior_mean_M(1),2.)/pow(prior_sd_M(1),2) + 
            pow(M_f - prior_mean_M(2),2.)/pow(prior_sd_M(2),2);     
}

void model_parameters::calcRepRate(void)
{
  for( int g = 1; g <= nFisheries - 2; g ++ )
  {
    dvariable tmpLenRep50;
    dvariable tmpLenRep95;
    dvariable step;
    tmpLenRep50 = lenRep50(g) * mfexp(log_epsLenRep50_g(g));
    tmpLenRep95 = lenRep95(g) * mfexp(log_epsLenRep95_g(g));
    step=  tmpLenRep50 - tmpLenRep95 ;
    for( int a = 1; a <= plusGroupAge; a++)
    {
      reportRate_ga_m(g,a) = 1. / (1. + exp( -(log(19)/step) * ( tmpLenRep50 - lenAge_m(a) ) ));
      reportRate_ga_f(g,a) = 1. / (1. + exp( -(log(19)/step) * ( tmpLenRep50 - lenAge_f(a) ) ));
    }
  }
}

void model_parameters::calcRepRatePrior(void)
{
  repRatePrior = sum(pow( log_epsLenRep95_g,2 ) + pow( log_epsLenRep50_g,2 ) );
}

void model_parameters::calc_pDirPrior(void)
{
  pDirPrior = 0.;
  // pDirPrior is centred at -3 in first year (5% directed)
  for( int g = 1; g <= nFisheries-2; g++ )
  {
    if( (ph_logit_pDir(g) > 0) &  (logit_pDir_g(g) != 0) )
      pDirPrior += 0.5*pow(((logit_pDir_g(g))/priorSD_logit_pDir),2);
    if( (ph_log_selAlphaDirDev(g) > 0) & (log_selAlphaDirDev_g(g) != 0) )
      pDirPrior += 0.5 * pow((log_selAlphaDirDev_g(g))/priorSD_seldev_pDir,2.0);
    if( (ph_log_selBetaDirDev(g) > 0) & (log_selBetaDirDev_g(g) != 0) )
      pDirPrior += 0.5 * pow((log_selBetaDirDev_g(g))/priorSD_seldev_pDir,2.0);
  }
  //for( int t = pDir_t1 + 1; t <= pDir_t2; t++)
  //{
    //// Take difference of subsequent points
    //dvariable tmpDiff = logit_propDirected_t(t) - logit_propDirected_t(t-1);
    //// treat logit pDir as a random walk in logit space
    //if(active(logit_propDirected_t) & tmpDiff > 0)
      //pDirPrior += 0.5 * pow((tmpDiff)/priorSD_logit_pDir,2.0);
    //// And now add selAlphapDirDevs
    //if(active(log_selAlphaDirDev_t) & (log_selAlphaDirDev_t(t) > 0) )
      //pDirPrior += 0.5 * pow((log_selAlphaDirDev_t(t))/priorSD_seldev_pDir,2.0);
    //if(active(log_selBetaDirDev_t) & (log_selBetaDirDev_t(t) > 0))
      //pDirPrior += 0.5 * pow((log_selBetaDirDev_t(t))/priorSD_seldev_pDir,2.0);
  //}
}

void model_parameters::makeSelPriorMtx(void)
{
  //cout << "Pulling together sel priors" << endl;
  // Loop to pull together the sel priors
  for( int gIdx = 1; gIdx <= 3; gIdx ++)
    for( int tIdx = firstSelPrior_g(gIdx); tIdx <= lastSelPrior_g(gIdx); tIdx ++)
    {
      //cout << "Looping, tIdx == " << tIdx << endl;
      switch( gIdx )
      { 
        case 1: // Trap
          p_alpha_gt(gIdx,tIdx) = p_alpha_1t(tIdx);
          p_beta_gt(gIdx,tIdx) = p_beta_1t(tIdx);
          sd_alpha_gt(gIdx,tIdx) = sd_alpha_1t(tIdx);
          sd_beta_gt(gIdx,tIdx) = sd_beta_1t(tIdx);
        break;
        case 2: // LL 
          p_alpha_gt(gIdx,tIdx) = p_alpha_2t(tIdx);
          p_beta_gt(gIdx,tIdx) = p_beta_2t(tIdx);
          sd_alpha_gt(gIdx,tIdx) = sd_alpha_2t(tIdx);
          sd_beta_gt(gIdx,tIdx) = sd_beta_2t(tIdx);
        break;
        case 3: // Trawl
          p_alpha_gt(gIdx,tIdx) = p_alpha_3t(tIdx);
          p_beta_gt(gIdx,tIdx) = p_beta_3t(tIdx);
          sd_alpha_gt(gIdx,tIdx) = sd_alpha_3t(tIdx);
          sd_beta_gt(gIdx,tIdx) = sd_beta_3t(tIdx);
        break;
      }
    }
}

void model_parameters::calc_pDir_tg(void)
{
  // Computes the proportion directed for each gear type at
  // each time step
  pDir_tg=1.0;
  // We can improve this with a sub-vector of the pDir mtx
  for( int g = 1; g <= nFisheries-2; g++)
    if( pDirOn(g) == 1 )
      for( int t = pDir_t1(g); t <= pDir_t2(g); t++)
      {
        pDir_tg(t,g) = pDirScale / ( 1 + exp(-1. * logit_pDir_g(g)));
      }
  // cout << "pDir_tg = " << pDir_tg << endl;
}

void model_parameters::calc_pR_gta(void)
{
  // Computes the proportion of catch-at-age released due to combination
  // of size limit and high-grading. hg50 is the high-grading curve length
  // at 50% release
  for( int t=sizeLimYear; t<=nT; t++ )
  {
    for( int a=1; a<=plusGroupAge; a++ )
    {
      // Get pR_gta_m and ..._f here by adjusting for high-grade
      // function
      for( int g=1; g<=nFisheries-2; g++ )
      {
        pR_gta_m(g,t,a) = pRet_m(t,a);
        pR_gta_f(g,t,a) = pRet_f(t,a);
        if( useHighgrading(g) == 1)
        {
          pR_gta_m(g,t,a) *= 1./(1.+exp( -log(19.)*(lenAge_m(a)-exp(4+log_hg50(g))) ));
          pR_gta_f(g,t,a) *= 1./(1.+exp( -log(19.)*(lenAge_f(a)-exp(4+log_hg50(g))) ));
        }
      }
    }
  }
}

void model_parameters::calc_sel_gta(void)
{
  // Adding undirected selectivity, identified with
  // directed selectivity, for non-trawl fleets, in case
  // we want to add an undirected/bycatch portion of F
  // for other fleets later on. U/Undir suffix used to 
  // denote variables
  int g;
  dvar_vector tmp(1,plusGroupAge); // temp vars for orderly computations
  dvar_vector tmp_lenAge(1,plusGroupAge); // temp vars for orderly computations
  dvar_matrix tmp_alpha_gt(1,nFisheries,1,nT); // temp vars for orderly computations
  dvar_matrix tmp_beta_gt(1,nFisheries,1,nT); // temp vars for orderly computations
  dvar_matrix tmp_alphaU_gt(1,nFisheries,1,nT); // temp vars for orderly computations of undirected sel
  dvar_matrix tmp_betaU_gt(1,nFisheries,1,nT); // temp vars for orderly computations of undirected sel
  dvar_vector tmpUndir(1,plusGroupAge);
  dvar_vector log_alpha_t(1,nT); alpha_gt_m.initialize();
  dvar_vector log_beta_t(1,nT);  beta_gt_m.initialize();
  // Loop over fisheries
  for( int g=1; g<=nFisheries; g++ )
  {
    if( verbose ) cout << "Calculating sel for fleet " << g << endl;
    // Re-initialize temp vectors
    log_alpha_t.initialize();
    log_beta_t.initialize();
    // Loop over years
    for( int t=1; t<=nT; t++ )
    {
      if( t==1 )
      {
        // Initialize first value (in log-cm units)
        alpha_g1(g)     = exp(log_alpha_g(g));
        beta_g1(g)      = exp(log_beta_g(g));
        //if( selType(g) != 3 )
        alpha_g1(g) += L1_f;
        log_alpha_t(1) = log_alpha_g(g);
        log_beta_t(1)  = log_beta_g(g);        
      }
      else
      {
        // Let selectivity vary for years with prior info
        // in the commercial fleets
        if( g <= 3 )
        {
          if( (t>=firstSelPrior_g(g)) && (t<=lastSelPrior_g(g)) )
          {
            log_alpha_t(t) = log_alpha_t(t-1) + log_alpha_dev_gt(g,t);
            log_beta_t(t)  = log_beta_t(t-1) + log_beta_dev_gt(g,t);                        
          }
          else
          {
            log_alpha_t(t) = log_alpha_t(t-1);
            log_beta_t(t)  = log_beta_t(t-1);            
          }
        }
        else
        {
          // Params are constant
          log_alpha_t(t) = log_alpha_t(t-1);
          log_beta_t(t)  = log_beta_t(t-1);
        }  
      }
      // add block deviations for commercial fisheries
      if( g <= nFisheries - 2)
      {
        // cout << "In block dev calculation, g = " << g << endl;
        // Check if phz is positive, if so, add block devs
        if( ph_alpha_blockDev(g) > 0)
        {
          // cout << "phz positive for g = "<< g << endl;
          // add block deviations
          for(int b = 1; b <= nSelBlocks(g); b ++ )
          {
            // cout << "block b = "<< b << " of nB = " << nSelBlocks(g) << endl;
            // Check gear type
            // Trap
            if(g == 1)
            {
              if( t == selBlockStart_g1(b) )
              {
                log_alpha_t(t) += log_alpha_blockDev_g(g,b);
                log_beta_t(t) += log_beta_blockDev_g(g,b);
              }
              if( t == selBlockEnd_g1(b) + 1 )
              {
                log_alpha_t(t) -= log_alpha_blockDev_g(g,b);
                log_beta_t(t) -= log_beta_blockDev_g(g,b); 
              }
            }
            // Hook
            if(g == 2)
            {
              if( t == selBlockStart_g2(b) )
              {
                log_alpha_t(t) += log_alpha_blockDev_g(g,b);
                log_beta_t(t) += log_beta_blockDev_g(g,b);
              }
              if( t == selBlockEnd_g2(b) + 1 )
              {
                log_alpha_t(t) -= log_alpha_blockDev_g(g,b);
                log_beta_t(t) -= log_beta_blockDev_g(g,b); 
              }
            }
            // Trawl
            if(g == 3)
            {
              if( t == selBlockStart_g3(b) )
              {
                log_alpha_t(t) += log_alpha_blockDev_g(g,b);
                log_beta_t(t) += log_beta_blockDev_g(g,b);
              }
              if( t == selBlockEnd_g3(b) + 1 )
              {
                log_alpha_t(t) -= log_alpha_blockDev_g(g,b);
                log_beta_t(t) -= log_beta_blockDev_g(g,b); 
              }
            }
          }
        }
      }
      if( verbose) cout << "Finished with block dev calcs... " << endl;
      // Final selectivity parameters (now in cm units)
      alpha_gt_m(g,t) = mfexp( log_alpha_t(t) );
      alpha_gt_m(g,t) += L1_m;
      alphaUndir_gt_m(g,t) = alpha_gt_m(g,t);
      alpha_gt_f(g,t) = mfexp( log_alpha_t(t) + log_alpha_sexdev_g(g) );
      alpha_gt_f(g,t) += L1_f;
      alphaUndir_gt_f(g,t) = alpha_gt_f(g,t);
      switch( selType(g) )
      {
        case 1: // beta is asymptotic L95_step
          beta_gt_m(g,t)  = alpha_gt_m(g,t) - mfexp( log_beta_t(t) );
          beta_gt_f(g,t)  = alpha_gt_f(g,t) - mfexp( log_beta_t(t) + log_beta_sexdev_g(g) );
        break;
        case 2: // beta is normal sd
          beta_gt_m(g,t)  = mfexp( log_beta_t(t) );
          beta_gt_f(g,t)  = mfexp( log_beta_t(t) + log_beta_sexdev_g(g) );
        break;
        case 3: // beta is gamma scale parameter
          beta_gt_m(g,t)  = mfexp( log_beta_t(t) );
          beta_gt_f(g,t)  = mfexp( log_beta_t(t) + log_beta_sexdev_g(g) );
        break;
      }
      // save beta as betaUndir, might need updating later
      betaUndir_gt_m(g,t) = beta_gt_m(g,t);
      betaUndir_gt_f(g,t) = beta_gt_f(g,t);
      if( verbose ) cout << "Finished with sex dev calcs for fleet " << g << endl;
      if( g <= nFisheries-2 )
        if( (t >= pDir_t1(g)) & (t <= pDir_t2(g)) )
        {
          if( verbose ) cout << "Computing undirected sel parameters for fleet " << g << endl;
          // Modify alpha parameter
          if( (pDirOn(g) == 1) )
          {
            alphaUndir_gt_m(g,t) *= exp(log_selAlphaDirDev_g(g) + selAlphaDirDev_base(g));
            alphaUndir_gt_f(g,t) *= exp(log_selAlphaDirDev_g(g) + selAlphaDirDev_base(g));
            // Modify beta parameter
            betaUndir_gt_m(g,t) *= exp(log_selBetaDirDev_g(g));
            betaUndir_gt_f(g,t) *= exp(log_selBetaDirDev_g(g));
          }
        }
      // Loops over sexes and age, compute selectivity      
      for( int l=1; l<=2; l++ )
      {      
        if( l==1 ) 
        {
          tmp_lenAge    = lenAge_m;
          tmp_alpha_gt  = alpha_gt_m;
          tmp_beta_gt   = beta_gt_m;
          tmp_alphaU_gt = alphaUndir_gt_m;
          tmp_betaU_gt  = betaUndir_gt_m;
        }
        else
        {
          tmp_lenAge    = lenAge_f;
          tmp_alpha_gt  = alpha_gt_f;
          tmp_beta_gt   = beta_gt_f;
          tmp_alphaU_gt = alphaUndir_gt_f;
          tmp_betaU_gt  = betaUndir_gt_f;
        }
        tmp.initialize();
        tmpUndir.initialize();
        for( int a=1; a<=plusGroupAge; a++ )
        {
          switch( selType(g) )
          {
            if( verbose ) cout << "Computing undirected sel curves for fleet " << g << endl;
            case 1: // asymptotic
              tmp(a) = exp( (-1.)*log(19.)*(tmp_lenAge(a)-tmp_beta_gt(g,t))/(tmp_alpha_gt(g,t) - tmp_beta_gt(g,t) ) );
              tmp(a) = 1./(1.+tmp(a));
              tmpUndir(a) = exp( (-1.)*log(19.)*(tmp_lenAge(a)-tmp_betaU_gt(g,t))/(tmp_alphaU_gt(g,t) - tmp_betaU_gt(g,t) ) );
              tmpUndir(a) = 1/(1 + tmpUndir(a));
            break;
            case 2: // normal
              tmp(a) = exp( (-0.5)*pow(tmp_lenAge(a)-tmp_alpha_gt(g,t),2.)/pow(tmp_beta_gt(g,t),2.) );               
              tmpUndir(a) = exp( (-0.5)*pow(tmp_lenAge(a)-tmp_alphaU_gt(g,t),2.)/pow(tmp_betaU_gt(g,t),2.) );               
            break;
            case 3: // gamma
              tmp(a)      = pow(tmp_lenAge(a),(tmp_alpha_gt(g,t) - 1)) * mfexp( -tmp_lenAge(a) / tmp_beta_gt(g,t));
              tmpUndir(a) = pow(tmp_lenAge(a),(tmp_alphaU_gt(g,t) - 1)) * mfexp( -tmp_lenAge(a) / tmp_betaU_gt(g,t));
          }
        } // end a-loop
        for( int a=1; a<=plusGroupAge; a++ )
        {
          if( g == 4 & a <= 2 )
          {
            tmp(a) = 0;
            tmpUndir(a) = 0;
          }
          // Adding undirected selectivity, identified with
          // directed selectivity for non-trawl fleets, in case
          // we want to add an undirected/bycatch portion of F
          // for other fleets later on
          // Normalise sel
          if( l==1 )
          {
            sel_gta_m(g,t,a) = tmp(a)/max(tmp);
            selUndir_gta_m(g,t,a) = tmpUndir(a)/max(tmpUndir);
          }
          else
          {
            sel_gta_f(g,t,a) = tmp(a)/max(tmp);
            selUndir_gta_f(g,t,a) = tmpUndir(a)/max(tmpUndir);
          }
          //if( g==1) cout << "t = "<<t<< " a = "<<a<<" sel_gta_f = "<< sel_gta_f(g,t,a) << endl;
        }
      }   // end l-loop
    }     // end t-loop
  }       // end g-loop
  //exit(1);
}

void model_parameters::calc_sel_prior(void)
{
  // Priors needed for
  // calc prior for years 32-48
  int g;
  selPrior = 0.;
  for( int g=1; g<=nFisheries-2; g++ )
  {
    for( int t=firstSelPrior_g(g); t<=lastSelPrior_g(g); t++ )
    {
      if(active(log_alpha_g(g)) )
      {
        // Males
        selPrior += 0.5*pow( alpha_gt_m(g,t)-p_alpha_gt(g,t),2. )/pow(sd_alpha_gt(g,t),2.);  
        // Females
        selPrior += 0.5*pow( alpha_gt_f(g,t)-p_alpha_gt(g,t),2. )/pow(sd_alpha_gt(g,t),2.);
      }
      if( active(log_beta_g(g))  )
      {
        // Males
        selPrior += 0.5*pow( beta_gt_m(g,t)-p_beta_gt(g,t),2. )/pow(sd_beta_gt(g,t),2.);  
        // Females
        selPrior += 0.5*pow( beta_gt_f(g,t)-p_beta_gt(g,t),2. )/pow(sd_beta_gt(g,t),2.);
      }
      // prior on selectivity std deviation deviations...
      if( active(log_alpha_dev_gt(g)) & (log_alpha_dev_gt(g,t) != 0) )
        selPrior += 0.5*pow(log_alpha_dev_gt(g,t)/priorSD_alpha_dev,2.);
      if( active(log_beta_dev_gt(g)) & (log_beta_dev_gt(g,t) != 0) )
        selPrior += 0.5*pow(log_beta_dev_gt(g,t)/priorSD_beta_dev,2.);
    }
    // Prior on sex-based selectivity devs
    if( (active(log_beta_sexdev_g(g))) & (log_beta_sexdev_g(g) != 0) )
      selPrior += 0.5*pow(log_beta_sexdev_g(g)/priorSD_beta_sexdev,2.);
    if( (active(log_beta_sexdev_g(g))) & (log_alpha_sexdev_g(g) != 0) )
      selPrior += 0.5*pow(log_alpha_sexdev_g(g)/priorSD_alpha_sexdev,2.);
    // Add block deviations to prior, with .1 sd
    if( active(log_alpha_blockDev_g(g) ) )
      selPrior  += 0.5*norm2(log_alpha_blockDev_g(g))/pow(blockDevPriorSd(g),2);
    if( active(log_beta_blockDev_g(g) ) )
      selPrior  += 0.5*norm2(log_beta_blockDev_g(g))/pow(blockDevPriorSd(g),2);
  }    
}

void model_parameters::calc_sel_prior2(void)
{
  // Priors needed for
  // calc prior for years 32-48
  int g;
  selPrior = 0.;
  for( int g=1; g<=nFisheries-2; g++ )
  {
    for( int t=firstSelPrior_g(g); t<=lastSelPrior_g(g); t++ )
    {
      if( active(log_alpha_g(g)) )
        selPrior += 0.5*pow( alpha_gt_m(g,t)-p_alpha_gt(g,t),2.)/pow(sd_alpha_gt(g,t),2.);
      if( active(log_beta_g(g)) )
        selPrior += 0.5*pow( beta_gt_m(g,t)-p_beta_gt(g,t),2.)/pow(sd_beta_gt(g,t),2.);
    }   
  }    
}

void model_parameters::calc_sel_prior3(void)
{
  // Priors for constant selectivity
  // in all three comm fisheries
  int g;
  selPrior = 0.;
  for( int g=1; g<=nFisheries-2; g++ )
  {
    if( active(log_alpha_g(g)) )
      selPrior += 0.5*pow( log(alpha_gt_m(g,1))-log(p_alpha_gt(g,firstSelPrior_g(g))),2. )/pow(sigma_alpha(g),2.);
    if( active(log_beta_g(g)) )
      selPrior += 0.5*pow( log(beta_gt_m(g,1))-log(p_beta_gt(g,firstSelPrior_g(g))),2. )/pow(sigma_beta(g),2.);
    // Add block deviations to prior, with .1 sd
    if( active(log_alpha_blockDev_g(g) ) )
      selPrior  += 0.5*norm2(log_alpha_blockDev_g(g))/pow(blockDevPriorSd(g),2);
    if( active(log_beta_blockDev_g(g) ) )
      selPrior  += 0.5*norm2(log_beta_blockDev_g(g))/pow(blockDevPriorSd(g),2);
  }    
  //exit(1);
}

void model_parameters::calc_rec_prior(void)
{
  recPrior  = 0.5*norm2(recDevs)/pow(sigma_R,2.);
}

void model_parameters::calc_h_prior(void)
{
  hPrior = 0.;
  hPrior = (prior_h(1)-1.)*log(h) + (prior_h(2)-1.)*log(1.-h);
  hPrior *= -1.;
}

void model_parameters::calc_release_likelihood(void)
{
  // Function computes likelihoods for at-sea releases for each gear type
  // using data from 1996+ (trawl) and 2006+ (trap, LL)
  dvariable res; validObs_rel.initialize();
  tauSquareReleases.initialize(); releaseLikelihood.initialize();
  int tMax = min(nT,likeYear);
  // cout << "Inside rel like" << endl;
  for( int g=1; g<=nFisheries-2; g++ )
  {
      ssR = 0.;
      for( int t=1; t<=tMax; t++ )
      {
        if( obs_relCtg(t,g) > 0. )
        {
          // Check if using the trawl fishery, or 
          //
          if(t < relDevInitYr(g) | t > relDevLastYr(g) )
          {
            res = log(obs_relCtg(t,g)) - log(relCtg(t,g));
            validObs_rel(g) +=1;
            ssR             += pow(res,2.);
          }
        } 
      }   
      if(validObs_rel(g) > 0)
      {
        tauSquareReleases(g) = ssR/validObs_rel(g);
        releaseLikelihood(g) = 0.5*(validObs_rel(g)-1.)*log(tauSquareReleases(g));
        // Added by SDNJ: an exponential prior on 
        // the release observation variance, attempting
        // to penalise bad fits - might have to play with
        // a weighting here.
        releaseLikelihood(g) += tauRelLambda(g) * sqrt(tauSquareReleases(g));
      }
      // Try adding a really stiff penalty for the residuals in the last 2 years
      // in the trawl fishery - this will be added to the likelihood in place
      // of the conditional normal MLEs above
      // cout << "Adding trawl release dev penalty in last 2 years" << endl;
      if( relDevInitYr(g) <= tMax )
      {
        for( int t = relDevInitYr(g); t <= relDevLastYr(g); t++)
        { 
          res = (log(obs_relCtg(t,g)) - log(relCtg(t,g)));
          if(res > 0)
            releaseLikelihood(g) += 0.5 * pow(res/relDevPenalty,2);
        }
      }
  }
  // cout << "Completed rel like" << endl;
}

void model_parameters::calc_index_likelihood(void)
{
  // Initialize likelihood function terms.
  int i,t,g;
  indexLikelihood.initialize(); validObs_idx.initialize();
  ss.initialize();
  dvariable zSum;       // cumulative residual
  dvar_vector z(1,nT);  // residuals
  // Calculations for MLEs of q.
  for ( i=1;i<=nIndexSeries;i++ )
  {
    z.initialize(); zSum=0.0; 
    for ( t=1;t<=likeYear;t++ )
    {
      if ( idxSeries(i,t) > 0. )
      { 
        g = idxIndex(i);
        z[t] = log(idxSeries(i,t)) - log( expBit(g,t)  );
        validObs_idx(i) += 1;
      }
    }
    // MLEe lnq 
    for( int t = 1; t <= nT; t++)
      lnq_it(i,t) = lnq(i);    
    // Sum-of-squares function.
    for( t=1; t<=likeYear; t++ )
    {
      // Add residual if asked for
      if( (t >= firstqdev_g(i)) & (t <= lastqdev_g(i)) & active(log_qdev_gt(i)))
        lnq_it(i,t) = lnq_it(i,t-1) + log_qdev_gt(i,t);
      if ( idxSeries(i,t) > 0. )
      {
        ss[i] += pow(z(t)-lnq_it(i,t),2.);        
      }
    }
    // Objective function value for this series.
    tauSquareIndex(i) = ss(i)/(validObs_idx(i)-1.); // index variance MLE
    indexLikelihood(i) = (validObs_idx(i)-1.)*log(tauSquareIndex(i)); // concentrated likelihood
    // Add a Jeffreys prior to sqrt(tauSquareIndex)
    indexLikelihood(i) += tauIndexLambda(i) * sqrt(tauSquareIndex(i));
  }
  //cout << "lnq_it = " << lnq_it << endl;
  //exit(1);
}

void model_parameters::calc_age_likelihood(void)
{
  dvariable meanDiff_m; dvariable meanDiff_f;
  //dvariable nYearsAges;
  dvariable etaSumSq_m; dvariable etaSumSq_f;
  n_residuals_m = 0; n_residuals_f = 0;  
  tauSquareAges_m.initialize(); tauSquareAges_f.initialize(); 
  ageLikelihood_m.initialize(); ageLikelihood_f.initialize();
  mnLL_m = 0.; mnLL_f = 0.;
  // If likeYear < nT, then use likeYear for retrospective
  int tMax = min(nT,likeYear);
  int g; 
  if(verbose == 1)
    cout << "Inside age likelihood" << endl;
  // Calculate predicted age-proportions for each index gear.
  for( int i=1;i<=nAgeSeries;i++ )
  {
    if(verbose == 1)
      cout << "Age series  = " << i << endl;
    g = fisheryAgeID(i);
    // males
    //cout<< "firstAge_m(i) = " <<  firstAge_m(i) << endl;
    etaSumSq_m = 0.; n_residuals_m = 0;
    for( int t=1; t<=nT; t++ )
    {
      if(verbose == 1)
        cout << "Time t = " << t << endl;
      // check whether this year has data
      dvariable accRes = 0.;
      if( firstAge_m(i,t) > 0 )
      {
        // extract temp vector of all age props for this year
        dvar_vector obsPat  = row(ageObsProp_m(i),t)(firstAge_m(i,t),lastAge_m(i,t));
        // cout << "t = " << t << endl;
        // cout << "firstAge_m(i,t) = " <<  firstAge_m(i,t) << endl;
        // cout << "lastAge_m(i,t) = " <<  lastAge_m(i,t) << endl;
        // cout << "obsPat = " << obsPat << endl;
        // NOTE: this needs to be summed up
        dvar_vector predPat = row(predObsProp_m(i),t);
        //dvar_vector predPat = row(uCgta_m(i),t)(firstAge_m(i,t),lastAge_m(i,t));
        pred_accBin_m(i,t) = 0.;
        for( int a=firstAge_m(i,t); a<=lastAge_m(i,t); a++ )
        {
          if( iz_m(i,t,a) == 0. )
          {
            // add predPat to accBin
            pred_accBin_m(i,t) += predPat(a);
            predPat(a) = 0.;
          }
        }
        // Likelihood calcs for males-------------------------
        // Residual function for this year
        int iRes = 0; dvariable sumRes = 0.;
        dvar_vector res(firstAge_m(i,t),lastAge_m(i,t));
        // Get residuals by age
        for( int a=firstAge_m(i,t); a<=lastAge_m(i,t); a++ )
        {
          if( iz_m(i,t,a) > 0. ) // Data exist and > prop_threshold
          {
            res(a)  = log(obsPat(a)) - log(predPat(a)); // residual
            sumRes   += res(a);
            iRes     += 1.;
          }
        }
        // If the accumulator bin was used, then get residual for that.
        if( pred_accBin_m(i,t) > 0. )
        {
          accRes = log(pred_accBin_m(i,t));
          sumRes += accRes; // add residual to total
          iRes   += 1.;     // add one more observation
        }
        // mean residual, SSQ, and n for MVN
        meanDiff_m     = sumRes/iRes;
        etaSumSq_m    += norm2( res - meanDiff_m );
        if( accRes > 0. ) 
          etaSumSq_m += pow(accRes - meanDiff_m,2.); // sumSq
        n_residuals_m += iRes;          // accumulate total n residuals    
        // cout << "etaSumSq_m = " << etaSumSq_m <<  endl;  
        // cout << "iRes = " << iRes <<  endl;  
        // cout << "sumRes = " << sumRes <<  endl;  
      }
    }
    if( verbose )
      cout << "Calculating tauSquareAges for males in gearid" << i << endl;
    // Final MLE of variance and likelihood
    tauSquareAges_m(i) = etaSumSq_m/n_residuals_m;
    ageLikelihood_m(i) = n_residuals_m*log(tauSquareAges_m(i))/2;
    ageLikelihood_m(i) += tauAgeLambda(i)*tauSquareAges_m(i);
    // cout << "i = " << i << endl;
    // cout << "etaSumSq_m = " << etaSumSq_m <<  endl;
    // cout << "n_residuals_m = " << n_residuals_m <<  endl;
    // exit(1);
    // females
    etaSumSq_f = 0.; n_residuals_f = 0;
    for( int t=1; t<=nT; t++ )
    {
      if(verbose)
        cout << "Female age observations in gear" << i << " at time " <<  t << endl;
      // check whether this year has data
      dvariable accRes = 0.;
      if( firstAge_f(i,t) > 0 )
      {
        // extract temp vector of all age props for this year
        dvar_vector obsPat  = row(ageObsProp_f(i),t)(firstAge_f(i,t),lastAge_f(i,t));
        // NOTE: this needs to be summed up
        dvar_vector predPat = row(predObsProp_f(i),t);
        //dvar_vector predPat = row(uCgta_f(i),t)(firstAge_f(i,t),lastAge_f(i,t));
        pred_accBin_f(i,t) = 0.;
        for( int a=firstAge_f(i,t); a<=lastAge_f(i,t); a++ )
        {
          if( iz_f(i,t,a) == 0. )
          {
            // add predPat to accBin
            pred_accBin_f(i,t) += predPat(a);
            predPat(a) = 0.;
          }
        }
        // Residual function for this year
        int iRes = 0; dvariable sumRes = 0.;
        dvar_vector res(firstAge_f(i,t),lastAge_f(i,t));
        // Get residuals by age in case some obsPat == 0
        for( int a=firstAge_f(i,t); a<=lastAge_f(i,t); a++ )
        {
          if( iz_f(i,t,a) > 0. ) // Data exist and > prop_threshold
          {
            res(a)  = log(obsPat(a)) - log(predPat(a)); // residual
            sumRes   += res(a);
            iRes     += 1.;
          }            
        }
        //cout<< "t = " << t << endl;
        //cout<< "res_f = " << res << endl;
        // If the accumulator bin was used, then get residual for that.
        if( pred_accBin_f(i,t) > 0. )
        {
          accRes = log(pred_accBin_f(i,t));
          sumRes += accRes; // add residual to total
          iRes   += 1.;     // add one more observation
        }
        // mean residual, SSQ, and n for MVN
        meanDiff_f     = sumRes/iRes;
        etaSumSq_f    += norm2( res - meanDiff_f );
        if( accRes > 0. ) 
          etaSumSq_f += pow(accRes - meanDiff_f,2.); // sumSq
        n_residuals_f += iRes;                    // total n residuals    
      }
    }
    if( verbose )
      cout << "Calculating tauSquareAges for females" << endl;
    // MLE of variance in female age-proportions.
    tauSquareAges_f(i) = etaSumSq_f/n_residuals_f;
    ageLikelihood_f(i) = n_residuals_f*log(tauSquareAges_f(i))/2;
    ageLikelihood_f(i) += tauAgeLambda(i)*tauSquareAges_f(i);
    //cout<< "i = " << i << endl;
    //cout<< "ageLike_f(i) = " << ageLikelihood_f(i) << endl;
  }
}

void model_parameters::calc_F_likelihood(void)
{
  fLike = 0.;
    /*dvariable ln_det;
    dvariable sign;
    dvar_vector rF  = column(log_Fdevs_gt,t)-row_means_F;
    dvar_vector tmp = solve(sigmaF,rF,ln_det,sign);
    fLike += (-.5)*ln_det-rF*tmp;
    */
    for( int g=1; g<=nFisheries-2; g++ )
    {
      for( int t=firstYearCatch(g)+1; t<=nT;t++ )
        fLike += 0.5*log_Fdevs_gt(g,t)*log_Fdevs_gt(g,t)/pow(priorSD_F,2.);
    /*if( last_phase() )
    {
      fLike += 1.*pow( mean(column(Ftg,g))-0.1,2. );
    }
    else
      fLike += 10000.*pow( mean(column(Ftg,g))-0.1,2. );
    */
    }
}

void model_parameters::calc_catch_likelihood(void)
{
  catLike = 0.;
  dvariable tmp;
  for( int t=1; t<=nT; t++ )
    for( int g=1; g<=nFisheries-2; g++ )
    {
      if( Ctg(t,g) > 0)
      {
        tmp = log(pCtg(t,g)) - log(Ctg(t,g));
        catLike += 0.5*pow(tmp/catError,2.); 
      }
    }
}

void model_parameters::initModelParameters(void)
{
  // Set likelihoods/priors to 0
  fLike = 0.; catLike = 0.; hPrior = 0.; ssbPrior = 0.;
  recPrior = 0.; selPrior = 0.;
  // Revise this order to match PARAMETER_SECTION
  // Natural mortality is constant over age and time
  // cout << "log_M = " <<  log_M << endl;
  M_m = mfexp( log_M(1) );
  M_f = mfexp( log_M(2) );
  // correlation in Fs
  rho = 1./(1.+exp(-logit_rho));
  // cout << "logit_rho = " <<  logit_rho << endl;
  // Random walk in Fs
  Ftg.initialize();
  if( !useBaranov )
  {
    // Ftg parameterized as random walk with estimated parameters
    //  log_Fg1,  log_Fg2, log_Fg3, log_Fdevs_1t,...3t
    log_Ftg.initialize(); 
    log_Fdevs_gt.initialize();
    // Fs in first year
    for(int g=1; g<=nFisheries-2; g++)
    {
      log_Ftg(firstYearCatch(g),g) = log_Fg1(g);
      // cout << "log_Fg1 = " << log_Fg1 << endl;
      Ftg(firstYearCatch(g),g)     = mfexp( log_Ftg(firstYearCatch(g),g) );
      switch( g )
      {
        case 1: // trap
          log_Fdevs_gt(g)(firstYearCatch(g)+1,nT) = log_Fdevs_1t;
        break;
        case 2: // normal
          log_Fdevs_gt(g)(firstYearCatch(g)+1,nT) = log_Fdevs_2t;
        break;
        case 3: // normal
          log_Fdevs_gt(g)(firstYearCatch(g)+1,nT) = log_Fdevs_3t;
        break;
      }
      for( int t=firstYearCatch(g)+1; t<=nT; t++ )
      {
        log_Ftg(t,g) = log_Ftg(t-1,g) + log_Fdevs_gt(g,t);
        Ftg(t,g)     = mfexp( log_Ftg(t,g) );
      }
    }
    /*// Fill ranwalk for remaining years
    for( int t=2; t<=nT; t++ )
    {
      for(int g=1; g<=nFisheries-2; g++)
      {
        log_Ftg(t,g) = log_Ftg(t-1,g) + log_Fdevs_gt(g,t);
        Ftg(t,g)     = mfexp( log_Ftg(t,g) );
      }
    }*/
    int nFs = log_Fdevs_gt.colmax()-log_Fdevs_gt.colmin();
    //cout << "rowmeans = " << rowsum(log_Fdevs_gt)/nFs << endl;
    row_means_F = rowsum(log_Fdevs_gt)/nFs;
  }
  /*ch_sigmaF.initialize();
  int ii=1;
  for (int i=1;i<=3;i++)
    for (int j=1;j<=i;j++)
      ch_sigmaF(i,j)=sigmaF_coeff(ii++);
  sigmaF = ch_sigmaF*trans(ch_sigmaF);
  for( int i=1; i<=3; i++ )
    sigmaF(i,i) += 0.001;
  */
  // get selectivity by gear and year
  calc_sel_gta();
  if( verbose ) cout << "calc_sel_gta() done..." << endl;
  // get proportion of catch-at-age released at sea due to
  // combination of size limit and high-grading 
  calc_pR_gta();
  //cout << "Done calc_sel_gta_mf " << endl;
  // get unfished spawning biomass per recruit for setting
  // up unfished recruitment and SR relationship
  calcPhi();
  // Average recruitment may be leading parameter
  avgR = mfexp( log_avgR );
  // Beverton-Holt steepness implementing 0.25-0.99 bounds
  // at logit_h=0, the default h=0.75
  h = 0.25 + (0.99-0.25)*mfexp(logit_h)/(1.+mfexp(logit_h));
  if( active(logit_h) )
    calc_h_prior();
  // Unfished SSB
  SSB0 = mfexp( log_SSB0 );
  // Unfished recruitment.
  R0 = SSB0/phiSSB;
  // Beverton-Holt a parameter.
  bh_a = 4.*h*R0/(SSB0*(1.-h));
  // Beverton-Holt b parameter.
  bh_b = (5.*h-1.)/(SSB0*(1.-h));
  // recruitment deviations
  omega_t.initialize(); // set all to 0
  // fill in the ones we are estimating
  for( int t=firstRecDev; t<=lastRecDev; t++ )
  {
    omega_t(t) = recDevs(t);
  } 
  // Overwrite recruitment deviations
  // for those years where they are fixed
  for( int i = 1; i <= nFixedRecDev; i ++)
    omega_t(tFixedRecDev(i)) = logFixedRecDev(i);
  //cout << R0 << endl;
}

void model_parameters::calcPhi(void)
{
  // Calculating unfished SSB-per-recruit. Borrowing
  // Nta and Bta for now. Initializing with 0.5 recruit
  // to get through calcs.
  Nta_m.initialize();SSBt.initialize();
  Nta_f.initialize();Bta_m.initialize();
  Bta_f.initialize();
  Nta_f(1,1) = 1.;
  for( int a=2; a<=(plusGroupAge-1); a++ )
  {
    Nta_f(1,a) = Nta_f(1,a-1)*mfexp( -M_f );
  }
  Nta_f(1,plusGroupAge) = Nta_f(1,plusGroupAge-1)*mfexp(-M_f)/(1.-mfexp(-M_f));
  Bta_f(1)(1,plusGroupAge) = elem_prod( row(Nta_f,1), wtAge_f );
  // Unfished spawning biomass-per-recruit
  phiSSB = sum( elem_prod( row(Bta_f,1),matAge_f ) );
  //cout << phiSSB << endl;
}

void model_parameters::popInit_SSB0(void)
{
  Bta_m.initialize();Bta_f.initialize(); SSBt.initialize();
  expBit.initialize(); legalBt.initialize(); subLegalBt.initialize();
  relCtg.initialize(); landedCtg.initialize(); 
  undirLandCtg.initialize(); undirRelCtg.initialize();
  dirLandCtg.initialize();   dirRelCtg.initialize();
  pCtg.initialize(); Nta_f.initialize(); Nta_m.initialize();
  //---------------------
  // Initialize from unfished
  Nta_m(1,1) = R0*mfexp( omega_t(1)  );    
  Nta_f(1,1) = Nta_m(1,1);        
  for( int a=2; a<=plusGroupAge-1; a++ )
  {
    Nta_m(1)(a) = Nta_m(1,a-1)*mfexp( -M_m );
    Nta_f(1)(a) = Nta_f(1,a-1)*mfexp( -M_f );
  }
  Nta_m(1)(plusGroupAge) = Nta_m(1)(plusGroupAge-1)*exp(-M_m);
  Nta_m(1)(plusGroupAge)/= 1.-exp(-M_m);
  Nta_f(1)(plusGroupAge) = Nta_f(1)(plusGroupAge-1)*exp(-M_f);
  Nta_f(1)(plusGroupAge)/= 1.-exp(-M_f);
  //cout << sum(row(Nta_f,1)) << endl;
  //exit(1);  
  // Biomass-at-age and SSB: SSBt(1) should match the value
  // at the end of calcPhi()
  Bta_m(1)(1,plusGroupAge) = elem_prod( row(Nta_m,1),wtAge_m );
  Bta_f(1)(1,plusGroupAge) = elem_prod( row(Nta_f,1),wtAge_f );
  SSBt(1)  = sum( elem_prod( row(Bta_f,1),matAge_f ) );
  Rt(1)    = Nta_f(1,1) + Nta_m(1,1);
  // Get total mortality rate
  Zta_m.initialize(); Zta_f.initialize();
  baranov_time = 1;
  if( useBaranov ) // solve catch equation to get Z each year
  {
    solveBaranov_multi_fleet();
  }
  else
  {
    // Fs are parameters: get all Zt
    for( int t=1; t<=nT; t++ )
    {
      // Initialize Z with only M 
      Zta_m(t)(1,plusGroupAge) = M_m; 
      Zta_f(t)(1,plusGroupAge) = M_f;
      // Add Fs
      for( int g=1; g<=nFisheries-2; g++ )
      {
        // subset the proportions released at age
        for( int a = 1; a <= plusGroupAge; a++)
        {
          // Directed portion of fishing mortality
          Zta_m(t,a) += pDir_tg(t,g) * sel_gta_m(g,t,a) * Ftg(t,g) * (pR_gta_m(g,t,a) + (1 - pR_gta_m(g,t,a)) * dM(g));
          Zta_f(t,a) += pDir_tg(t,g) * sel_gta_f(g,t,a) * Ftg(t,g) * (pR_gta_f(g,t,a) + (1 - pR_gta_f(g,t,a)) * dM(g));
          // Add undirected portion of fishing mortality
          Zta_m(t,a) += (1 - pDir_tg(t,g)) * selUndir_gta_m(g,t,a) * Ftg(t,g) * (pR_gta_m(g,t,a) + (1 - pR_gta_m(g,t,a)) * dM(g));
          Zta_f(t,a) += (1 - pDir_tg(t,g)) * selUndir_gta_f(g,t,a) * Ftg(t,g) * (pR_gta_f(g,t,a) + (1 - pR_gta_f(g,t,a)) * dM(g));
        }
        // cout << "t = "<<t<<" g = "<<g<< endl;
        // cout << "sel_gta(g,t) = " << row(sel_gta_f(g),t) << endl;
        //Zta_m(t)(1,plusGroupAge) += row(sel_gta_m(g),t)*Ftg(t,g);
        //Zta_f(t)(1,plusGroupAge) += row(sel_gta_f(g),t)*Ftg(t,g);
      }
    }
  // exit(1);
  }
  //-------------------------
  // End of Biomass dynamics
  //-------------------------
  //-------------------------------------------------
  // Observation models for indices and length comps
  calc_observation_models();
  //-------------------------------------------------
}

void model_parameters::popDynamics_MF(void)
{
  int g;
  //---------------------
  // Biomass dynamics
  //---------------------
  // Population abundances-at-age were initialized 
  // in popInit_SSB0 for year 1
  for( int t=2; t<=nT; t++ )
  {
    // Age-1 recruitment: m=male, f=female 
    // assumes 50:50 sex ratio at birth
    Nta_m(t,1) = bh_a*SSBt(t-1)*mfexp( omega_t(t) )/(1.+bh_b*SSBt(t-1)); 
    Nta_f(t,1) = Nta_m(t,1);             
    // Age-structured model equations: total mortality Zta calculated in popInit_SSB0
    // for all years by adding Ftg and M
    for( int a=2; a<=plusGroupAge-1; a++ )
    {
      Nta_m(t,a) = Nta_m(t-1,a-1)*mfexp( -Zta_m(t-1,a-1) );
      Nta_f(t,a) = Nta_f(t-1,a-1)*mfexp( -Zta_f(t-1,a-1) );
    }
    // The plus group:
    // Males
    Nta_m(t,plusGroupAge) = Nta_m(t-1,plusGroupAge-1)*mfexp( -Zta_m(t-1,plusGroupAge-1) ) + 
                            Nta_m(t-1,plusGroupAge)*mfexp(   -Zta_m(t-1,plusGroupAge)    );
    // Females
    Nta_f(t,plusGroupAge) = Nta_f(t-1,plusGroupAge-1)*mfexp( -Zta_f(t-1,plusGroupAge-1) ) + 
                            Nta_f(t-1,plusGroupAge)*mfexp(   -Zta_f(t-1,plusGroupAge)    );
    // Biomass-at-age and SSB
    Bta_m(t) = elem_prod( row(Nta_m,t), wtAge_m );
    Bta_f(t) = elem_prod( row(Nta_f,t), wtAge_f );
    SSBt(t)  = sum( elem_prod( row(Bta_f,t), matAge_f ) );
    Rt(t)    = Nta_f(t,1) + Nta_m(t,1);
    baranov_time = t;
    if( useBaranov ) // get Z for this year
    {
      solveBaranov_multi_fleet();
    }
    //---------------------
    // End biomass dynamics
    //---------------------
    //-----------------------------------------------------
    // Observation models
    //-----------------------------------------------------
    calc_observation_models();
  if( verbose == 1 )
    cout << "Done with observation models for time " << t << endl;
  } // end loop over year
  //exit(1);
  if( sd_phase() ) 
  {
    D2015 = SSBt(nT)/SSB0;
    B2015 = SSBt(nT);
  }
  if ( last_phase() ) neval+=1; 
  else neval=-2.;  // start neval counter at -2 at start of last phase so it equals admb screen output
}

void model_parameters::calc_observation_models(void)
{
  int t = baranov_time;
  int g;
  //-----------------------------------------------------
  // Observation models for indices
  //-----------------------------------------------------
  // predicted exploitable biomass by fishery (g=1)
  // and survey (g=2,3)
  for( int g=1; g<=nFisheries; g++ )
  {
    expBit(g,t)  = sum( elem_prod( elem_prod( row(Bta_m,t),row(sel_gta_m(g),t) ), 
                                              mfexp(-fracYearFleet(g)*Zta_m(t)) ) );      
    expBit(g,t) += sum( elem_prod( elem_prod( row(Bta_f,t),row(sel_gta_f(g),t) ), 
                                              mfexp(-fracYearFleet(g)*Zta_f(t)) ) );      
  }
  //-----------------------------------------------------
  // Observation models for age-composition
  //-----------------------------------------------------          
  dvar_vector tmp(1,plusGroupAge);
  dvariable sumP;
  dvariable fracYrSurvey;
  for( int i=1; i<=nAgeSeries; i++ )
  {
    if( verbose == 1 )
      cout << "Inside calc_observation_models for fleet " << i << endl;
    int g = fisheryAgeID(i);
    if( verbose)
      cout << "fisheryAgeID(" << i << ") = " << g << endl;
    //cout <<"t = "<< t <<" g = "<< g << endl;
    //cout <<  "sum_m = " << sum( row(ageObsProp_m(i),t) ) << endl;
    //cout <<  "sum_f = " << sum( row(ageObsProp_f(i),t) ) << endl;
    // Male age comps
    if( sum( row(ageObsProp_m(i),t) ) > 0. )
    {
      tmp.initialize();
      sumP=0.;
      // calc raw proportions-at-age
      for( int a=minAge(i); a<=plusGroupAge; a++ )
      {
        tmp(a) = Nta_m(t,a)*sel_gta_m(g,t,a)*exp(-fracYearAges(i)*Zta_m(t,a));
        if( (g < nFisheries-2) & last_phase() )
          tmp(a) *= reportRate_ga_m(g,a);
        sumP  += tmp(a);
      }
      // normalize proportions
      for( int a=minAge(i); a<=plusGroupAge; a++ )
      {
        uCgta_m(i,t,a) = tmp(a)/sumP;
      }
      if( useAgeError_m(i) )
        predObsProp_m(i,t)(minAge(i),maxAge(i)) = Q*row(uCgta_m(i),t);
      else
        predObsProp_m(i,t)(minAge(i),maxAge(i)) = row(uCgta_m(i),t);
      //cout<< "predObsProp_m = " <<  predObsProp_m(i,t) << endl;
      //cout<< "reportRate_ga_m = " <<  reportRate_ga_m(g) << endl;
      //if( g <= nFisheries-2 )
        //predObsProp_m(i,t) = predObsProp_m(i,t) * reportRate_ga_m(g);
    } // close if(obs>0)
    if( verbose == 1 )
      cout << "Done with male age obs" << endl;
    // Female age comps
    if( sum( row(ageObsProp_f(i),t) ) > 0. )
    {
      tmp.initialize();
      sumP=0.;
      for( int a=minAge(i); a<=plusGroupAge; a++ )
      {
        tmp(a) = Nta_f(t,a)*sel_gta_f(g,t,a)*exp(-fracYearAges(i)*Zta_f(t,a));
        if( (g < nFisheries-2) &  last_phase() )
          tmp(a) *= reportRate_ga_f(g,a);
        sumP  += tmp(a);
      }
      for( int a=minAge(i); a<=plusGroupAge; a++ )
      {
        uCgta_f(i,t,a) = tmp(a)/sumP;
      }
      if( useAgeError_f(i) )
        predObsProp_f(i,t)(minAge(i),maxAge(i)) = Q*row(uCgta_f(i),t);
      else
        predObsProp_f(i,t)(minAge(i),maxAge(i)) = row(uCgta_f(i),t);
      //cout<< "predObsProp_f = " <<  predObsProp_f(i,t) << endl;
      //cout<< "reportRate_ga_f = " <<  reportRate_ga_f(g) << endl;
      //if( g <= nFisheries-2 )
        //predObsProp_f(i,t) = predObsProp_f(i,t) * reportRate_ga_f(g);
    } // close if(obs>0)    
    if( verbose == 1 )
      cout << "Done with female age obs at time " << t << endl;
  }
  //-----------------------------------------------------
  // Observation models for sub-legal releases by gear
  //-----------------------------------------------------
  if( verbose == 1 )
      cout << "Calculating sublegal releases for year " << t << endl;
  dvar_vector tmpDirCatch_m(1,plusGroupAge);
  dvar_vector tmpDirCatch_f(1,plusGroupAge);
  dvar_vector tmpUndirCatch_m(1,plusGroupAge);
  dvar_vector tmpUndirCatch_f(1,plusGroupAge);
  dvar_vector tmpDirRel(1,plusGroupAge);
  dvar_vector tmpDirLand(1,plusGroupAge);
  dvar_vector tmpUndirRel(1,plusGroupAge);
  dvar_vector tmpUndirLand(1,plusGroupAge);
  for( int g=1;g<=nFisheries-2;g++ )
  {
    // Initialize released and landed catch
    relCtg(t,g) = 0.; landedCtg(t,g) = 0.; 
    undirLandCtg(t,g) = 0.;
    dirLandCtg(t,g) = 0.;
    undirRelCtg(t,g) = 0.;
    dirRelCtg(t,g) = 0.;
    // accumulate legal catch, discards, and sublegal total discards
    tmpDirCatch_m.initialize(); tmpDirCatch_f.initialize();
    tmpUndirCatch_m.initialize(); tmpUndirCatch_f.initialize();
    tmpDirRel.initialize(); tmpDirLand.initialize();
    tmpUndirRel.initialize(); tmpUndirLand.initialize();
    dvar_vector pR_m = row(pRet_m,t);
    dvar_vector pR_f = row(pRet_f,t);
    for( int a=1; a<=plusGroupAge; a++ )
    {
        // get onboard catch-at-age
        tmpDirCatch_m(a) = pDir_tg(t,g)*Bta_m(t,a)*sel_gta_m(g,t,a)*(Ftg(t,g))*(1.-exp(-Zta_m(t,a)));
        tmpDirCatch_m(a) /= Zta_m(t,a);
        tmpDirCatch_f(a) = pDir_tg(t,g)*Bta_f(t,a)*sel_gta_f(g,t,a)*Ftg(t,g)*(1.-exp(-Zta_f(t,a)));
        tmpDirCatch_f(a) /= Zta_f(t,a);
        // determine how much onboard catch is released at-sea
        // Directed catch
        tmpDirRel(a) = tmpDirCatch_m(a)*(1.-pR_m(a)) + tmpDirCatch_f(a)*(1.-pR_f(a));
        tmpDirLand(a) = value(tmpDirCatch_m(a)*pR_m(a) + tmpDirCatch_f(a)*pR_f(a) );
        // Undirected catch  
        //   Separate fleet F into a directed and undirected portion (as below)
        //   and deviate peak selectivity downwards in years where we have data
        //
        tmpUndirCatch_m(a) = (1-pDir_tg(t,g)) * Bta_m(t,a)*selUndir_gta_m(g,t,a)*Ftg(t,g)*(1.-exp(-Zta_m(t,a)));
        tmpUndirCatch_m(a) /= Zta_m(t,a);
        tmpUndirCatch_f(a) = (1-pDir_tg(t,g)) * Bta_f(t,a)*selUndir_gta_f(g,t,a)*Ftg(t,g)*(1.-exp(-Zta_f(t,a)));
        tmpUndirCatch_f(a) /= Zta_f(t,a);
        tmpUndirRel(a) = tmpUndirCatch_m(a)*(1.-pR_m(a)) + tmpUndirCatch_f(a)*(1.-pR_f(a));
        tmpUndirLand(a) = value(tmpUndirCatch_m(a)*pR_m(a) + tmpUndirCatch_f(a)*pR_f(a) );   
    }        
    landedCtg(t,g) = sum(tmpDirLand + tmpUndirLand);
    relCtg(t,g) = sum(tmpDirRel + tmpUndirRel);
  }
  //-----------------------------------------------------
  // Observation models for landed catch by year/gear
  //-----------------------------------------------------
  // males                               females
  dvar_vector Bprime_m(1,plusGroupAge);  dvar_vector Bprime_f(1,plusGroupAge);
  dvar_vector Bprime_mU(1,plusGroupAge); dvar_vector Bprime_fU(1,plusGroupAge);
  dvar_vector pR_m = row(pRet_m,t);      dvar_vector pR_f = row(pRet_f,t);
  dvar_vector Za_m = row(Zta_m,t);       dvar_vector Za_f = row(Zta_f,t);
  dvar_vector tmp_m(1,plusGroupAge);     dvar_vector tmp_f(1,plusGroupAge);
  dvar_vector tmp_mU(1,plusGroupAge);    dvar_vector tmp_fU(1,plusGroupAge);
  // Predicted combined male/female catch by fishery
  for( int g=1;g<=nFisheries-2;g++ )
  {
    // compute predicted total catch
    // Fully-selected total biomass of males and females
    Bprime_m = elem_prod( row(sel_gta_m(g),t),row(Bta_m,t) );
    Bprime_f = elem_prod( row(sel_gta_f(g),t),row(Bta_f,t) );
    Bprime_mU = elem_prod( row(selUndir_gta_m(g),t),row(Bta_m,t) );
    Bprime_fU = elem_prod( row(selUndir_gta_f(g),t),row(Bta_f,t) );
    tmp_m = elem_div( elem_prod( Bprime_m*Ftg(t,g),1.-mfexp(-Za_m) ), Za_m );
    tmp_m = pDir_tg(t,g) * elem_prod( tmp_m, row(pR_gta_m(g),t) );
    tmp_f = elem_div( elem_prod( Bprime_f*Ftg(t,g),1.-mfexp(-Za_f) ), Za_f );
    tmp_f = pDir_tg(t,g) * elem_prod( tmp_f, row(pR_gta_f(g),t) );
    tmp_mU = elem_div( elem_prod( Bprime_mU*Ftg(t,g),1.-mfexp(-Za_m) ), Za_m );
    tmp_mU = (1 - pDir_tg(t,g) ) * elem_prod( tmp_m, row(pR_gta_m(g),t) );
    tmp_fU = elem_div( elem_prod( Bprime_fU*Ftg(t,g),1.-mfexp(-Za_f) ), Za_f );
    tmp_fU = (1 - pDir_tg(t,g) ) * elem_prod( tmp_f, row(pR_gta_f(g),t) );
    pCtg(t,g) = sum( tmp_m + tmp_f + tmp_mU + tmp_fU );
    //if( g==2 )
      //cout << pCtg(t,2) << endl;
  }
  //cout << row(pCtg,t) << endl;
  //if( t==51 ) exit(1);
  // At end, get legal and sublegal biomasses
  if( last_phase() )
  {
    legalBt(t) = 0.; subLegalBt(t) = 0.;
    for( int a=1; a<=plusGroupAge; a++ )
    {
      //legalBt(t)     += value( Bta_m(t,a)*pRet_m(t,a) + Bta_f(t,a)*pRet_f(t,a) );
      //subLegalBt(t)  += value( Bta_m(t,a)*(1.-pRet_m(t,a)) + Bta_f(t,a)*(1.-pRet_f(t,a)) );
      legalBt(t)     += Bta_m(t,a)*pRet_m(t,a) + Bta_f(t,a)*pRet_f(t,a);
      subLegalBt(t)  += Bta_m(t,a)*(1.-pRet_m(t,a)) + Bta_f(t,a)*(1.-pRet_f(t,a));
      // save sdreport_numbers
      legalHR    = elem_div( rowsum(landedCtg), legalBt(1,nT) );
      LHR2015    = legalHR(nT);
      subLegalHR = elem_div( rowsum(relCtg), subLegalBt(1,nT) );
      SLHR2015   = subLegalHR(nT);        
    }
  }
}

void model_parameters::solveBaranov_multi_fleet(void)
{
  RETURN_ARRAYS_INCREMENT();
  int t = baranov_time;
  dvariable f;
  dvariable J;
  dvar_vector Bprime_m(1,plusGroupAge);
  dvar_vector Bprime_f(1,plusGroupAge);
  dvar_vector tmp_m(1,plusGroupAge);
  dvar_vector tmp_f(1,plusGroupAge);
  dvar_vector Za_m(1,plusGroupAge); 
  dvar_vector Za_f(1,plusGroupAge); 
  dvar_vector ZaNew_m(1,plusGroupAge);  
  dvar_vector ZaNew_f(1,plusGroupAge);  
  dvar_vector pR_m = row(pRet_m,t);
  dvar_vector pR_f = row(pRet_f,t);
  // Initialize Z to current vector of M
  ZaNew_m.initialize(); ZaNew_f.initialize(); 
  ZaNew_m = M_m;
  ZaNew_f = M_f; 
  // Initial approximation of F and Z...
  for( int g=1;g<=nFisheries-2;g++ )
  {  
    // Fully-selected total biomass of legal males and females
    Bprime_m = elem_prod(elem_prod( row(sel_gta_m(g),t),row(Bta_m,t) ), row(pR_gta_m(g),t) );
    Bprime_f = elem_prod(elem_prod( row(sel_gta_f(g),t),row(Bta_f,t) ), row(pR_gta_f(g),t) );
    /*cout << "g=" << g << endl;
    cout << Bprime_m << endl;
    cout << Bprime_f << endl;
    */
    Ftg(t,g) = Ctg(t,g)/(sum(Bprime_m) + sum(Bprime_f) );
    ZaNew_m += elem_prod(sel_gta_m(g,t)*Ftg(t,g),row(pR_gta_m(g),t)+(1.-row(pR_gta_m(g),t))*dM(g));
    ZaNew_f += elem_prod(sel_gta_f(g,t)*Ftg(t,g),row(pR_gta_f(g),t)+(1.-row(pR_gta_f(g),t))*dM(g));
    /*cout <<"t = "<< t <<" g = "<< g << endl;
    cout << " sum(Bprime_m) = "<< sum(Bprime_m) <<" sum(Bprime_f) = "<< sum(Bprime_f) << endl;   
    cout << " Ftg = "<< Ftg(t, g) <<" Ctg = "<< Ctg(t, g) << endl;   
    //cout << " ZaNew_m = " << ZaNew_m << endl;
    //cout << " ZaNew_f = " << ZaNew_f << endl;  
    cout << " pR_m = " << pR_m << endl;
    cout << " pR_f = " << pR_f << endl;
    */      
  }
  // refine F for major fisheries (surveys not critical)
  for( int i=1; i<=baranovIter; i++ )
  {
    // Save old Z value and then revise new
    Za_m = ZaNew_m; ZaNew_m = M_m;
    Za_f = ZaNew_f; ZaNew_f = M_f;
    for( int g=1;g<=nFisheries-2;g++ )
    {
      f = 0.;
      // compute predicted total catch
      // Fully-selected total biomass of legal males and females
      Bprime_m = elem_prod( row(sel_gta_m(g),t),row(Bta_m,t) );
      Bprime_f = elem_prod( row(sel_gta_f(g),t),row(Bta_f,t) );
      /*
      cout << "g=" << g << endl;
      cout << Bprime_m << endl;
      cout << Bprime_f << endl;
      */
      tmp_m = elem_div( elem_prod( Bprime_m*Ftg(t,g),1.-mfexp(-Za_m) ), Za_m );
      tmp_m = elem_prod( tmp_m, row(pR_gta_m(g),t) );
      tmp_f = elem_div( elem_prod( Bprime_f*Ftg(t,g),1.-mfexp(-Za_f) ), Za_f );
      tmp_f = elem_prod( tmp_f, row(pR_gta_f(g),t) );
      // Function value: difference of pred - obs catch
      f =  sum(tmp_m + tmp_f) - Ctg(t,g);
      pCtg(t,g) = sum( tmp_m + tmp_f );
      //cout <<"iter = "<< i << " Ftg = "<< Ftg(t, g) << " pCtg = " << pCtg(t,g) << endl;   
      tmp_m=0.;tmp_f=0.;
      // Jacobian: males
      dvar_vector tmp1 = elem_prod(row(pR_gta_m(g),t),Bprime_m);
      dvar_vector tmp2 = elem_prod(sel_gta_m(g,t),row(pR_gta_m(g),t)+(1.-row(pR_gta_m(g),t))*dM(g));
      dvar_vector tmp3 = 1.-mfexp(-Za_m);
      dvar_vector tmp4 = elem_prod(tmp1,tmp3);
      dvar_vector tmp5 = elem_div(elem_prod(elem_prod(Ftg(t,g)*tmp1,mfexp(-Za_m)),tmp2),Za_m);
      dvar_vector tmp6 = elem_div(elem_prod(elem_prod(Ftg(t,g)*tmp1,tmp3),tmp2),elem_prod(Za_m,Za_m));
      dvar_vector dfm  = elem_div(tmp4,Za_m) + tmp5 - tmp6;
      tmp1 = 0.; tmp2 = 0.; tmp3 = 0.; tmp4 = 0.; tmp5 = 0.; tmp6 = 0.;
      // Jacobian: females
      tmp1 = elem_prod(row(pR_gta_f(g),t),Bprime_f);
      tmp2 = elem_prod(sel_gta_f(g,t),row(pR_gta_f(g),t)+(1.-row(pR_gta_f(g),t))*dM(g));
      tmp3 = 1.-mfexp(-Za_f);
      tmp4 = elem_prod(tmp1,tmp3);
      tmp5 = elem_div(elem_prod(elem_prod(Ftg(t,g)*tmp1,mfexp(-Za_f)),tmp2),Za_f);
      tmp6 = elem_div(elem_prod(elem_prod(Ftg(t,g)*tmp1,tmp3),tmp2),elem_prod(Za_f,Za_f));
      dvar_vector dff  = elem_div(tmp4,Za_f) + tmp5 - tmp6;
      J = sum(dfm+dff);
      Ftg(t,g) -= f/J;
      tmp1 = 0.; tmp2 = 0.; tmp3 = 0.; tmp4 = 0.; tmp5 = 0.; tmp6 = 0.;
      ZaNew_m += elem_prod( row(sel_gta_m(g),t)*Ftg(t,g),row(pR_gta_m(g),t)+(1.-row(pR_gta_m(g),t))*dM(g));
      ZaNew_f += elem_prod( row(sel_gta_f(g),t)*Ftg(t,g),row(pR_gta_f(g),t)+(1.-row(pR_gta_f(g),t))*dM(g));
      /*cout <<"iter = "<< i <<" g = "<<g << " t = " << t << " phase = "<< current_phase() <<endl;
      cout << "pR_m = " << pR_m << endl;
      cout << "pR_f = " << pR_f << endl;
      cout << "Ftg = "<< Ftg(t, 1) << endl;   
      cout <<"f = "<< f << " J = "<< J << " dfm = " << sum(dfm) << " dff = " << sum(dff) << endl;   
      */
    }
  }
  //exit(1);
  for( int g=1;g<=nFisheries-2;g++ )
  {
    f = 0.;
    // compute predicted total catch
    // Fully-selected total biomass of males and females
    Bprime_m = elem_prod( row(sel_gta_m(g),t),row(Bta_m,t) );
    Bprime_f = elem_prod( row(sel_gta_f(g),t),row(Bta_f,t) );
    tmp_m = elem_div( elem_prod( Bprime_m*Ftg(t,g),1.-mfexp(-Za_m) ), Za_m );
    tmp_m = elem_prod( tmp_m, row(pR_gta_m(g),t) );
    tmp_f = elem_div( elem_prod( Bprime_f*Ftg(t,g),1.-mfexp(-Za_f) ), Za_f );
    tmp_f = elem_prod( tmp_f, row(pR_gta_f(g),t) );
    pCtg(t,g) = sum( tmp_m + tmp_f );    
    tmp_m=0.;tmp_f=0.;
  }
  Zta_m(t) = ZaNew_m;
  Zta_f(t) = ZaNew_f;
  RETURN_ARRAYS_DECREMENT();
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  time(&finish);
  elapsed_time=difftime(finish,start);
  hour=long(elapsed_time)/3600;
  minute=long(elapsed_time)%3600/60;
  second=(long(elapsed_time)%3600)%60;
  report<<endl<<endl<<"## *******************************************"<<endl;
  report<<"## --Start time: "<<ctime(&start)<<endl;
  report<<"## --Finish time: "<<ctime(&finish)<<endl;
  report<<"## --Runtime: ";
  report<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
  report<<"## *******************************************"<<endl;
  // Minimization performance details
  report << endl;
  report << "## Minimization performance" << endl;
  report << "# validObs_idx" << endl;           report << validObs_idx << endl;
  report << "# n_residuals_f" << endl;          report << n_residuals_f << endl;    
  report << "# n_residuals_m" << endl;          report << n_residuals_m << endl;    
  report << "# indexLikelihood" << endl;        report << indexLikelihood << endl;
  report << "# ageLikelihood_m" << endl;        report << ageLikelihood_m << endl;
  report << "# ageLikelihood_f" << endl;        report << ageLikelihood_f << endl;
  report << "# releaseLikelihood" << endl;      report << releaseLikelihood << endl;
  report << "# dataLikelihood" << endl;         report << dataLikelihood << endl;
  report << "# fLike" << endl;                 report << fLike << endl;
  report << "# catLike" << endl;                 report << catLike << endl;
  report << "# mPrior" << endl;                 report << mPrior << endl;
  report << "# ssbPrior" << endl;               report << ssbPrior << endl;
  report << "# selPrior" << endl;               report << selPrior << endl;
  report << "# recPrior"<< endl;                report << recPrior << endl;
  report << "# qdevPrior"<< endl;                report << qdevPrior << endl;
  report << "# hPrior"<< endl;                  report << hPrior << endl;
  report << "# pDirPrior"<< endl;                  report << pDirPrior << endl;
  report << "# dirCatPrior"<< endl;                  report << dirCatPrior << endl;
  report << "# undirCatPrior"<< endl;                  report << undirCatPrior << endl;
  report << "# objFun"  << endl;    report << *objective_function_value::pobjfun << endl;
  report << "# maxGrad" << endl;    report << objective_function_value::gmax << endl;
  report << "# exitCode"<< endl;    report << iexit << endl;
  report << "# funEvals"<< endl;    report << neval << endl;
  // Derived likelihood components
  report << "## Standard error estimates" << endl;
  report << "# tauIndex" << endl;   report << sqrt(tauSquareIndex) << endl;
  report << "# tauAge_m" << endl;   report << sqrt(tauSquareAges_m) << endl;
  report << "# tauAge_f" << endl;   report << sqrt(tauSquareAges_f) << endl;
  report << "# tauRel"   << endl;   report << sqrt(tauSquareReleases) << endl;
  // Parameter estimates
  report << "## Parameter estimates" << endl;
  report << "# avgR"   << endl;     report << avgR << endl;
  report << "# h"   << endl;        report << h << endl;
  report << "# SSB0"    << endl;    report << SSB0 << endl;
  report << "# D2015"    << endl;   report << D2015 << endl;
  report << "# R0"      << endl;    report << R0 << endl;
  report << "# M_m" << endl;        report << M_m << endl;
  report << "# M_f" << endl;        report << M_f << endl;
  report << "# sigma_R" << endl;    report << sigma_R << endl;
  report << "# rho" << endl;        report << rho << endl;
  report << "# recDevs" << endl;    report << omega_t << endl;
  report << "# log_Fg1"   << endl;      report << log_Fg1 << endl;
  report << "# log_Fdevs_gt"   << endl; report << log_Fdevs_gt << endl;
  // selectivity parameters in first year
  report << "# alpha_g1" << endl;    report << alpha_g1 << endl;
  report << "# beta_g1"  << endl;    report << beta_g1  << endl;
  // selectivity parameters for all years
  report << "# alpha_gt_m" << endl;    report << alpha_gt_m << endl;
  report << "# beta_gt_m"  << endl;    report << beta_gt_m  << endl;
  report << "# alpha_gt_f" << endl;    report << alpha_gt_f << endl;
  report << "# beta_gt_f"  << endl;    report << beta_gt_f  << endl;
  // Undirected selectivity parameters
  report << "# alphaU_gt_m" << endl;    report << alphaUndir_gt_m << endl;
  report << "# alphaU_gt_f" << endl;    report << alphaUndir_gt_f << endl;
  report << "# betaU_gt_m" << endl;     report << betaUndir_gt_m << endl;
  report << "# betaU_gt_f" << endl;     report << betaUndir_gt_f << endl;
  report << "# qg"       << endl;   report << mfexp( lnq ) << endl;
  for( int i=1; i <= nIndexSeries; i ++ )
  {
    int g=idxIndex[i];
    report << "# q_"<<g<<"t"<< endl; 
    report << exp(row(lnq_it,i))<< endl;
  }
  report << "# reportRate_gaf"  << endl;
  report << reportRate_ga_f << endl;
  report << "# reportRate_gam"  << endl;
  report << reportRate_ga_m << endl;
  // High-grading 50% length
  report << "# hg50" << endl;    report << mfexp(log_hg50) << endl;
  // Proportion directed
  report << "# pDir_t1" << endl;    report << pDir_t1 << endl;  
  report << "# pDir_t2" << endl;    report << pDir_t2 << endl;  
  report << "# pDir_tg" << endl;    report << pDir_tg << endl;  
  report << "# logit_pDir_g" << endl;    report << logit_pDir_g << endl;  
  report << "# ph_logit_pDir" << endl;    report << ph_logit_pDir << endl;  
  report << "# priorSD_logit_pDir" << endl;    report << priorSD_logit_pDir << endl;  
  report << "# priorSD_seldev_pDir" << endl;    report << priorSD_seldev_pDir << endl;  
  report << "# selAlphaDirDev_base " << endl;   report << selAlphaDirDev_base << endl;
  // Output concentrated obs error quantities
  report << "# validObs_idx" << endl; report << validObs_idx << endl;
  report << "# validObs_rel" << endl; report << validObs_rel << endl;
  // Derived variables
  report << endl;
  report << "## Derived variables" << endl;
  report << "# bh_a" << endl;      report << bh_a << endl;
  report << "# bh_b" << endl;      report << bh_b << endl;  
  report << "# projDep" << endl;    report << SSBt(nT+1)/SSBt(1) << endl;
  report << "# projSSB" << endl;    report << SSBt(nT+1) << endl;
  report << "# projLegalB" << endl; report << legalBt(nT+1) << endl;
  report << "# projSublegalB" << endl; report << subLegalBt(nT+1) << endl;  
  report << "# SSBt" << endl;       report << SSBt(1,nT) << endl;
  report << "# Rt" << endl;         report << Rt(1,nT) << endl;
  report << "# legalBt" << endl;    report << legalBt(1,nT) << endl;
  report << "# sublegalBt" << endl; report << subLegalBt(1,nT) << endl;
  // Legal and subLegal total harvest rates
  //legalHR    = elem_div( rowsum(landedCtg), legalBt(1,nT) );
  report << "# legalHR" << endl;    report << legalHR << endl;
  //subLegalHR = elem_div( value(rowsum(relCtg)), subLegalBt(1,nT) );
  report << "# subLegalHR" << endl; report << subLegalHR << endl;
  report <<"## Total Observed landings by fishery"<< endl;
  for( int g=1;g<=nFisheries-2;g++ )
  {
    report << "# Ctg"<<g<< endl; 
    report << column(Ctg,g)<< endl;
  }
  report <<"## Fishing mortality rates by fishery"<< endl;
  for( int g=1;g<=nFisheries-2;g++ )
  {
    report << "# Ftg"<<g<< endl; 
    report << column(Ftg,g)<< endl;
  }
  report <<"## Legal landings by fishery"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# landedCtg"<<g<< endl; 
    report << column(landedCtg,g)<< endl;
  }
  report <<"## Total model discards by fishery"<< endl;
  for( int g=1;g<=nFisheries-2;g++ )
  {
    report << "# relCtg"<<g<< endl; 
    report << column(relCtg,g)<< endl;
  }
  report <<"## Total observed discards by fishery"<< endl;
  for( int g=1;g<=nFisheries-2;g++ )
  {
    report << "# obs_relCtg"<<g<< endl; 
    report << column(obs_relCtg,g)<< endl;
  }
  report <<"## Undirected landings and releases by fishery"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# undirLandCtg"<<g<< endl; 
    report << column(undirLandCtg,g)<< endl;
  }
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# undirRelCtg"<<g<< endl; 
    report << column(undirRelCtg,g)<< endl;
  }
  report <<"## Directed landings and releases by fishery"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# dirLandCtg"<<g<< endl; 
    report << column(dirLandCtg,g)<< endl;
  }
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# dirRelCtg"<<g<< endl; 
    report << column(dirRelCtg,g)<< endl;
  }
  report <<"## Predicted total landings by fishery"<< endl;
  for( int g=1;g<=nFisheries-2;g++ )
  {
    report << "# pCtg" << g << endl; 
    report << column(pCtg,g)<< endl;
  }
  /*report <<"## Legal (dead) discards by fishery"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# legalDtg"<<g<< endl; 
    report << row(legalDtg,g)<< endl;
  }
  */
  report <<"## SubLegal (dead) discards by fishery"<< endl;
  dvector tmp_D(1,nT);
  for( int g=1;g<=nFisheries-2;g++ )
  {
    tmp_D = value(column(relCtg,g)*(1.-exp(-dM(g))));
    report << "# sublegalDtg"<<g<< endl; 
    report << tmp_D << endl;
  }
  report <<"## Exploitable biomass by fishery (vals < 0 are missing...)"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# expBit"<<g<< endl; 
    report << row(expBit,g) << endl;
  }
  report <<"## TotalMortality-at-age for males and females"<< endl;
  report << "# Zta_m" << endl;      report << Zta_m << endl;
  report << "# Zta_f" << endl;      report << Zta_f << endl;
  report <<"## Numbers-at-age for males and females"<< endl;  
  report << "# Nta_m" << endl;      report << Nta_m << endl;
  report << "# Bta_m" << endl;      report << Bta_m << endl;
  report << "# Nta_f" << endl;      report << Nta_f << endl;
  report << "# Bta_f" << endl;      report << Bta_f << endl;
  report << "# pRet_m" << endl;      report << pRet_m << endl;
  report << "# pRet_f" << endl;      report << pRet_f << endl;
  // Selectivity
  report << "# selType" << endl;    report << selType << endl;
  report <<"## Selectivity-at-age by fishery: males"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# sel_gta_m"<<g<< endl;
    report << sel_gta_m(g) << endl;      
  }  
  report <<"## Selectivity-at-age by fishery: females"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# sel_gta_f"<<g<< endl; 
    report << sel_gta_f(g)<< endl;
  }
  report <<"## Undirected Selectivity-at-age by fishery: males"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# selU_gta_m"<<g<< endl;
    report << selUndir_gta_m(g) << endl;      
  }  
  report <<"## Selectivity-at-age by fishery: females"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# selU_gta_f"<<g<< endl; 
    report << selUndir_gta_f(g)<< endl;
  }
  report << "## Selectivity blocks start and end timesteps" << endl;
  report << "# selBlockStart_g1" << endl;
  report << selBlockStart_g1 << endl;
  report << "# selBlockEnd_g1" << endl;
  report << selBlockEnd_g1 << endl;
  report << "# selBlockStart_g2" << endl;
  report << selBlockStart_g2 << endl;
  report << "# selBlockEnd_g2" << endl;
  report << selBlockEnd_g2 << endl;
  report << "# selBlockStart_g3" << endl;
  report << selBlockStart_g3 << endl;
  report << "# selBlockEnd_g3" << endl;
  report << selBlockEnd_g3 << endl;
  // Life History Schedules
  report << "# ages "    << endl;   report << ages << endl;
  report << "# lenAge_m" << endl;   report << lenAge_m << endl;
  report << "# lenAge_f" << endl;   report << lenAge_f << endl;
  report << "# wtAge_m"  << endl;   report << wtAge_m << endl;
  report << "# wtAge_f"  << endl;   report << wtAge_f << endl;
  report << "# matAge_f" << endl;   report << matAge_f << endl;
  report << "# sizeLim"<< endl;   report << sizeLim << endl;
  report << "# sizeLimYear"<< endl; report << sizeLimYear << endl;
  report <<"## Scaled CPUE indices by fishery (vals < 0 are missing...)"<< endl;
  for( int i=1;i<=nIndexSeries;i++ )
  {
    int g=idxIndex[i];
    report << "# Itg"<<g<< endl; 
    report << elem_div(row(idxSeries,i),exp(row(lnq_it,i)))<< endl;
  }
  report <<"## Predicted male age proportions by fishery"<< endl;
  for( int i=1;i<=nAgeSeries;i++ )
  {
    report << "# uCgta_m"<<fisheryAgeID(i)<< endl; 
    report << predObsProp_m(i) << endl;
  }
  report <<"## Predicted female age proportions by fishery"<< endl;
  for( int i=1;i<=nAgeSeries;i++ )
  {
    report << "# uCgta_f"<<fisheryAgeID[i]<< endl; 
    report << predObsProp_f[i] << endl;
  }
  //--------------------------------------------------------------------------//
  // Begin echo of inputs for guiMOd()   ARK (03-Sep-10)                      //
  //--------------------------------------------------------------------------//
  report << "## sableOpMod Report" << endl;
  // Echo operating model parameters.
  // Echo data inputs from opModCatch.dat.
  report << "# nT" << endl;
  report << nT << endl;
  report << "# nFisheries" << endl;
  report << nFisheries << endl;
  report << "# landCatchMatrix" << endl;
  report << landCatchMatrix << endl;
  report << "# releaseMatrix" << endl;
  report << releaseMatrix << endl;
  // Echo data inputs from opModIndex.dat.
  report << "# nIndexSeries" << endl;
  report << nIndexSeries << endl;
  report << "# idxIndex" << endl;
  report << idxIndex << endl;
  report << "# idxLikeWeight" << endl;
  report << idxLikeWeight << endl;
  report << "# idxFirstYear" << endl;
  report << idxFirstYear << endl;
  report << "# idxLastYear" << endl;
  report << idxLastYear << endl;
  report << "# idxSeries" << endl;
  report << idxSeries << endl;
  // Echo data inputs from opModAges.dat.
  report << "# ages" << endl;
  report << ages << endl;
  report << "# minAge" << endl;
  report << minAge << endl;
  report << "# plusGroupAge" << endl;
  report << plusGroupAge << endl;
  report << "# nAgeSeries" << endl;
  report << nAgeSeries << endl;
  report << "# fisheryAgeID" << endl;
  report << fisheryAgeID << endl;
  report << "## ageObsProp_m (written as nAgeSeries slices)" << endl;
  for ( int i=1; i<=nAgeSeries;i++ )
  {
    report << "# ageObsProp_m" << i << endl;
    report << ageObsProp_m[i] << endl;
  }
  report << "## ageObsProp_f (written as nAgeSeries slices)" << endl;
  for ( int i=1; i<=nAgeSeries;i++ )
  {
    report << "# ageObsProp_f" << i << endl;
    report << ageObsProp_f[i] << endl;
  }
  // Echo data from opModLifeScheds.dat.
  report << "# sizeLim" << endl;
  report << sizeLim << endl;
  report << "# dM" << endl;
  report << dM << endl;
  report << "# lenAge_m" << endl;
  report << lenAge_m << endl;
  report << "# lenAge_f" << endl;
  report << lenAge_f << endl;
  report << "# wtAge_m" << endl;
  report << wtAge_m << endl;
  report << "# wtAge_f" << endl;
  report << wtAge_f << endl;
  report << "# matAge_f" << endl;
  report << matAge_f << endl;
  // Echo fixed inputs from control file
  // Echo selectivity priors
  report << "p_alpha_gt" << endl;
  report << p_alpha_gt << endl;
  report << "p_beta_gt" << endl;
  report << p_beta_gt << endl;
  report << "sd_alpha_gt" << endl;
  report << sd_alpha_gt << endl;
  report << "sd_beta_gt" << endl;
  report << sd_beta_gt << endl;
  // Echo ageing error matrix
  report << "# Q" << endl;
  report << Q << endl;
  //--------------------------------------------------------------------------//
  // End echo of inputs for guiMod()                                          //
  //--------------------------------------------------------------------------//  
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_CMPDIF_BUFFER_SIZE(25000000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
