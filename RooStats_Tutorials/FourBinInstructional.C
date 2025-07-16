#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"
 
using namespace RooFit;
using namespace RooStats;
 
void FourBinInstructional(bool doBayesian = false, bool doFeldmanCousins = false, bool doMCMC = false)
{
 
   // let's time this challenging example
   TStopwatch t;
   t.Start();
 
   // set RooFit random seed for reproducible results
   RooRandom::randomGenerator()->SetSeed(4357);
 
   // make model
   RooWorkspace *wspace = new RooWorkspace("wspace");
   wspace->factory("Poisson::on(non[0,1000], sum::splusb(s[40,0,100],b[100,0,300]))");
   wspace->factory("Poisson::off(noff[0,5000], prod::taub(b,tau[5,3,7],rho[1,0,2]))");
   wspace->factory("Poisson::onbar(nonbar[0,10000], bbar[1000,500,2000])");
   wspace->factory("Poisson::offbar(noffbar[0,1000000], prod::lambdaoffbar(bbar, tau))");
   wspace->factory("Gaussian::mcCons(rhonom[1.,0,2], rho, sigma[.2])");
   wspace->factory("PROD::model(on,off,onbar,offbar,mcCons)");
   wspace->defineSet("obs", "non,noff,nonbar,noffbar,rhonom");
 
   wspace->factory("Uniform::prior_poi({s})");
   wspace->factory("Uniform::prior_nuis({b,bbar,tau, rho})");
   wspace->factory("PROD::prior(prior_poi,prior_nuis)");
 
   // ----------------------------------
   // Control some interesting variations
   // define parameers of interest
   // for 1-d plots
   wspace->defineSet("poi", "s");
   wspace->defineSet("nuis", "b,tau,rho,bbar");
   // for 2-d plots to inspect correlations:
   //  wspace->defineSet("poi","s,rho");
 
   // test simpler cases where parameters are known.
   //  wspace->var("tau")->setConstant();
   //  wspace->var("rho")->setConstant();
   //  wspace->var("b")->setConstant();
   //  wspace->var("bbar")->setConstant();
 
   // inspect workspace
   //  wspace->Print();
 
   // ----------------------------------------------------------
   // Generate toy data
   // generate toy data assuming current value of the parameters
   // import into workspace.
   // add Verbose() to see how it's being generated
   std::unique_ptr<RooDataSet> data{wspace->pdf("model")->generate(*wspace->set("obs"), 1)};
   //  data->Print("v");
   wspace->import(*data);
 
   // ----------------------------------
   // Now the statistical tests
   // model config
   ModelConfig *modelConfig = new ModelConfig("FourBins");
   modelConfig->SetWorkspace(*wspace);
   modelConfig->SetPdf(*wspace->pdf("model"));
   modelConfig->SetPriorPdf(*wspace->pdf("prior"));
   modelConfig->SetParametersOfInterest(*wspace->set("poi"));
   modelConfig->SetNuisanceParameters(*wspace->set("nuis"));
   wspace->import(*modelConfig);
   // wspace->writeToFile("FourBin.root");
 
   // -------------------------------------------------
   // If you want to see the covariance matrix uncomment
   //  wspace->pdf("model")->fitTo(*data);
 
   // use ProfileLikelihood
   ProfileLikelihoodCalculator plc(*data, *modelConfig);
   plc.SetConfidenceLevel(0.95);
   LikelihoodInterval *plInt = plc.GetInterval();
   RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   plInt->LowerLimit(*wspace->var("s")); // get ugly print out of the way. Fix.
   RooMsgService::instance().setGlobalKillBelow(msglevel);
 
   // use FeldmaCousins (takes ~20 min)
   FeldmanCousins fc(*data, *modelConfig);
   fc.SetConfidenceLevel(0.95);
   // number counting: dataset always has 1 entry with N events observed
   fc.FluctuateNumDataEntries(false);
   fc.UseAdaptiveSampling(true);
   fc.SetNBins(40);
   PointSetInterval *fcInt = NULL;
   if (doFeldmanCousins) {                          // takes 7 minutes
      fcInt = (PointSetInterval *)fc.GetInterval(); // fix cast
   }
 
   // use BayesianCalculator (only 1-d parameter of interest, slow for this problem)
   BayesianCalculator bc(*data, *modelConfig);
   bc.SetConfidenceLevel(0.95);
   SimpleInterval *bInt = NULL;
   if (doBayesian && wspace->set("poi")->getSize() == 1) {
      bInt = bc.GetInterval();
   } else {
      cout << "Bayesian Calc. only supports on parameter of interest" << endl;
   }
 
   // use MCMCCalculator  (takes about 1 min)
   // Want an efficient proposal function, so derive it from covariance
   // matrix of fit
   std::unique_ptr<RooFitResult> fit{wspace->pdf("model")->fitTo(*data, Save())};
   ProposalHelper ph;
   ph.SetVariables((RooArgSet &)fit->floatParsFinal());
   ph.SetCovMatrix(fit->covarianceMatrix());
   ph.SetUpdateProposalParameters(true); // auto-create mean vars and add mappings
   ph.SetCacheSize(100);
   ProposalFunction *pf = ph.GetProposalFunction();
 
   MCMCCalculator mc(*data, *modelConfig);
   mc.SetConfidenceLevel(0.95);
   mc.SetProposalFunction(*pf);
   mc.SetNumBurnInSteps(500); // first N steps to be ignored as burn-in
   mc.SetNumIters(50000);
   mc.SetLeftSideTailFraction(0.5); // make a central interval
   MCMCInterval *mcInt = NULL;
   if (doMCMC)
      mcInt = mc.GetInterval();
 
   // ----------------------------------
   // Make some  plots
   TCanvas *c1 = (TCanvas *)gROOT->Get("c1");
   if (!c1)
      c1 = new TCanvas("c1");
 
   if (doBayesian && doMCMC) {
      c1->Divide(3);
      c1->cd(1);
   } else if (doBayesian || doMCMC) {
      c1->Divide(2);
      c1->cd(1);
   }
 
   LikelihoodIntervalPlot *lrplot = new LikelihoodIntervalPlot(plInt);
   lrplot->Draw();
 
   if (doBayesian && wspace->set("poi")->getSize() == 1) {
      c1->cd(2);
      // the plot takes a long time and print lots of error
      // using a scan it is better
      bc.SetScanOfPosterior(20);
      RooPlot *bplot = bc.GetPosteriorPlot();
      bplot->Draw();
   }
 
   if (doMCMC) {
      if (doBayesian && wspace->set("poi")->getSize() == 1)
         c1->cd(3);
      else
         c1->cd(2);
      MCMCIntervalPlot mcPlot(*mcInt);
      mcPlot.Draw();
   }
 
   // ----------------------------------
   // query intervals
   cout << "Profile Likelihood interval on s = [" << plInt->LowerLimit(*wspace->var("s")) << ", "
        << plInt->UpperLimit(*wspace->var("s")) << "]" << endl;
   // Profile Likelihood interval on s = [12.1902, 88.6871]
 
   if (doBayesian && wspace->set("poi")->getSize() == 1) {
      cout << "Bayesian interval on s = [" << bInt->LowerLimit() << ", " << bInt->UpperLimit() << "]" << endl;
   }
 
   if (doFeldmanCousins) {
      cout << "Feldman Cousins interval on s = [" << fcInt->LowerLimit(*wspace->var("s")) << ", "
           << fcInt->UpperLimit(*wspace->var("s")) << "]" << endl;
      // Feldman Cousins interval on s = [18.75 +/- 2.45, 83.75 +/- 2.45]
   }
 
   if (doMCMC) {
      cout << "MCMC interval on s = [" << mcInt->LowerLimit(*wspace->var("s")) << ", "
           << mcInt->UpperLimit(*wspace->var("s")) << "]" << endl;
      // MCMC interval on s = [15.7628, 84.7266]
   }
 
   t.Print();
}
