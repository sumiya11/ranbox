
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"

#include "TMath.h"


// -----------------------------------------------------------------------------

// The (approximation) upper bound on the number of signal events we expect to see
// if there is no signal.
// Nexp: expected events in the signal box
// Nobs: observed events in the signal box
// -----------------------------------------------------------------
// In our problem, we see Nexp events in the sideband, and set Nobs = Nexp 
double Sup_approx(double Nexp, double Nobs, double alpha=0.05) {
    return TMath::Max(
      0.,
      0.5*TMath::ChisquareQuantile(1 - alpha, 2*(Nobs + 1)) - Nexp
    );
}

// The "True" upper bound on the number of events,
// given that we see Nexp events in the sideband
//
// Approximately, Sup(Nexp) == Sup_approx(Nexp, Nexp)
double Sup(double Nexp, double alpha = 0.05, int iters=10)
{
  // S_up(N_exp) = SUM_{i=0}^{inf} { Pois(i|N_exp) * [ max(0, 0.5*F_chi2(1-alpha, 2(i+1))-N_exp] }
  double p = 0, ans = 0;
  // for (int i = 0; i < Nexp + 5 * sqrt(Nexp); i++)
  for (int i = 0; i < iters; i++)
  {
    p = 0.5 * TMath::ChisquareQuantile(1 - alpha, 2 * (i + 1)) - Nexp;
    ans += TMath::Poisson(i, Nexp) * TMath::Max(0., p);
  }
  return ans;
}

// -----------------------------------------------------------------------------

// The test statistic optimizing the upper limit on the number of events
// in the signal box, given
// sup - upper limit on the number of events and
// eff - efficiency of the signal box
double UpperLimitTS(double sup, double eff, double luminosity = 1)
{
  return sup / (eff * luminosity);
}


void Compare_2(int start = 1, int end = 500) {

    int n = end - start;
    double xs[n], sup[n], approx[n];

    for (int i = 0; i < n; i += 10) {
        xs[i]    = start + i;
        sup[i]   = Sup(start + i, 0.05, 2*(end + i)); 
        approx[i] = Sup_approx(start + i, start + i);
    
        cout << xs[i] << " " << sup[i] << " " << approx[i] << '\n';
    }

    

}