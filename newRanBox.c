
#include <vector>
#include <iostream>
#include <string>
#include <tuple>

using namespace std;

// Separate file
const string output_path = "RANBOX";
const string apath = "RANBOX";
const string ascii_path = "RANBOX";

/*
    Setup filenames for output based on the parameters used
*/
/*
std::tuple<string, string, string> setup_output()
{
    std::stringstream sstr;

    if (dataset == 0)
    {
        sstr << "/RB22_HEPMASS_" << id << "_Ntr" << Ntrials
             << "_NS" << Nsignal << "_NB" << Nbackground << "_A" << Algorithm;
    }

    if (useSB)
    {
        sstr << "_SB";
    }
    else
    {
        sstr << "_Vol";
    }

    string rootfile = outputPath + sstr.str() + ".root";
    string asciifile = asciiPath + sstr.str() + ".asc";
    string summaryfile = asciiPath + "/SummaryRB22_HEPMASS.asc";

    return std::make_tuple(rootfile, asciifile, summaryfile);
}

void log_header(ofstream &, ofstream &results, ofstream &summary)
{
    // Header Printout
    // ---------------
    cout << endl;
    cout << "  ------------------------------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << "                                           R  A  N  B  O  X    2 2" << endl;
    cout << "                                           -----------------------" << endl;
    cout << endl;

    // Set up ascii output
    // -------------------

    results << "  ----------------------------- " << endl;
    results << "           R A N B O X   2 2    " << endl;
    results << "  ----------------------------- " << endl
            << endl
            << endl;
    results << "  Parameters: " << endl;
    results << "  ----------------------------- " << endl
            << endl;
    results << "  Dataset       = ";
    results << "HEPMASS data" << endl;
    results << "  Nsignal      = " << Nsignal << endl;
    results << "  Nbackground  = " << Nbackground << endl;
    results << "  NAD          = " << NAD << endl;
    results << "  Nvar         = " << Nvar << endl;
    results << "  Algorithm    = " << Algorithm << endl;
    results << "  Speedup      = " << speedup << endl;
    results << "  useSB        = " << useSB << endl;
    results << "  useZPL       = " << useZPL << endl;
    results << "  Root file    = " << rootfile << endl;
    results << "  maxBoxVolume = " << maxBoxVolume << endl;
    results << "  maxGDLoops   = " << maxGDLoops << endl;
    results << "  Nremoved     = " << Nremoved << endl;
    results << "  id           = " << id << endl;
    results << endl;
}
*/

// Cholesky-Banachiewicz decomposition of covariance matrix, used to generate a multivariate
// Gaussian distribution in N dimensions.
// NB It returns false (without finishing the calculation of L) if A is not pos.def., thereby
// it can be directly used to check that A is indeed positive definite.
// ------------------------------------------------------------------------------------------

bool Cholesky(const vector<vector<double>>& A,
                vector<vector<double>>& L)
{
    // We take a covariance matrix A, of dimension NxN, and we find a lower triangular
    // matrix L such that LL^T = A. See https://en.wikipedia.org/wiki/Cholesky_decomposition
    // -------------------------------------------------------------------------------------
    double sum1;
    double sum2;
    int N = A.size();

    if (N != A[0].size()) return false;

    L = vector<vector<double>>(N, vector<double>(N, 0));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            sum1 = 0;
            sum2 = 0;
            if (j > 0)
            {
                for (int k = 0; k < j - 1; k++)
                {
                    sum1 += L[i][k] * L[j][k];
                }
                for (int k = 0; k < j - 1; k++)
                {
                    sum2 += L[j][k] * L[j][k];
                }
            }
            if (i == j)
            {
                L[j][j] = sqrt(A[j][j] - sum2);
            }
            else
            { // i>j
                L[i][j] = (A[i][j] - sum1) / L[j][j];
            } // all others remain zero
            if (L[i][j] != L[i][j])
                return false;
        }
    }
    return true;
}

void RanBox(int dataset = 0, int Ntrials = 1000, int Nsignal = 200, int Nbackground = 9800,
            int Algorithm = 5, bool useSB = true, int Nvar = 8)
{

    /*
        Here we select
            Algorithm = 0,
            PCA = false,
            useSB = true,
            Nremoved = 0,
            speedup = 1,
            NH0 = 1,
            useZPL=false
    */

    // we want to specialize on hepmass first
    bool mock = false;

    // Ivar_fix - ?
    std::vector<int> Ivar_fix = {17, 21, 9, 5, 10, 13, 14, 18, 26, 6, 24, 11};

    bool PCA = false;
    bool RemoveHighCorr = false;

    if (dataset == 0)
    {
        // NAD = 27;              // NAD dimensional feature space
        // NSEL = NAD - Nremoved; // == NAD
    }

    // NAD should be >= Nvar
    // Nvar should be >= 2

    // maxHalfMuRange -- ?
    double maxHalfMuRange = 0.35;
    // Gaussian_dims -- ?
    int Gaussian_dims = 15;

    bool plots = false;

    // FlatFrac: fraction of flat events in toy
    // (the rest are multivariate normal)
    double FlatFrac = Nbackground / (Nbackground + Nsignal);

    // setup for mock = false, dataset = 0
    /*
    auto filenames = setup_output();
    string rootfile = std::get<0>(filenames);
    string asciifile = std::get<1>(filenames);
    string summaryfile = std::get<2>(filenames);
    */
    string rootfile, asciifile, summaryfile; 

    // setup arrays for storing data (Create big array to store event kinematics)
    int ND = 50;           // total number of considered space dimensions
    int maxEvents = 10000; // hmmm..

    vector<float> feature_all(maxEvents * ND);
    vector<float> order_ind_all(maxEvents * ND);
    // These are esentially the same ?
    vector<vector<float>> feature(ND); // feature[i] points to feature_all
    vector<vector<int>> order_ind(ND); // order_ind[i] points to order_ind_all

    // double and floats ?
    vector<double> xtmp(maxEvents);
    vector<bool> isSignal(maxEvents);
    vector<int> cum_ind(maxEvents);

    // Arrays used for cluster search
    // ------------------------------
    vector<int> Closest(maxEvents);
    int maxClosest = 1; // ?????
    vector<int> PB_all(maxEvents * maxClosest);
    vector<vector<int>> PointedBy(maxClosest); // PointedBy[i] points to PB_all
    vector<int> AmClosest(maxEvents);

    // Logging header
    ofstream results; // output file
    ofstream summary; // summary file
    results.open(asciifile);
    summary.open(summaryfile, std::ios::out | std::ios::app);
    // log_header(results, summary);

    // Set some control ad-hoc variables

    // We omit plot constructing at this stage
    // Set plots and canvases

    // TString vs std::string
    // Not used, to delete
    vector<string> varname = {"V00", "V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08", "V09",
                           "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19",
                           "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29",
                           "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39",
                           "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49"};

    std::string progress = "[--10%--20%--30%--40%--50%--60%--70%--80%--90%-100%]";
    int currchar;
    int block;

    // here again, we have NH0 to be 1

    // not mock
    // read real data

    // fill feature[dim][goodevents] with lines of the dataset
    // isSignal[goodevents] = true for signal
    // Try to read Nsignal, Nbackground events and fill
    // NSread;
    // NBread;
    // goodevents(counts NSread and NBread);
    // read_data();

    // No PCA, data is pretransformed
    // (integral + compactify + remove high correlation + sort)

    for (int trial = 0; trial < Ntrials; trial++)
    {

        // case : NseedTrials == 1
        // set
        // Ivar[k] used[Ivar[k]]

        // Choose initial box

        // Box initializtion with Algorithm 0

        // Debug info
    }

} // end RanBox