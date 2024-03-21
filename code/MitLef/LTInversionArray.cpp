#include <mex.h>
#include <matrix.h>
#include <limits> // Inf
#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <cmath>  // M_PI
#include <algorithm>
#include <vector>
#include <iostream>
#include <complex>

// alexander.pletzer@nesi.org.nz
// this is a translation of LTInversion function defined in ml.m to C++

// Constants
const double Inf = std::numeric_limits<double>::max();
const double pi = M_PI;
const double eps = 2.2204460492503131e-16; // MATLAB's eps
const std::complex<double> imag(0., 1.);


/* 
 * Comparison functor to sort an index array given the values of another array
 */
class LessThan {
public:
    
    /*
     Constructor
     @param array the "other" array 
     */
    LessThan(const std::vector<double>& array) {
        this->array = array;
    }
    
    /*
     * Comparison operator
     * @param ia
     * @param ib
     * @return true if array[ia] < array[ib], false otherwise
     */
    bool operator()(int ia, int ib) const { 
        return(this->array[ia] < this->array[ib]);
    }
private:
    std::vector<double> array; 
};

void EvalPoles(const std::complex<double>& lambda, double alpha, double beta, double gama, // input
              std::vector< std::complex<double> >& s_star, std::vector<double>& phi_s_star) { // output
    
    // Evaluation of the relevant poles
    // theta = angle(lambda) ;
    // kmin = ceil(-alpha/2 - theta/2/pi) ;
    // kmax = floor(alpha/2 - theta/2/pi) ;

    const double theta = std::arg(lambda);
    const int kmin = (int) std::ceil(-alpha/2 - theta/2/pi) ;
    const int kmax = (int) std::floor(alpha/2 - theta/2/pi) ;
    
    // k_vett = kmin : kmax ;
    // s_star = abs(lambda)^(1/alpha) * exp(1i*(theta+2*k_vett*pi)/alpha) ;
    std::vector<double> phi_s_star_old(kmax - kmin + 1);
    std::vector< std::complex<double> > s_star_old(kmax - kmin + 1);
    std::vector<int> index_s_star(phi_s_star_old.size());

    for (int k_vett = kmin; k_vett < kmax + 1; ++k_vett) {
        
        std::complex<double> ss = std::pow( std::abs(lambda), 1.0/alpha ) * std::exp(imag*(theta + 2.0*k_vett*pi)/alpha);
        s_star_old[k_vett - kmin] = ss;
        
        // Evaluation of phi(s_star) for each pole
        // phi_s_star = (real(s_star)+abs(s_star))/2 ;
        phi_s_star_old[k_vett - kmin] = ( ss.real() + std::abs(ss) )/2.0;
        
        index_s_star[k_vett - kmin] = k_vett - kmin; // 0-based indexing
    }
        
    // Sorting of the poles according to the value of phi(s_star)
    //[phi_s_star , index_s_star ] = sort(phi_s_star) ;
    LessThan lt(phi_s_star_old);
    std::sort(index_s_star.begin(), index_s_star.end(), lt);
    
    // Set s_star and phi_s_star
    s_star.resize(0);
    phi_s_star.resize(0);
    // Inserting the origin in the set of the singularities
    // s_star = [0, s_star] ;
    // phi_s_star = [0, phi_s_star] ;
    s_star.push_back( std::complex<double>(0., 0.) );
    phi_s_star.push_back(0.);
    
    for (auto i = 0; i < index_s_star.size(); ++i) {

        int index_new = index_s_star[i];
        
        double pss = phi_s_star_old[index_new];
        // Do not include zero poles
        if (pss > 1.e-15) {
            s_star.push_back(s_star_old[index_new]);
            phi_s_star.push_back(phi_s_star_old[index_new]);
        }
    }    
}

/*
 * Find the minimum value and the corresponding index
 * @param array array
 * @param min_element minimum (output)
 * @param min_index index of iminimum (ouput)
 */
void findMinValueAndIndex(const std::vector<std::size_t>& array, 
                         std::size_t& min_elem, std::size_t& min_index) {

    std::vector<std::size_t>::const_iterator it = std::min_element(array.begin(), array.end());
    min_elem = *it;
    min_index = std::distance(array.begin(), it) + 1; // MATLAB indices start at 1
}

void OptimalParam_RU(double t, double phi_s_star_j, double pj, double log_epsilon,
    double& muj, double& hj, std::size_t& Nj) {

    // Constants
    const double Inf = std::numeric_limits<double>::max();
    const double pi = M_PI;
    const double eps = 2.2204460492503131e-16; // MATLAB's eps

    // Evaluation of the starting values for sq_phi_star_j
    double sq_phi_s_star_j = sqrt(phi_s_star_j);
    double phibar_star_j;
    if (phi_s_star_j > 0) {
        phibar_star_j = phi_s_star_j*1.01;
    } else {
        phibar_star_j = 0.01;
    }
    double sq_phibar_star_j = sqrt(phibar_star_j) ;

    // Definition of some constants
    double f_min = 1, f_max = 10, f_tar = 5;

    // Iterative process to look for fbar in [f_min,f_max]
    bool stop = false;
    double phi_t, A, sq_muj, log_eps_phi_t, fbar;
    while (!stop) {
        phi_t = phibar_star_j*t;
        log_eps_phi_t = log_epsilon/phi_t ;
        Nj = (std::size_t) ceil(phi_t/pi*(1 - 3*log_eps_phi_t/2 + sqrt(1-2*log_eps_phi_t))) ;
        A = pi*Nj/phi_t;
        sq_muj = sq_phibar_star_j*std::abs(4.0 - A)/std::abs(7.0 - sqrt(1.0 + 12.0*A));
        fbar = pow((sq_phibar_star_j - sq_phi_s_star_j)/sq_muj, -pj);
        stop = (pj < 1.0e-14) || (f_min < fbar && fbar < f_max);
        if (!stop) {
            sq_phibar_star_j = pow(f_tar, -1/pj)*sq_muj + sq_phi_s_star_j;
            phibar_star_j = sq_phibar_star_j *  sq_phibar_star_j;
        }
    }

    muj = sq_muj * sq_muj;
    hj = (-3*A - 2 + 2*sqrt(1+12*A))/(4-A)/ (double) Nj;
    
    // Adjusting integration parameters to keep round-off errors under control
    double log_eps = log(eps);
    double Q, u, w;
    double threshold = (log_epsilon - log_eps)/t;
    if (muj > threshold) {
        if (std::abs(pj) < 1.0e-14) {
            Q = 0;
        } else {
            Q = pow(f_tar, -1/pj)*sqrt(muj);
        }
        phibar_star_j = pow(Q + sqrt(phi_s_star_j), 2);
        if (phibar_star_j < threshold) {
            w = sqrt(log_eps/(log_eps-log_epsilon)) ;
            u = sqrt(-phibar_star_j*t/log_eps) ;
            muj = threshold ;
            Nj = (std::size_t) ceil(w*log_epsilon/2/pi/(u*w-1)) ;
            hj = sqrt(log_eps/(log_eps - log_epsilon))/Nj ;
        }
        else {
            Nj = (std::size_t) Inf;
            hj = 0;
        }
    }
}

void OptimalParam_RB(double t, double phi_s_star_j, double phi_s_star_j1, double pj, double qj, double log_epsilon,
                    double& muj, double& hj, std::size_t& Nj) {

    // Constants
    const double Inf = std::numeric_limits<double>::max();
    const double pi = M_PI;
    // const double eps = 2.2204460492503131e-16; // MATLAB's eps

    // Definition of some constants
    double log_eps = -36.043653389117154; // log(eps)
    double fac = 1.01;
    bool conservative_error_analysis = false;

    // Maximum value of fbar as the ratio between tolerance and round-off unit
    double f_max = exp(log_epsilon - log_eps);

    // Evaluation of the starting values for sq_phi_star_j and sq_phi_star_j1
    double sq_phi_star_j = sqrt(phi_s_star_j);
    double threshold = 2*sqrt((log_epsilon - log_eps)/t);
    double sq_phi_star_j1 = std::min(sqrt(phi_s_star_j1), threshold - sq_phi_star_j);

    // Zero or negative values of pj and qj
    double sq_phibar_star_j, sq_phibar_star_j1; // NEED TO INITIALIZE?
    bool adm_region; // NEED TO INITIALIZE?
    if (pj < 1.0e-14 && qj < 1.0e-14) {
        sq_phibar_star_j = sq_phi_star_j;
        sq_phibar_star_j1 = sq_phi_star_j1;
        adm_region = true;
    }

    double f_min, fq, fp, f_bar;

    // Zero or negative values of just pj
    if (pj < 1.0e-14 && qj >= 1.0e-14) {
        sq_phibar_star_j = sq_phi_star_j;
        if (sq_phi_star_j > 0) {
            f_min = fac* pow(sq_phi_star_j/(sq_phi_star_j1-sq_phi_star_j), qj);
        }
        else {
            f_min = fac;
        }
        if (f_min < f_max) {
            f_bar = f_min + f_min/f_max*(f_max-f_min);
            fq = pow(f_bar, -1/qj);
            sq_phibar_star_j1 = (2*sq_phi_star_j1-fq*sq_phi_star_j)/(2+fq);
            adm_region = true;
        }
        else {
            adm_region = false;
        }
    }

    // Zero or negative values of just qj
    if (pj >= 1.0e-14 && qj < 1.0e-14) {
        sq_phibar_star_j1 = sq_phi_star_j1;
        f_min = fac*pow(sq_phi_star_j1/(sq_phi_star_j1-sq_phi_star_j), pj);
        if (f_min < f_max) {
            f_bar = f_min + f_min/f_max*(f_max-f_min);
            fp = pow(f_bar, -1/pj);
            sq_phibar_star_j = (2*sq_phi_star_j+fp*sq_phi_star_j1)/(2-fp) ;
            adm_region = 1 ;
        }
        else {
            adm_region = false;
        }
    }

    // Positive values of both pj and qj
    double w, den;
    if (pj >= 1.0e-14 && qj >= 1.0e-14) {
        f_min = fac*(sq_phi_star_j+sq_phi_star_j1) /
            pow(sq_phi_star_j1-sq_phi_star_j, std::max(pj,qj));
        if (f_min < f_max) {
            f_min = std::max(f_min,1.5);
            f_bar = f_min + f_min/f_max*(f_max-f_min);
            fp = pow(f_bar, -1/pj);
            fq = pow(f_bar, -1/qj);
            if (!conservative_error_analysis) {
                w = -phi_s_star_j1*t/log_epsilon;
            }
            else {
                w = -2*phi_s_star_j1*t/(log_epsilon-phi_s_star_j1*t);
            }
            den = 2+w - (1+w)*fp + fq;
            sq_phibar_star_j = ((2+w+fq)*sq_phi_star_j + fp*sq_phi_star_j1)/den;
            sq_phibar_star_j1 = (-(1+w)*fq*sq_phi_star_j 
                + (2+w-(1+w)*fp)*sq_phi_star_j1)/den;
            adm_region = true;
        }
        else {
            adm_region = false;
        }
    }

    if (adm_region) {
        log_epsilon = log_epsilon  - log(f_bar);
        if (!conservative_error_analysis) {
            w = -pow(sq_phibar_star_j1, 2)*t/log_epsilon;
        }
        else {
            w = -2*pow(sq_phibar_star_j1,2)*t/(log_epsilon- pow(sq_phibar_star_j1,2)*t) ;
        }
        muj = pow(((1+w)*sq_phibar_star_j + sq_phibar_star_j1)/(2+w), 2);
        hj = -2*pi/log_epsilon*(sq_phibar_star_j1-sq_phibar_star_j)
            /((1+w)*sq_phibar_star_j + sq_phibar_star_j1);
        Nj = (std::size_t) ceil(sqrt(1-log_epsilon/t/muj)/hj);
    }
    else {
        muj = 0;
        hj = 0;
        Nj = (std::size_t) Inf;
    }
}

void FindRegions(const std::vector<int>& admissible_regions, int J1, double t, // input
                 const std::vector<double>& phi_s_star, const std::vector<double>& p, // input
                 const std::vector<double>& q, double log_epsilon, // input
                 double& mu, double& h, std::size_t& N, std::size_t& iN) { // output 

    std::size_t n = admissible_regions.size();
    std::size_t JJ1 = admissible_regions[n - 1];

    std::vector<double> mu_vett(JJ1, Inf);
    std::vector<std::size_t> N_vett(JJ1, 0);
    std::vector<double> h_vett(JJ1, Inf);

    bool find_region = false;
    double muj, hj;
    std::size_t Nj;
    while (!find_region) {
        for (auto j1 : admissible_regions) {
            auto j1m1 = (std::size_t)(j1 - 1);
            if (j1 < J1) {
                OptimalParam_RB(
                    t, phi_s_star[j1m1], phi_s_star[(std::size_t) j1], p[j1m1], q[j1m1],log_epsilon, // input
                    muj, hj, Nj); // output
            }
            else {
                OptimalParam_RU(
                    t, phi_s_star[j1m1], p[j1m1],log_epsilon, // input
                    muj, hj, Nj) ; // output
            }
            mu_vett[j1m1] = muj;
            h_vett[j1m1] = hj;
            N_vett[j1m1] = Nj;
        }
                
        std::vector<std::size_t>::iterator N_vett_min_itr = std::min_element(N_vett.begin(), N_vett.end());
        
        if (*N_vett_min_itr > 200) {
            log_epsilon += log(10.0);
        }
        else {
            find_region = true;
        }
    }
    // Selection of the admissible region for integration which involves the
    // minimum number of nodes

    findMinValueAndIndex(
        N_vett, // input
        N, iN); // output    

    mu = mu_vett[iN - 1];
    h = h_vett[iN - 1];
}

void InverseLaplace(std::size_t N, double h, double mu, double t, // input
                    const std::complex<double>& lambda, // input
                    double alpha, double beta, double gama, // input
                    std::complex<double>& Integral) { // output
    
    Integral = std::complex<double>(0., 0.);
    int N_int = (int) N;
    for (int k = -N_int; k < N_int + 1; ++k) {
        double u = h*k;
        std::complex<double> z = mu * (imag*u + 1.0)*(imag*u + 1.0) ;
        std::complex<double> zd = 2*mu*(imag - u) ;
        std::complex<double> zexp = exp(z*t) ;
        std::complex<double> F = zd * std::pow(z, alpha*gama-beta) /  std::pow(std::pow(z, alpha) - lambda, gama) ;
        std::complex<double> S = zexp*F ;
        Integral += S ;
    }
    Integral *= h/2.0/pi/imag;
}

void  EvalResidues(std::size_t iN, double t, double alpha, double beta, // input
                   const std::vector< std::complex<double> >& s_star, // input
                   std::complex<double>& Residues) { // output
    
    const std::size_t n = s_star.size();
    
    // Evaluation of residues
    Residues = std::complex<double>(0., 0.);
    for (auto i = iN; i < n; ++i) {
        std::complex<double> ss = s_star[i];
        Residues += 1.0 / alpha * (std::pow(ss, 1-beta) * std::exp(t*ss));
    }
}


/* =========================================================================
 * Evaluation of the ML function by Laplace transform inversion
 * =========================================================================
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    // function E = LTInversionArray(z, alpha, beta, gama, log_epsilon)
    
    // Might want to add some checks here (nrhs, nlhs, types, etc)
    
    // Assume real z
    const double* z = mxGetPr(prhs[0]);
    const double alpha = mxGetScalar(prhs[1]);
    const double beta = mxGetScalar(prhs[2]);
    const double gama = mxGetScalar(prhs[3]);
    const double log_epsilon = mxGetScalar(prhs[4]);
    
    const double t = 1.0;

    std::size_t nelem = mxGetNumberOfElements(prhs[0]);
    mxArray* EMatrix = mxCreateDoubleMatrix(1, nelem, mxREAL);
    double* E = mxGetPr(EMatrix);

// spawn $OMP_NUM_THREADS threads. OMP_NUM_THREADS is set in SLURM when you set
// --cpus-per-task=N (i.e. export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}). When setting up
// parpool, make sure to set --cpus-per-task to num_parpool_works * $OMP_NUM_THREADS.
#pragma omp parallel for
    for (std::size_t ielem = 0; ielem < nelem; ++ielem) {
        
        if (std::abs(z[ielem]) < 1.e-15) {
            E[ielem] = std::tgamma(1.0/beta);
            continue;
        }

        const std::complex<double> lambda(z[ielem], 0.0);

        // Evaluation of the relevant poles
        std::vector< std::complex<double> > s_star; // 0-sized array
        std::vector<double> phi_s_star; // 0-sized array


        EvalPoles(lambda, alpha, beta, gama, // in
                  s_star, phi_s_star); // out

        //J1 = length(s_star) ; J = J1 - 1 ;
        std::size_t J1 = s_star.size();
        std::size_t J = J1 - 1;

        // Strength of the singularities
        // p = [ max(0,-2*(alpha*gama-beta+1)) , ones(1,J)*gama ]  ;
        // q = [ ones(1,J)*gama , +Inf] ;
        std::vector<double> p(J1);
        std::vector<double> q(J1);
        p[0] = std::max(0.0, -2.0*(alpha*gama - beta + 1.0));
        for (auto i = 0; i < J; ++i) {
            p[i + 1] = 1.0*gama;
            q[i] = 1.0*gama;
        }
        q[J] = +Inf;

        // phi_s_star = [phi_s_star, +Inf] ;
        phi_s_star.push_back(+Inf);

        // Looking for the admissible regions with respect to round-off errors
        // admissible_regions = find( ...
        // (phi_s_star(1:end-1) < (log_epsilon - log(eps))/t) & ...
        // (phi_s_star(1:end-1) < phi_s_star(2:end))) ;
        std::vector<int> admissible_regions;
        admissible_regions.reserve(phi_s_star.size());
        for (auto i = 0; i < J1; ++i) {
            if ( (phi_s_star[i] < (log_epsilon - log(eps))/t) && (phi_s_star[i] < phi_s_star[i+1]) ) {
                admissible_regions.push_back(i + 1); // 1-based indexing
            }
        }

        double mu, h;
        std::size_t N, iN;
        FindRegions(admissible_regions, J1, t, phi_s_star, p, q, log_epsilon,
                   mu, h, N, iN);

        std::complex<double> Integral;
        InverseLaplace(N, h, mu, t, lambda, alpha, beta, gama,
                       Integral);

        std::complex<double> Residues;
        EvalResidues(iN, t, alpha, beta, s_star,
                     Residues);

        // Assume real z (?)
        E[ielem] = (Integral + Residues).real(); 
    }

    // Load the output values
    plhs[0] = EMatrix;
}
