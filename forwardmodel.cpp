#include<bits/stdc++.h>
#include<vector>
#include<Eigen/SparseLU>
#include<Eigen/Eigenvalues>
#include<Eigen/SparseCore>
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<cmath>
#include<algorithm>
#include<fstream>
#define SIZE 960
using namespace std;
using namespace Eigen;
double rspd;
double omga;
void calculateGeneralizedEigenValueAndVector(Eigen::MatrixXd&, Eigen::MatrixXd&);
// Define a function to calculate vgk
double calculateVgk(const MatrixXd& x, double d, const MatrixXd& Ka1, const MatrixXd& Ka2, double omga, const MatrixXd& M, const MatrixXd& Cm, int Nnfo, int Nn) {
    int size = Nnfo + 2 * Nn;

    MatrixXd x_sub = x.block(0, 0, size, 1);
    MatrixXd x_transpose = x_sub.transpose();

    MatrixXd term1 = x_transpose * ((2 * d * Ka1) - Ka2) * VectorXd::LinSpaced(size, 1, size);

    MatrixXd M_sub = M.block(0, 0, size, size);
    MatrixXd term2 = 2 * omga * (x_transpose * M_sub * x_sub);

    MatrixXd term3 = x_transpose * Cm * x_sub;

    double vgk = term1.sum() / (term2(0, 0) - term3(0, 0));
    return vgk;
}

std::vector<std::complex<double>> findRoots(const std::vector<double>& coefficients) {
    // Degree of the polynomial
    int degree = 3;

    // Construct the companion matrix
    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();
    companion.block(1, 0, degree - 1, degree - 1).setIdentity();
    companion.col(degree - 1) = -Eigen::Map<const Eigen::VectorXd>(coefficients.data(), degree) / coefficients[degree];

    // Solve for the eigenvalues of the companion matrix
    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    Eigen::MatrixXcd eigenvalues = solver.eigenvalues();

    // Extract the roots
    std::vector<std::complex<double>> roots;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i).imag()) < 1e-4) {
            roots.push_back(eigenvalues(i).real());
        }
    }

    return roots;
}

tuple<double,double,double,vector<double>> raylee_lysmer(int Nn,vector<double> vsv,vector<double> vpv,vector<double> rhov,double f,vector<int> hv,double modn,int Nnf,vector<double> vpfv,vector<double> rhofv,vector<int> hfv){
    double kk=0.0, vpk=0.0, vgk=0.0;
    vector<double> ev;
    return make_tuple(kk,vpk,vgk,ev);

}
int main(){
    int Nn = 240;
    vector<int> hv(Nn,250);

    int Nnf = 0;
    vector<int> hfv(Nnf,100);

    int Nf = 56;
    double fmin = 0.1;
    double df = 0.01;
    vector<double> fks(Nf);
    for (int i = 0; i < Nf; i++) {
    fks[i] = fmin+(df*i);
    }
    Nf=109;
    vector<double> modnv(Nf);
    // Assign mode numbers to modnv
    for (int i = 0; i < Nf / 2; ++i) {
        modnv[i] = 1.0;
    }
    for (int i = Nf / 2; i < Nf; ++i) {
        modnv[i] = 2.0;
    }

    int size = fks.size();
    for(int i=3 ; i<size; i++){
        fks.push_back(fks[i]);
    }

    vector<double> vtypv(Nf, 0);

    //medium parameters
    double vpvsr = 1.7321;                 //Vp/Vs ratio
    double gardc = 309.6;                  // constant in Gardner relation
    double powr = 0.25;                    //exponent in Gardner relation

    int layrth1 = 5;   // thickness in elements, layrth1*h = thickness in meters
    int layrth2 = 10;  // thickness in elements, layrth2*h = thickness in meters
    int layrth3 = 50;

    double vplay1 = 4000; double vslay1 = vplay1 / vpvsr; double rholay1 = gardc * pow(vplay1, powr);
    double vplay2 = 3396; double vslay2 = vplay2 / vpvsr; double rholay2 = gardc * pow(vplay2, powr);
    double vplay3 = 4500; double vslay3 = vplay3 / vpvsr; double rholay3 = gardc * pow(vplay3, powr);
    double vplay4 = 6000; double vslay4 = vplay4 / vpvsr; double rholay4 = gardc * pow(vplay4, powr);

    vector<double> vpv;
    vector<double> vsv;
    vector<double> rhov;

    for(int i=0;i<Nn;i++){
        if(i<layrth1){vpv.push_back(vplay1);}
        else if(i< layrth1 + layrth2){vpv.push_back(vplay2);}
        else if(i< layrth1 + layrth2 + layrth3){vpv.push_back(vplay3);}
        else{vpv.push_back(vplay4);}

        if(i<layrth1){vsv.push_back(vslay1);}
        else if(i< layrth1 + layrth2){vsv.push_back(vslay2);}
        else if(i< layrth1 + layrth2 + layrth3){vsv.push_back(vslay3);}
        else{vsv.push_back(vslay4);}

        if(i<layrth1){rhov.push_back(rholay1);}
        else if(i< layrth1 + layrth2){rhov.push_back(rholay2);}
        else if(i< layrth1 + layrth2 + layrth3){rhov.push_back(rholay3);}
        else{rhov.push_back(rholay4);}
    }
    //make the fluid part of the model

    vplay4 = 1500;
    rholay4 = 1000;
    //the true model in the fluid
    vector<double> vpfv;
    vector<double> rhofv;

    int modn;
    vector<double> vout;
    int countr = 0;

    int st=0;

    Eigen::MatrixXd M (480,480);
    Eigen::MatrixXd L1(4,4);
    Eigen::MatrixXd L2(4,4);
    Eigen::MatrixXd L3(4,4);
    Eigen::MatrixXd M1(4,4);
    M1.setZero();

    L1.setZero();
    L2.setZero();
    L3.setZero();
    M.setZero();

    for(int i=0;i<1;i++){ // from here raylee lysmer starts
    //for(int i=0;i<fks.size();i++){  temporary only
        modn = modnv[countr++];
        double kk;
        double vpk;
        //double vgk;
        vector<double> kappafv;
        //the number of nodes in the fluid, based on the number of elements
        int Nnfo;
        if(Nnf>0){Nnfo=Nnf+1;}else{Nnfo=0;}
        MatrixXd Ka1 (Nnfo+(2*Nn),Nnfo+(2*Nn));
        MatrixXd Ka2 (Nnfo+(2*Nn),Nnfo+(2*Nn));
        MatrixXd Ka3 (Nnfo+(2*Nn),Nnfo+(2*Nn));
        Ka1.setZero();
        Ka2.setZero();
        Ka3.setZero();

        omga = 2*M_PI*fks[i];

        //make mu and lambda
        vector<double> muv(rhov.size());
        vector<double> lamdav(rhov.size());

        for(int i=0;i<rhov.size();i++){
            muv[i] = rhov[i]*vsv[i]*vsv[i];
            lamdav[i] = rhov[i]*vpv[i]*vpv[i]-2*muv[i];
        }

        // Eigen sparse matrices for M, Ka1, Ka2, Ka3

        for(int ii =0;ii<Nn;ii++)
        {

            int h = hv[ii];
            double mu = muv[ii];
            double lamda = lamdav[ii];

            //make elemental mass matrix

            M1(0,0)=h*rhov[ii]/2.0;
            M1(1,1)=h*rhov[ii]/2.0;
            M1(2,2)=h*rhov[ii]/2.0;
            M1(3,3)=h*rhov[ii]/2.0;




            // some alternate variables from Lysmer
            double alph = ((2*mu)+lamda)/6;
            double bet = mu/6;
            double theta = (mu+lamda)/4;
            double psi = (mu-lamda)/4;

            L1(0,0) = 2*alph*h;
            L3(0,0) = (6*bet/h);
            L2(0,1) = 2*psi;
            L1(0,2) = alph*h;
            L3(0,2) = -(6*bet/h);
            L2(0,3) = 2*theta;
            L2(1,0) = L2(0,1);
            L1(1,1) = 2*bet*h;
            L3(1,1) = (6*alph/h);
            L2(1,2) = -2*theta;
            L1(1,3) = bet*h;
            L3(1,3) = -(6*alph/h);
            L1(2,0) = L1(0,2);
            L3(2,0) = L3(0,2);
            L2(2,1) = L2(1,2);
            L1(2,2) = L1(0,0);
            L3(2,2) = L3(0,0);
            L2(2,3) = -2*psi;
            L2(3,0) = L2(0,3);
            L1(3,1) = L1(1,3);
            L3(3,1) = L3(1,3);
            L2(3,2) = L2(2,3);
            L1(3,3) = L1(1,1);
            L3(3,3) = L3(1,1);

        // Assemble mass and stiffness matrices from elemental matrices
        if (ii == (Nn - 1)) {
            M.block(st, st, 2, 2) += M1.block(0, 0, 2, 2);
            Ka1.block(st, st, 2, 2) += L1.block(0, 0, 2, 2);
            Ka2.block(st, st, 2, 2) += L2.block(0, 0, 2, 2);
            Ka3.block(st, st, 2, 2) += L3.block(0, 0, 2, 2);
        } else {
            M.block(st, st, 4, 4) += M1;
            Ka1.block(st, st, 4, 4) += L1;
            Ka2.block(st, st, 4, 4) += L2;
            Ka3.block(st, st, 4, 4) += L3;
        }
        st+=2;
    }


        MatrixXd Cm(Nnfo+(2*Nn),Nnfo+(2*Nn));
        Cm.setZero();


///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%
//% find the rayleigh/scholte wave speed which would exist if the solid model
//% were a halfspace with the minimum model velocity
//%
//% this is a lower bound on the velocity that can be passed to EIGS, based
//% on ARPACK
//%
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// Find the minimum value and its location
    auto min_element = std::min_element(vsv.begin(), vsv.end());
    double msval = *min_element;
    int msloc = std::distance(vsv.begin(), min_element);

    // Assign values
    double vsmay = msval;
    double vpmay = vpv[msloc];

// coefficients of rayleigh's polynomial

    double t1 = 1 / std::pow(vsmay, 6);
    double t2 = -8 / std::pow(vsmay, 4);
    double t3 = (24 / std::pow(vsmay, 2)) - (16 / std::pow(vpmay, 2));
    double t4 = -16 * (1 - std::pow(vsmay / vpmay, 2));

    // Coefficients of the polynomial
    std::vector<double> coefficients = {t4,t3,t2,t1};

    // Find the roots
    std::vector<std::complex<double>> roots = findRoots(coefficients);

    // Find the minimum element
    auto minElement = std::min_element(roots.begin(), roots.end(),
        [](const std::complex<double>& a, const std::complex<double>& b) {
            return std::abs(a) < std::abs(b);
        });

        rspd =  pow(real((*minElement)),0.5);

    MatrixXd lu(480,480);
    lu.setZero();

    MatrixXd ru = MatrixXd::Identity(480,480);


    MatrixXd ld1 (480,480);
    ld1<<(M.array() * omga * omga);

    MatrixXd ld2 (480,480);
    ld2<<Ka3;

    MatrixXd ld3 (480,480);
    ld3<<(Cm.array() * omga);

    MatrixXd l (480,960);
    l<< ld1 - ld2 - ld3 , Ka2;

    MatrixXd A (960,960);
    A<< lu,ru,l;

    MatrixXd nlu = MatrixXd::Identity(480,480);

    MatrixXd nru(480,480);
    nru.setZero();

    MatrixXd nld(480,480);
    nld.setZero();

    MatrixXd nrd(480,480);
    nrd<<Ka1;

    MatrixXd nu (480,960);
    nu<< nlu, nru;

    MatrixXd nd (480,960);
    nd<< nld, nrd;

    MatrixXd B (960,960);
    B<< nu, nd;

    // Get matrix size
    int size = A.cols();

//tuple<double,double> forward_values = calculateGeneralizedEigenValueAndVector(B,A);

    kk=0.0, vpk=0.0;// vgk=0.0;
    double ev[SIZE], ev_new[SIZE];
    double temp, lambda_new, lambda_old, error;
    int step = 1;
    int n = 960;
    error = 1;

    /* Reading Initial Guess Vector */
    for (i = 0; i < n; i++) {
        ev[i]=1;
    }

    /* Initializing Lambda_Old */
    lambda_old = omga/rspd;

    /* Multiplication */
    /* Setting label for repetition */
up:
    for (i = 0; i < n; i++) {
        temp = 0.0;
        for (int j = 0; j < n; j++) {
            temp = temp + A(i,j) * ev[j];
        }
        ev_new[i] = temp;
    }

    /* Replacing x by x_new */
    for (i = 0; i < n; i++) {
        ev[i] = ev_new[i];
    }

    /* Finding largest value from x */
    lambda_new = fabs(ev[1]);
    for (i = 1; i < n; i++) {
        if (fabs(ev[i]) > lambda_new) {
            lambda_new = fabs(ev[i]);
        }
    }

    /* Normalization */
    for (i = 0; i < n; i++) {
        ev[i] = ev[i] / B(i, i);
    }
    /* Checking Accuracy */
    if (fabs(lambda_new - lambda_old) > error) {
        lambda_old = lambda_new;
        step++;
        goto up;
    }

    /* Display */
    cout << endl << endl << "STEP-" << step << endl;
    cout << "Generalized Eigen Value = " << lambda_new << endl;
    cout << "Generalized Eigen Vector: " << endl;
    cout << "[";
    for (i = 1; i <= n; i++) {
        cout << ev[i] << "\t";
    }
    cout << "\b\b\b]"; /* \b is for backspace */
    // Print the result:
    //std::cout << "vgk = \n" << vgk << std::endl;
    kk = lambda_new;
    vpk = omga/kk;
    double d = kk;
    const int size_x = 480; // Size of the array x.
    double vgk = 0.0;
    for (int i = 0; i < size_x; ++i) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int j = 0; j < size_x; ++j) {
            sum1 += ev[j] * (2 * d * Ka1(i,j)) - Ka2(i,j) * (j + 1);
            sum2 += ev[j] * M(i,i) * ev[j];
        }
        vgk += sum1 / ((2 * omga * sum2) - (ev[i] * Cm(i,i)));
    }
    cout<<"\nkk: "<<kk<<endl;
    cout<<"\nvpk: "<<vpk<<endl;
    cout<<"\nvgk: "<<vgk<<endl;
}
vector<double> vflg(109,0);
vector<double> fks_orig(109,0);
vector<double> modnv_orig(109,0);
// Copy the elements from sourceVector to destinationVector
std::copy(fks.begin(), fks.end(), fks_orig.begin());
std::copy(modnv.begin(), modnv.end(), modnv_orig.begin());
for(int i=0;i<fks.size();i++){cout<<i+1<<" : "<<fks_orig[i]<<endl;}

    // if (vflg[0] == 1) {
    //     fks = {fks_orig[0] * 0.999, fks_orig[0], fks_orig[0] * 1.001};
    //     modnv = {modnv_orig[0], modnv_orig[0], modnv_orig[0]};
    // } else {
    //     fks.push_back(fks_orig[0]);
    //     modnv.push_back(modnv_orig[0]);
    // }

    // for (int ii = 1; ii < fks_orig.size(); ii++) {
    //     if (vflg[ii] == 1) {
    //         fks.push_back(fks_orig[ii] * 0.999);
    //         fks.push_back(fks_orig[ii]);
    //         fks.push_back(fks_orig[ii] * 1.001);
    //         modnv.push_back(modnv_orig[ii]);
    //         modnv.push_back(modnv_orig[ii]);
    //         modnv.push_back(modnv_orig[ii]);
    //     } else {
    //         fks.push_back(fks_orig[ii]);
    //         modnv.push_back(modnv_orig[ii]);
    //     }
    // }

    // // Print the values of fks and modnv
    // int z=0;
    // std::cout << "fks: ";
    // for (const auto& element : fks) {
    //     std::cout<<z++<<" : "  << element << endl;
    // }
    // std::cout << std::endl;

    // z=0;
    // std::cout << "modnv: ";
    // for (const auto& element : modnv) {
    //     std::cout<<z++<<" : "  << element << endl;
    // }
    // std::cout << std::endl;

//for(int i=0;i<fks.size();i++){cout<<i+1<<" : "<<fks_orig[i]<<endl;}
vector<double> vp (fks.size(),0);
vector<double> U (fks.size(),0);
vector<vector<double>> xm (2*Nn,vector<double>(fks.size()));
int length = 2 * Nn;
for (int countr = 0; countr < fks.size(); countr++) {
        double f = fks[countr];

        double kk, vpk, vgk;
        std::vector<double> x;

        tuple<double,double,double,vector<double>> result = raylee_lysmer(Nn, vsv, vpv, rhov, f, hv, modnv[countr + 1], Nnf, vpfv, rhofv, hfv);

        vp[countr] = get<1>(result);
        U[countr] = get<2>(result);
        for (int i = 0; i < length; i++) {
            xm[i][countr] = get<3>(result)[i];
        }
    }


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% construct phase sensitivity kernels
//%
//% these are derivative matrices and are extremely sparse, so the
//% necessary matrix-vector multiplications are hardcoded
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//% for all depths (except the bottom element) and frequencies
vector<vector<double>> snsmf (Nn,vector<double>(fks.size()));
vector<vector<double>> snsmflam (Nn,vector<double>(fks.size()));
vector<vector<double>> snsmfrho (Nn,vector<double>(fks.size()));
vector<vector<double>> snsmfh (Nn,vector<double>(fks.size()));

for(int ii=0;ii<Nn;ii++){
    double h = hv[ii];
    for (int countr = 0; countr < fks.size(); countr++) {
            double f = fks[countr];
            double kk, vpk, vgk;
            std::vector<double> x;

        //  % density sensitivity
            snsmfrho[ii][countr] = -(vp[countr] / (2 * U[countr])) *
                (xm[2 * ii - 1][countr] * (h / 2) * xm[2 * ii - 1][countr] +
                xm[2 * ii + 0][countr] * (h / 2) * xm[2 * ii + 0][countr] +
                xm[2 * ii + 1][countr] * (h / 2) * xm[2 * ii + 1][countr] +
                xm[2 * ii + 2][countr] * (h / 2) * xm[2 * ii + 2][countr]) *
                rhov[ii];

        //% lambda sensitivity
        //% k^2 term
            snsmflam[ii][countr] = (1 / (2 * vp[countr] * U[countr])) *
                (xm[2 * ii - 1][countr] * (h / 3) * xm[2 * ii - 1][countr] +
                 xm[2 * ii + 1][countr] * (h / 3) * xm[2 * ii + 1][countr] +
                 xm[2 * ii + 1][countr] * (h / 6) * xm[2 * ii - 1][countr] +
                 xm[2 * ii - 1][countr] * (h / 6) * xm[2 * ii + 1][countr]) *
                (rhov[ii] * (vpv[ii] * vpv[ii] - 2 * vsv[ii] * vsv[ii]));

        double pi = std::acos(-1.0); // Pi value
        //%k^1 term
            snsmflam[ii][countr] += (1 / (2 * (2 * pi * fks[countr]) * U[countr])) *
                (xm[2 * ii - 1][countr] * (1.0 / 2) * xm[2 * ii][countr] +
                 xm[2 * ii - 1][countr] * (-1.0 / 2) * xm[2 * ii + 2][countr] +
                 xm[2 * ii][countr] * (1.0 / 2) * xm[2 * ii - 1][countr] +
                 xm[2 * ii][countr] * (1.0 / 2) * xm[2 * ii + 1][countr] +
                 xm[2 * ii + 1][countr] * (1.0 / 2) * xm[2 * ii][countr] +
                 xm[2 * ii + 1][countr] * (-1.0 / 2) * xm[2 * ii + 2][countr] +
                 xm[2 * ii + 2][countr] * (-1.0 / 2) * xm[2 * ii - 1][countr] +
                 xm[2 * ii + 2][countr] * (-1.0 / 2) * xm[2 * ii + 1][countr]) *
                (rhov[ii] * (vpv[ii] * vpv[ii] - 2 * vsv[ii] * vsv[ii]));

        //% k^0 term
                snsmflam[ii][countr] += (vp[countr] / (2 * std::pow(2 * pi * fks[countr], 2) * U[countr])) *
                    (xm[2 * ii + 0][countr] * (1.0 / h) * xm[2 * ii + 0][countr] +
                    xm[2 * ii + 2][countr] * (1.0 / h) * xm[2 * ii + 2][countr] +
                    xm[2 * ii + 0][countr] * (-1.0 / h) * xm[2 * ii + 2][countr] +
                    xm[2 * ii + 2][countr] * (-1.0 / h) * xm[2 * ii + 0][countr]) *
                    (rhov[ii] * (vpv[ii] * vpv[ii] - 2 * vsv[ii] * vsv[ii]));

        //% mu sensitivity
        //% k^2 term
                snsmf[ii][countr] = (1 / (2 * vp[countr] * U[countr])) *
                    (xm[2 * ii - 1][countr] * (2 * h / 3) * xm[2 * ii - 1][countr] +
                    xm[2 * ii + 1][countr] * (2 * h / 3) * xm[2 * ii + 1][countr] +
                    xm[2 * ii + 1][countr] * (h / 3) * xm[2 * ii - 1][countr] +
                    xm[2 * ii - 1][countr] * (h / 3) * xm[2 * ii + 1][countr] +
                    xm[2 * ii][countr] * (h / 3) * xm[2 * ii][countr] +
                    xm[2 * ii + 2][countr] * (h / 3) * xm[2 * ii + 2][countr] +
                    xm[2 * ii][countr] * (h / 6) * xm[2 * ii + 2][countr] +
                    xm[2 * ii + 2][countr] * (h / 6) * xm[2 * ii][countr]) *
                    (rhov[ii] * (vsv[ii] * vsv[ii]));

        // % k^1 term
            snsmf[ii][countr] += (1 / (2 * (2 * pi * fks[countr]) * U[countr])) *
                (xm[2 * ii - 1][countr] * (-1.0 / 2) * xm[2 * ii][countr] +
                 xm[2 * ii - 1][countr] * (-1.0 / 2) * xm[2 * ii + 2][countr] +
                 xm[2 * ii][countr] * (-1.0 / 2) * xm[2 * ii - 1][countr] +
                 xm[2 * ii][countr] * (1.0 / 2) * xm[2 * ii + 1][countr] +
                 xm[2 * ii + 1][countr] * (1.0 / 2) * xm[2 * ii][countr] +
                 xm[2 * ii + 1][countr] * (1.0 / 2) * xm[2 * ii + 2][countr] +
                 xm[2 * ii + 2][countr] * (-1.0 / 2) * xm[2 * ii - 1][countr] +
                 xm[2 * ii + 2][countr] * (1.0 / 2) * xm[2 * ii + 1][countr]) *
                (rhov[ii] * (vsv[ii] * vsv[ii]));

        //% k^0 term
            snsmf[ii][countr] += (vp[countr] / (2 * std::pow(2 * pi * fks[countr], 2) * U[countr])) *
                (xm[2 * ii + 0][countr] * (2.0 / h) * xm[2 * ii + 0][countr] +
                 xm[2 * ii + 2][countr] * (2.0 / h) * xm[2 * ii + 2][countr] +
                 xm[2 * ii + 0][countr] * (-2.0 / h) * xm[2 * ii + 2][countr] +
                 xm[2 * ii + 2][countr] * (-2.0 / h) * xm[2 * ii + 0][countr] +
                 xm[2 * ii - 1][countr] * (1.0 / h) * xm[2 * ii - 1][countr] +
                 xm[2 * ii + 1][countr] * (1.0 / h) * xm[2 * ii + 1][countr] +
                 xm[2 * ii - 1][countr] * (-1.0 / h) * xm[2 * ii + 1][countr] +
                 xm[2 * ii + 1][countr] * (-1.0 / h) * xm[2 * ii - 1][countr]) *
                (rhov[ii] * (vsv[ii] * vsv[ii]));

        //% thickness sensitivity
        //% omega^2 term
            snsmfh[ii][countr] = -(vp[countr] / (2 * U[countr])) *
                    (xm[2 * ii - 1][countr] * (rhov[ii] / 2) * xm[2 * ii - 1][countr] +
                    xm[2 * ii + 0][countr] * (rhov[ii] / 2) * xm[2 * ii + 0][countr] +
                    xm[2 * ii + 1][countr] * (rhov[ii] / 2) * xm[2 * ii + 1][countr] +
                    xm[2 * ii + 2][countr] * (rhov[ii] / 2) * xm[2 * ii + 2][countr]) *
                    h;
        //% k^2 term
            double pmod = rhov[ii] * vpv[ii] * vpv[ii];
            double smod = rhov[ii] * vsv[ii] * vsv[ii];

            snsmfh[ii][countr] += (1 / (2 * vp[countr] * U[countr])) *
                (xm[2 * ii - 1][countr] * (pmod / 3) * xm[2 * ii - 1][countr] +
                 xm[2 * ii + 1][countr] * (pmod / 3) * xm[2 * ii + 1][countr] +
                 xm[2 * ii + 1][countr] * (pmod / 6) * xm[2 * ii - 1][countr] +
                 xm[2 * ii - 1][countr] * (pmod / 6) * xm[2 * ii + 1][countr] +
                 xm[2 * ii][countr] * (smod / 3) * xm[2 * ii][countr] +
                 xm[2 * ii + 2][countr] * (smod / 3) * xm[2 * ii + 2][countr] +
                 xm[2 * ii][countr] * (smod / 6) * xm[2 * ii + 2][countr] +
                 xm[2 * ii + 2][countr] * (smod / 6) * xm[2 * ii][countr]) *
                h;
        // % k^0 term
            snsmfh[ii][countr] += (vp[countr] / (2 * std::pow(2 * pi * fks[countr], 2) * U[countr])) *
                (xm[2 * ii + 0][countr] * (pmod) * (-1.0 / (h * h)) * xm[2 * ii + 0][countr] +
                 xm[2 * ii + 2][countr] * (pmod) * (-1.0 / (h * h)) * xm[2 * ii + 2][countr] +
                 xm[2 * ii + 0][countr] * (-pmod) * (-1.0 / (h * h)) * xm[2 * ii + 2][countr] +
                 xm[2 * ii + 2][countr] * (-pmod) * (-1.0 / (h * h)) * xm[2 * ii + 0][countr] +
                 xm[2 * ii - 1][countr] * (smod) * (-1.0 / (h * h)) * xm[2 * ii - 1][countr] +
                 xm[2 * ii + 1][countr] * (smod) * (-1.0 / (h * h)) * xm[2 * ii + 1][countr] +
                 xm[2 * ii - 1][countr] * (-smod) * (-1.0 / (h * h)) * xm[2 * ii + 1][countr] +
                 xm[2 * ii + 1][countr] * (-smod) * (-1.0 / (h * h)) * xm[2 * ii - 1][countr]) *
                h;
    }//fks loop ends
}//ii loop ends

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% the shear velocity phase sensitivity kernel for frequencies of interest
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//% sensitivity for Vp fixed or Vp/Vs fixed

double pratioflag; //comes from function call
vector<vector<double>> snsmf_vs (Nn,vector<double>(fks.size()));
if (pratioflag == 0) {
        for (int ii = 0; ii < Nn; ii++) {
            for (int countr = 0; countr < fks.size(); countr++) {
                double pratio = vsv[ii] * vsv[ii] / (vpv[ii] * vpv[ii] - 2 * vsv[ii] * vsv[ii]);
                snsmf_vs[ii][countr] = 2 * snsmf[ii][countr] - 4 * snsmflam[ii][countr] * pratio;
            }
        }
    }
    else if (pratioflag == 1) {
        for (int ii = 0; ii < Nn; ii++) {
            for (int countr = 0; countr < fks.size(); countr++) {
                snsmf_vs[ii][countr] = 2 * (snsmf[ii][countr] + snsmflam[ii][countr]);
            }
        }
    }
    else {
        // Handle the case when pratioflag is neither 0 nor 1
    }

vector<double> vfi(109,0);
//make a vector of the frequency of interest
if(vflg[0]==0){
    vfi[0]=1;
}else{
    vfi[0]=2;
}
for(int ii=1;ii<vflg.size();ii++){
    if (vflg[ii] == 1 & vflg[ii-1] == 0)
        vfi[ii] = vfi[ii-1] + 2;

    else if (vflg[ii] == 1 & vflg[ii-1] == 1)
        vfi[ii] = vfi[ii-1] + 3;

    else if (vflg[ii] == 0 & vflg[ii-1] == 0)
        vfi[ii] = vfi[ii-1] + 1;

    else
        vfi[ii] = vfi[ii-1] + 2;
}
countr = 0;;
vector<double> Uf(109);
vector<vector<double>> snsmf_vstotf (Nn,vector<double>(fks.size()));
vector<vector<double>> snsmf_htotf (Nn,vector<double>(fks.size()));
vector<vector<double>> snsmf_vstot (Nn,vector<double>(fks.size()));
vector<vector<double>> snsmf_htot (Nn,vector<double>(fks.size()));
int ii;
    for (double f : fks_orig) {
        countr++;

        if (vflg[countr - 1] == 0) {
            Uf[countr - 1] = vp[vfi[countr - 1]];

            for (ii = 0; ii < Nn; ii++) {
                snsmf_vstotf[ii][countr - 1] = vp[vfi[countr - 1]] * snsmf_vs[ii][vfi[countr - 1]] / vsv[ii];
                snsmf_htotf[ii][countr - 1] = vp[vfi[countr - 1]] * snsmfh[ii][vfi[countr - 1]] / hv[ii];
            }
        }
        else {
            if (vfi[countr - 1] > 0 && vfi[countr - 1] < fks.size() - 1) {
                double df = fks[vfi[countr - 1] + 1] - fks[vfi[countr - 1] - 1];
                double derivative = (snsmf_vs[ii][vfi[countr - 1] + 1] - snsmf_vs[ii][vfi[countr - 1] - 1]) / (df * 2 * M_PI);

                for (int ii = 0; ii < Nn; ii++) {
                    snsmf[ii][countr - 1] = snsmf_vs[ii][vfi[countr - 1]] +
                        (U[vfi[countr - 1]] / vp[vfi[countr - 1]]) * (2 * M_PI * fks[vfi[countr - 1]] * derivative);

                    snsmf_vstot[ii][countr - 1] = snsmf[ii][countr - 1];
                    snsmf_htot[ii][countr - 1] = snsmfh[ii][vfi[countr - 1]];

                    snsmf_vstotf[ii][countr - 1] = vp[vfi[countr - 1]] * snsmf_vstot[ii][countr - 1] / vsv[ii];
                    snsmf_htotf[ii][countr - 1] = U[vfi[countr - 1]] * snsmf_htot[ii][countr - 1] / hv[ii];
                }
            }
            else {
                for (int ii = 0; ii < Nn; ii++) {
                    snsmf[ii][countr - 1] = snsmf_vs[ii][vfi[countr - 1]];
                    snsmf_vstot[ii][countr - 1] = snsmf[ii][countr - 1];
                    snsmf_htot[ii][countr - 1] = snsmfh[ii][vfi[countr - 1]];
                    snsmf_vstotf[ii][countr - 1] = vp[vfi[countr - 1]] * snsmf_vstot[ii][countr - 1] / vsv[ii];
                    snsmf_htotf[ii][countr - 1] = U[vfi[countr - 1]] * snsmf_htot[ii][countr - 1] / hv[ii];
                }
            }
        }
    }
return 0;
}

