//
// Created by afrotscher on 1/7/19.
//

#pragma once

#include <vector>
#include <iostream>
#include "TGraph.h"
#include <mutex>
#include "TH3D.h"
#include "TH2D.h"
using std::vector;

struct TMinosClust{
    std::vector<double> x_mm;
    std::vector<double> y_mm;
    std::vector<double> t_ns;
    std::vector<double> z_mm;
    std::vector<double> Chargemax;
    std::vector<double> n_Cluster;
    std::vector<double> n_Pads;
    std::vector<double> Chi2;

    void add(double xmm, double ymm, double tns, double zmm, double chargemax,
             double ncluster, double npads, double chi2){
        x_mm.push_back(xmm);
        y_mm.push_back(ymm);
        t_ns.push_back(tns);
        z_mm.push_back(zmm);
        Chargemax.push_back(chargemax);
        n_Cluster.push_back(ncluster);
        n_Pads.push_back(npads);
        Chi2.push_back(chi2);
    }

    void replace(double xmm, double ymm, double tns, double zmm, double chargemax,
                 double ncluster, double npads, double chi2, double index){
        if(index >= x_mm.size()){
            std::cout << "Replacing in TMinosClust Event " << index + 1
                      << " but only " << x_mm.size() << " Elements!!." <<std::endl;
            std::__throw_invalid_argument("Index out of bounds\n");
        }
        x_mm.at(index) = xmm;
        y_mm.at(index) = ymm;
        t_ns.at(index) = tns;
        z_mm.at(index) = zmm;
        Chargemax.at(index) = chargemax;
        n_Cluster.at(index) = ncluster;
        n_Pads.at(index) = npads;
        Chi2.at(index) = chi2;
    }

    TMinosClust() = default;
    void reset(){
        x_mm = y_mm = t_ns = z_mm = Chargemax = n_Cluster = n_Pads = Chi2 = {0};
    }
};

struct TMinosResult{
    std::vector<double> x_mm;
    std::vector<double> y_mm;
    std::vector<double> z_mm;
    std::vector<double> Chargemax;
    std::vector<double> n_Cluster;
    std::vector<double> n_Pads;
    std::vector<double> Chi2;

    void add(double xmm, double ymm, double zmm, double chargemax,
             double ncluster, double npads, double chi2){
        x_mm.push_back(xmm);
        y_mm.push_back(ymm);
        z_mm.push_back(zmm);
        Chargemax.push_back(chargemax);
        n_Cluster.push_back(ncluster);
        n_Pads.push_back(npads);
        Chi2.push_back(chi2);
    }

    void replace(double xmm, double ymm, double zmm, double chargemax,
                 double ncluster, double npads, double chi2, double index){
        if(index >= x_mm.size()){
            std::cout << "Replacing in TMinosResult Event " << index + 1
                      << " but only " << x_mm.size() << " Elements!!." <<std::endl;
            std::__throw_invalid_argument("Index out of bounds\n");
        }
        x_mm.at(index) = xmm;
        y_mm.at(index) = ymm;
        z_mm.at(index) = zmm;
        Chargemax.at(index) = chargemax;
        n_Cluster.at(index) = ncluster;
        n_Pads.at(index) = npads;
        Chi2.at(index) = chi2;
    }

    TMinosResult() = default;

    void reset(){
        x_mm = y_mm = z_mm = Chargemax = n_Cluster = n_Pads = Chi2 = {0};
    }
};

struct TMinosPass{
    double r_vertex = 0;            // distance between reaction vertex and beam
    double x_vertex = 0;            // x-position of vertex
    double y_vertex = 0;            // y-position of vertex
    double z_vertex = 0;            // z-position of vertex
    vector<double> theta{};         // beam angles for each track
    vector<double> lambda{};        // interangles between tracks
    int trackNbr = 0;               // number of tracks pre-3D-hough
    int trackNbr_final = 0;         // final number of tracks

    vector<double> phi2d{};         // projected angles to x-y plane
    vector<double> chargeweight{};  // mean charge per track (pre length)
    vector<double> vertexdist{};    // distacne between the 3 vertices
    vector<double> phi2dE{};        // projected angles error
    vector<double> thetaerr{};      // theta error

    TMinosPass() = default;
};

struct minos_internal{  // MInos iNTernal values
    std::array<double, 4> pStart{0,1,0,1};   // start parameters for the 2D fit
    std::array<double, 4> parFit;   // fitted parameters for the 2D fit
    std::array<double, 4> parFit_r;   // fitted parameters for the 2D fit
    std::array<double, 4> err;      // uncertainties on all parameters
    std::array<double, 2> chi2;     // chi² for the fit
    double covxy;                   // covariance of slope Ax and Ay
    std::array<int, 2> fitStatus;   // goodness of the fit
    std::array<double, 10> arglist;
    
    minos_internal() = default;
};

class minosana{
public:
    minosana(int filled_, double TShaping_, double TimeBinElec_,
             double DelayTrigger_, double VDrift_,
             vector<vector<double>> *minostrackxy_,
             vector<vector<double>> *minoscalibvalues_,
             vector<vector<double>> *minostime_, vector<double> *xpad_,
             vector<double> *ypad_, vector<double> *qpad_, int threadno_,
             vector<TH2C> &minossingleevent_):
             filled(filled_), Tshaping(TShaping_), TimeBinElec(TimeBinElec_),
             DelayTrigger(DelayTrigger_), VDrift(VDrift_),
             minostrackxy(*minostrackxy_), minoscalibvalues(*minoscalibvalues_),
             minostime(*minostime_), Xpad(*xpad_), Ypad(*ypad_),Qpad(*qpad_),
             threadno(threadno_), minossingleevent(minossingleevent_){
    }
    TMinosResult getTMinosResult(){return dataresult;}
    TMinosPass analyze();

    static int Obertelli_filter(vector<double> &x,vector<double> &y,vector<double> &q,
                                vector<double> &x_out,vector<double> &y_out,
                                vector<double> &q_out, vector<bool> &ringbool);

private:
    static void Hough_filter(vector<double> &x,vector<double> &y,vector<double> &z,
                      vector<double> &q,vector<double> &x_out,vector<double> &y_out,
                      vector<double> &z_out,vector<double> &q_out);
    void FindStart(minos_internal &mi_nt, TGraph &grxz, TGraph &gryz);
    static void vertex(const std::array<double, 4> &p, const std::array<double, 4> &pp,
                       double &xv, double &yv, double &zv);
    static double distancelineline(const std::array<double,4> &l1,
                                   const std::array<double,4> &l2);
    void debug();
    static std::array<double, 4> rotatesp(const double &rot,
                                   const std::array<double,4> &initialvector);
    
    int filled;
    double Tshaping;
    double TimeBinElec;
    double DelayTrigger;
    double VDrift;
    vector<vector<double>> minostrackxy;
    vector<vector<double>> minoscalibvalues;
    vector<vector<double>> minostime;
    vector<double> Xpad, Ypad, Qpad, Xpadnew,Ypadnew, Qpadnew, Zpadnew;
    
    vector<minos_internal> mint;
    
    int threadno=0;
    vector<TH2C> &minossingleevent;
    TMinosClust fitdata;
    TMinosResult dataresult;
};

// Ugly part outside of the class, because member-functions are trickier than
// regular functions

double distancelinepoint(double x, double y, double z, double *p);

thread_local inline TMinosResult tmr;
thread_local inline int mode = 0;