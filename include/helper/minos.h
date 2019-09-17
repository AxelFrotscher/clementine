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
    double r_vertex;                 // distance between reaction vertex and beam
    vector<double> thetaz;           // beam angles for each track
    vector<double> lambda;           // interangles between tracks
    double trackNbr;                 // number of tracks pre-3D-hough
    double trackNbr_final;           // final number of tracks
    double z_vertex;                 // z-position of vertex
    vector<double> phi2d;         // projected angles to x-y plane
    vector<double> chargeweight;     // mean charge per track (pre length)
    vector<double> vertexdist;       // distacne between the 3 vertices
    vector<double> phi2dE;        // projected angles error
    vector<double> thetaerr;         // theta error

    TMinosPass(double r_vertex_, vector<double> &thetaz_,
               vector<double> &lambda_, double trackNbr_, double trackNbr_Final,
               double z_vertex_, vector<double> &phi2d_,
               vector<double> &chargeweight_, vector<double> &vertexdist_,
               vector<double> &phi2dE_, vector<double> &thetaerr_):
            r_vertex(r_vertex_), thetaz(thetaz_), lambda(lambda_),
            trackNbr(trackNbr_), trackNbr_final(trackNbr_Final),
            z_vertex(z_vertex_), phi2d(phi2d_), chargeweight(chargeweight_),
            vertexdist(vertexdist_), phi2dE(phi2dE_), thetaerr(thetaerr_){}
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
    void FindStart(vector<double> &pStart, vector<double> &chi, vector<int> &fitstatus,
                             TGraph &grxz, TGraph &gryz);
    static void vertex(const vector<double> &p, const vector<double> &pp, double &xv,
                          double &yv, double &zv);
    static double distancelineline(vector<double> &l1, vector<double> &l2);
    void debug();
    static vector<double> rotatesp(double &rot, vector<double> &initialvector);

    double z_vertex =0, x_vertex=0, y_vertex =0, r_vertex =0;
    
    int filled =0 ;
    double Tshaping;
    double TimeBinElec = 0;
    double DelayTrigger;
    double VDrift = 0;
    vector<vector<double>> minostrackxy;
    vector<vector<double>> minoscalibvalues;
    vector<vector<double>> minostime;
    vector<double> Xpad, Ypad, Qpad, Xpadnew,Ypadnew, Qpadnew, Zpadnew;
    
    int threadno=0;
    vector<TH2C> &minossingleevent;
    vector<double> theta{0,0,0};
    vector<double> thetaerr{0,0,0};
    vector<double> chargeweight;
    vector<double> lambda={};
    TMinosClust fitdata;
    TMinosResult dataresult;
    //inline static std::mutex minos5;
};

// Ugly part outside of the class, because member-functions are trickier than
// regular functions

double distancelinepoint(double x, double y, double z, double *p);

thread_local inline TMinosResult tmr;
thread_local inline int mode = 0;