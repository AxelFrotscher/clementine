//
// Created by afrotscher on 1/7/19.
//

#pragma once

#include <vector>
#include <iostream>
#include "TGraph.h"
#include <mutex>
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
    double r_vertex;
    double thetaz1;
    double thetaz2;
    double phi_vertex;
    double trackNbr;
    double trackNbr_final;

    TMinosPass(double r_vertex_, double thetaz1_, double thetaz2_,
               double phi_vertex_, double trackNbr_, double trackNbr_Final):
               r_vertex(r_vertex_), thetaz1(thetaz1_), thetaz2(thetaz2_),
               phi_vertex(phi_vertex_),trackNbr(trackNbr_),
               trackNbr_final(trackNbr_Final) {}
};

class minosana{
public:
    minosana(int filled_, double TShaping_, double TimeBinElec_,
             double DelayTrigger_, double VDrift_,
             vector<vector<double>> minostrackxy_,
             vector<vector<double>> minoscalibvalues_, vector<double> xpad_,
             vector<double> ypad_, vector<double> qpad_, int threadno_,
             vector<TH2D> &minossingleevent_):
             filled(filled_), Tshaping(TShaping_), TimeBinElec(TimeBinElec_),
             DelayTrigger(DelayTrigger_), VDrift(VDrift_),
             minostrackxy(minostrackxy_), minoscalibvalues(minoscalibvalues_),
             Xpad(xpad_), Ypad(ypad_),Qpad(qpad_), threadno(threadno_),
             minossingleevent(minossingleevent_){

    }
    TMinosResult getTMinosResult(){return dataresult;}
    TMinosPass analyze();

private:
    int Obertelli_filter(vector<double> &x,vector<double> &y,vector<double> &q,
                         vector<double> &x_out,vector<double> &y_out,
                         vector<double> &q_out, vector<bool> &ringbool);
    void Hough_filter(vector<double> &x,vector<double> &y,vector<double> &z,
                      vector<double> &q,vector<double> &x_out,vector<double> &y_out,
                      vector<double> &z_out,vector<double> &q_out);
    void FindStart(vector<double> pStart, vector<double> chi, vector<int> fitstatus,
                             TGraph *grxz, TGraph *gryz);
    void vertex(vector<double> &p, vector<double> &pp, double &xv,
                          double &yv, double &zv);
    void debug();

    vector<double> Xpad, Ypad, Qpad, Xpadnew,Ypadnew, Qpadnew, Zpadnew;

    double z_vertex =0, x_vertex=0, y_vertex =0, r_vertex =0, phi_vertex=0,
            thetaz1=0, thetaz2=0;

    int filled =0 ;
    double Tshaping =0;
    double TimeBinElec = 0;
    double DelayTrigger =0;
    double VDrift = 0;
    double threadno=0;
    vector<vector<double>> minoscalibvalues;
    vector<vector<double>> minostrackxy;
    vector<TH2D> &minossingleevent;

    TMinosClust fitdata;
    TMinosResult dataresult;
    static std::mutex minos5;
};

// Ugly part outside of the class, because member-functions are trickier than
// regular functions

double FitFunction(double *x, double *p);
void SumDistance2(int &, double *, double &sum, double *par, int);
void SumDistance1(int &, double *, double &sum, double *par, int);
double distancelinepoint(double &x, double &y, double &z, double *p);
double conv_fit(double *x, double *p);
extern TMinosResult tmr;