#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <TGraph.h>
#include "TArtEventStore.hh"

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

void generatetree(std::__cxx11::string infile, std::__cxx11::string output);
//void progressbar(int currevent, int totevent, int offset, int barwidth=30);

int Obertelli_filter(std::vector<double> &x,std::vector<double> &y,std::vector<double> &q,
                     std::vector<double> &x_out,std::vector<double> &y_out,
                     std::vector<double> &q_out, std::vector<bool> &ringbool);
double conv_fit(double *x, double *p);

void Hough_filter(std::vector<double> &x,std::vector<double> &y,std::vector<double> &z,
                  std::vector<double> &q,std::vector<double> &x_out,std::vector<double> &y_out,
                  std::vector<double> &z_out,std::vector<double> &q_out);

double FitFunction(double *x, double *p);

void FindStart(std::vector<double> pStart, std::vector<double> chi,
               std::vector<int> fitstatus, TGraph *grxz, TGraph *gryz);

void SumDistance1(int &, double *, double &sum, double *par, int);

double distancelinepoint(double &x, double &y, double &z, double *p);

void SumDistance2(int &, double *, double &sum, double *par, int);

void vertex(std::vector<double> &p, std::vector<double> &pp, double &xv,
            double &yv, double &zv);

uint getset(TArtEventStore &estore);