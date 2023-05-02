//
// Created by Filippo Casari on 01/05/23.
//
#include <vector>
#include "Particle.h"
#include <fstream>
#include <iostream>
#ifndef AS3PM_RDF_H
#define AS3PM_RDF_H
using namespace std;
void computeRdf(vector<double> &rdf, vector<Particle> particles, double L, double rc, int n_bins){
    //vector<double> rdf(n_bins);
    double r = 0;
    int count = 0;
    double d_r =rc/n_bins;
    rdf[0] = 0;
    count++;
    r += d_r;
    while(r<=rc){
        rdf[count] = 0;
        for(int i=0; i<particles.size(); i++){
            for(int j=0; j<particles.size(); j++){
                if(i!=j){
                    double dx = particles[i].x - particles[j].x;
                    double dy = particles[i].y - particles[j].y;
                    if(dx > L/2){
                        dx -= L;
                    }else if(dx < -L/2){
                        dx += L;
                    }
                    if(dy > L/2){
                        dy -= L;
                    }else if(dy < -L/2){
                        dy += L;
                    }
                    double d = sqrt(dx*dx + dy*dy);
                    if(d<=r+d_r && d>r){
                        rdf[count] += 1;
                    }
                }
            }
        }

        rdf[count] /= (r *r* M_PI* d_r * (particles.size()/ L*L));

        count++;
        r += d_r;
    }

    //return &rdf;
}
#endif //AS3PM_RDF_H
