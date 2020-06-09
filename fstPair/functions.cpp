//
//  functions.cpp
//  fstCalculation
//
//  Created by Gang Peng on 6/8/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#include "functions.hpp"
#include "normal.h"

#include <fstream>
#include <iostream>

using namespace std;

std::map<std::string, int>  getSampleSize(std::string filename){
    map<string, int> rlt;
    
    ifstream fin(filename.c_str());
    if(!fin.is_open()){
        cout << "Cannot open file " << filename << ". Please check the name of sample size file." << endl;
        rlt.clear();
        return rlt;
    }
    
    string header;
    getline(fin, header);
    //vector<string> vsHeader = split(header, "\t");
    
    while(!fin.eof()){
        string fline;
        getline(fin, fline);
        if(fline.size() < 2){
            break;
        }
        vector<string> vsLine = split(fline, "\t");
        rlt[vsLine[0]] = stoi(vsLine[2]);
    }
    fin.close();
    
    return rlt;
}


double fstHudson(double p1, double p2){
    double N = (p1-p2)*(p1-p2);
    double D = p1*(1.0-p2) + p2*(1.0-p1);
    if(D==0){
        return 0;
    } else {
        return N/D;
    }
}

double fstHudson(double p1, double p2, int s1, int s2){
    double D = p1*(1.0-p2) + p2*(1.0-p1);
    if(D==0){
        return 0;
    }
    double N = (p1-p2)*(p1-p2) - p1*(1-p1)/(s1-1.0) - p2*(1-p2)/(s2-1.0);
    return N/D;
}

double fstNHudson(double p1, double p2){
    return (p1-p2)*(p1-p2);
}

double fstNHudson(double p1, double p2, int s1, int s2){
    return (p1-p2)*(p1-p2) - p1*(1-p1)/(s1-1.0) - p2*(1-p2)/(s2-1.0);
}

double fstDHudson(double p1, double p2){
    return p1*(1.0-p2) + p2*(1.0-p1);
}

