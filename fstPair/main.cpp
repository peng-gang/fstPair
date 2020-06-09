//
//  main.cpp
//  fstPair
//
//  Created by Gang Peng on 6/8/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>

#include "normal.h"
#include "functions.hpp"

using namespace std;

int main(int argc, char ** argv) {
/*
 * calculation fst between pairwise race groups
 * -i input.txt // input file containing MAF informaiton
 * -s sampleSize.txt // text file including sample size
 * -o output.txt // output file
 * -r race1 race2 race3 ... // race groups to compare. these groups should be in the first row (header) the input file. -r all: comparing all pairs
 */
    
    vector<string> mustOptions = {"-i", "-s", "-o", "-r"};
    vector<string> allOptions = {"-i", "-s", "-o", "-r"};
    
    map<string, vector<string> > cmLine = parseCMLine(argc, argv, allOptions, mustOptions);
    if(cmLine.size()==0){
        return -1;
    }
    
    string fnameIn = cmLine["-i"][0];
    string fnameOut = cmLine["-o"][0];
    string fnameSS = cmLine["-s"][0];
    
    ifstream fin(fnameIn.c_str());
    if(!fin.is_open()){
        cout << "Cannot open file " << fnameIn << ". Please check the input file name." << endl;
        return -1;
    }
    
    ofstream fout(fnameOut.c_str());
    if(!fout.is_open()){
        cout << "Cannot open file " << fnameOut << ". Please check the output file name." << endl;
        return -1;
    }
    
    map<string, int> sampleSize = getSampleSize(fnameSS);
    if(sampleSize.size()==0){
        return -1;
    }
    
    string header;
    getline(fin, header);
    
    //check header and race in the header
    vector<string> vsHeader = split(header, "\t");
    
    vector<size_t> idxRace;
    if(cmLine["-r"].size()==1){
        if(cmLine["-r"][0] == "all"){
            idxRace.clear();
            // There is a bug here
            // I assume race info from 7 to 11 // 6-10
            for(size_t i=6; i<11; i++){
                idxRace.push_back(i);
            }
        } else {
            cout << "-r must include more than 1 race group to take pairwise comparision unless you use \"-r all\"" << endl;
            fin.close();
            return -1;
        }
    } else {
        for(size_t i=0; i<cmLine["-r"].size(); i++){
            vector<string>::iterator it = find(vsHeader.begin(), vsHeader.end(), cmLine["-r"][i]);
            if(it!=vsHeader.end()){
                idxRace.push_back(distance(vsHeader.begin(), it));
            } else {
                cout << "Cannot find race \"" << cmLine["-r"][i] << "\" in the header of the input file. Plese take a look at the header of the input file." << endl;
                fin.close();
                return -1;
            }
        }
    }
    
    size_t numCompare = idxRace.size()*(idxRace.size()-1)/2;
    
    //find sample size for each race
    vector<int> ss;
    for(size_t i=0; i<idxRace.size(); i++){
        map<string, int>::iterator it = sampleSize.find(vsHeader[idxRace[i]].substr(3));
        if(it != sampleSize.end()){
            ss.push_back(it->second);
        } else {
            cout << "Cannot find sample size for " << vsHeader[idxRace[i]] << ". Please check sample size file." << endl;
            fin.close();
            return -1;
        }
    }
    
    cout << "Start" <<endl;
    long long numVar = 0;
    
    // results for allele exist in at least one race ethnicity group
    vector<double> NAll1 (numCompare, 0);
    vector<double> N2All1 (numCompare, 0);
    vector<double> DAll1 (numCompare, 0);
    vector<double> fstAll1 (numCompare, 0);
    vector<double> fst2All1 (numCompare, 0);
    vector<int> count1(numCompare, 0);
    
    
    // results for allele exist in all groups
    vector<double> NAllAll (numCompare, 0);
    vector<double> N2AllAll (numCompare, 0);
    vector<double> DAllAll (numCompare, 0);
    vector<double> fstAllAll (numCompare, 0);
    vector<double> fst2AllAll (numCompare, 0);
    vector<int> countAll(numCompare, 0);

    // results for allele exist in afr
    vector<double> NAllAfr (numCompare, 0);
    vector<double> N2AllAfr (numCompare, 0);
    vector<double> DAllAfr (numCompare, 0);
    vector<double> fstAllAfr (numCompare, 0);
    vector<double> fst2AllAfr (numCompare, 0);
    vector<int> countAfr(numCompare, 0);
    
    
    
    //std::mt19937 gen(123); //seed
    srand(123);
    while(!fin.eof()){
        string fline;
        getline(fin, fline);
        if(fline.size() < 2){
            break;
        }
        
        numVar++;
        if(numVar % 1000000 == 0){
            cout << numVar << endl;
        }
        
        vector<string> vsLine = split(fline, "\t");
        
        
        //std::uniform_real_distribution<> dis(0, 1.0);
        //double tmp = dis(gen);
        double tmp = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        bool flagRecord = false;
        if(tmp < 0.011){
            flagRecord = true;
        }
            
        if(flagRecord){
            fout << vsLine[0] << "\t" << vsLine[1] << "\t" << vsLine[11] << "\t" << vsLine[13];
        }
        
        vector<double> maf;
        for(size_t i=0; i<idxRace.size(); i++){
            if(vsLine[idxRace[i]] == "NA"){
                maf.push_back(-1.0);
            } else {
                maf.push_back(stof(vsLine[idxRace[i]]));
            }
        }
        
        // allele exist in at least one race ethnicity group
        bool flag1 = false;
        for(size_t i=0; i<maf.size(); i++){
            if(maf[i] > 0){
                flag1=true;
                break;
            }
        }
        
        if(flag1){
            int idx=0;
            for(size_t i=0; i<(idxRace.size()-1); i++){
                for(size_t j=(i+1); j<idxRace.size(); j++){
                    if(vsLine[idxRace[i]] != "NA" && vsLine[idxRace[j]] != "NA"){
                        double maf1 = stof(vsLine[idxRace[i]]);
                        double maf2 = stof(vsLine[idxRace[j]]);
                        
                        double N = fstNHudson(maf1, maf2);
                        double N2 = fstNHudson(maf1, maf2, ss[i], ss[j]);
                        double D = fstDHudson(maf1, maf2);
                        
                        double fst = 0;
                        double fst2 = 0;
                        if(D!=0){
                            fst = N/D;
                            fst2 = N2/D;
                        }
                        
                        NAll1[idx] = NAll1[idx] + N;
                        N2All1[idx] = N2All1[idx] + N2;
                        DAll1[idx] = DAll1[idx] + D;
                        fstAll1[idx] = fstAll1[idx] + fst;
                        fst2All1[idx] = fst2All1[idx] + fst2;
                        count1[idx]++;
                        idx++;
                        
                        if(flagRecord){
                            fout << "\t" << N << "\t" << N2 << "\t" << D;
                        }
                       
                    } else {
                        if(flagRecord){
                            fout <<"\tNA\tNA\tNA";
                        }
                    }
                }
            }
        } else {
            if(flagRecord){
                for(size_t i=0; i<numCompare; i++){
                    fout <<"\tNA\tNA\tNA";
                }
            }
        }
        
        // allele exist in all groups
        bool flagAll = true;
        for(size_t i=0; i<maf.size(); i++){
            if(maf[i] == 0){
                flagAll = false;
                break;
            }
        }
        
        
        if(flagAll){
            int idx=0;
            for(size_t i=0; i<(idxRace.size()-1); i++){
                for(size_t j=(i+1); j<idxRace.size(); j++){
                    if(vsLine[idxRace[i]] != "NA" && vsLine[idxRace[j]] != "NA"){
                        double maf1 = stof(vsLine[idxRace[i]]);
                        double maf2 = stof(vsLine[idxRace[j]]);
                        
                        double N = fstNHudson(maf1, maf2);
                        double N2 = fstNHudson(maf1, maf2, ss[i], ss[j]);
                        double D = fstDHudson(maf1, maf2);
                        
                        double fst = 0;
                        double fst2 = 0;
                        if(D!=0){
                            fst = N/D;
                            fst2 = N2/D;
                        }
                        
                        NAllAll[idx] = NAllAll[idx] + N;
                        N2AllAll[idx] = N2AllAll[idx] + N2;
                        DAllAll[idx] = DAllAll[idx] + D;
                        fstAllAll[idx] = fstAllAll[idx] + fst;
                        fst2AllAll[idx] = fst2AllAll[idx] + fst2;
                        countAll[idx]++;
                        idx++;
                        
                        if(flagRecord){
                            fout << "\t" << N << "\t" << N2 << "\t" << D;
                        }
                       
                    } else {
                        if(flagRecord){
                            fout <<"\tNA\tNA\tNA";
                        }
                    }
                }
            }
        } else {
            if(flagRecord){
                for(size_t i=0; i<numCompare; i++){
                    fout <<"\tNA\tNA\tNA";
                }
            }
        }
        
        // allele exist in afr
        bool flagAfr = true;
        if(vsLine[6] != "NA"){
            if(stof(vsLine[6])==0){
                flagAfr = false;
            }
        }
        
        if(flagAfr){
            int idx=0;
            for(size_t i=0; i<(idxRace.size()-1); i++){
                for(size_t j=(i+1); j<idxRace.size(); j++){
                    if(vsLine[idxRace[i]] != "NA" && vsLine[idxRace[j]] != "NA"){
                        double maf1 = stof(vsLine[idxRace[i]]);
                        double maf2 = stof(vsLine[idxRace[j]]);
                        
                        double N = fstNHudson(maf1, maf2);
                        double N2 = fstNHudson(maf1, maf2, ss[i], ss[j]);
                        double D = fstDHudson(maf1, maf2);
                        
                        double fst = 0;
                        double fst2 = 0;
                        if(D!=0){
                            fst = N/D;
                            fst2 = N2/D;
                        }
                        
                        NAllAfr[idx] = NAllAfr[idx] + N;
                        N2AllAfr[idx] = N2AllAfr[idx] + N2;
                        DAllAfr[idx] = DAllAfr[idx] + D;
                        fstAllAfr[idx] = fstAllAfr[idx] + fst;
                        fst2AllAfr[idx] = fst2AllAfr[idx] + fst2;
                        countAfr[idx]++;
                        idx++;
                        
                        if(flagRecord){
                            fout << "\t" << N << "\t" << N2 << "\t" << D;
                        }
                       
                    } else {
                        if(flagRecord){
                            fout <<"\tNA\tNA\tNA";
                        }
                    }
                }
            }
        } else {
            if(flagRecord){
                for(size_t i=0; i<numCompare; i++){
                    fout <<"\tNA\tNA\tNA";
                }
            }
        }
        
        if(flagRecord){
            fout<<endl;
        }
    }
    
    fin.close();
    fout.close();
    
    // 1
    cout << endl << endl;
    for(size_t i=0; i<numCompare; i++){
        cout << setprecision(10) << NAll1[i] << "\t" << N2All1[i] << "\t" << DAll1[i] << "\t" << fstAll1[i] << "\t" << fst2All1[i] << "\t" << count1[i] << endl;
    }
    
    
    // All
    cout << endl << endl;
    for(size_t i=0; i<numCompare; i++){
        cout << setprecision(10) << NAllAll[i] << "\t" << N2AllAll[i] << "\t" << DAllAll[i] << "\t" << fstAllAll[i] << "\t" << fst2AllAll[i] << "\t" << countAll[i] << endl;
    }
    
    
    // Afr
    cout << endl << endl;
    for(size_t i=0; i<numCompare; i++){
        cout << setprecision(10) << NAllAfr[i] << "\t" << N2AllAfr[i] << "\t" << DAllAfr[i] << "\t" << fstAllAfr[i] << "\t" << fst2AllAfr[i] << "\t" << countAfr[i] << endl;
    }
    
    
    cout << "DONE!" << endl;
    return 0;
    
}
