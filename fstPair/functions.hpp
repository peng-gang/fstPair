//
//  functions.hpp
//  fstCalculation
//
//  Created by Gang Peng on 6/8/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>
#include <map>
#include <string>

/*
 * read sample size file
 * format
 * POPULATION    DESCRIPTION    GENOMES
 * afr    African/African-American    21042
 * ami    Amish    450
 * amr    Latino/Admixed American    6835
 * asj    Ashkenazi Jewish    1662
 * eas    East Asian    1567
 * fin    Finnish    5244
 * nfe    Non-Finnish European    32299
 * sas    South Asian    1526
 * oth    Other (population not assigned)    1077
 * return a map rlt
 * rlt["afr"] = 21042
 * rlt["ami"] = 450
 * ...
 */
std::map<std::string, int>  getSampleSize(std::string filename);

/*
 calculaiton fst with Hudson method (biased estimation)
 */
double fstHudson(double p1, double p2);
/*
 calculaiton fst with Hudson method (unbiased estimation)
 */
double fstHudson(double p1, double p2, int s1, int s2);

double fstNHudson(double p1, double p2);

double fstNHudson(double p1, double p2, int s1, int s2);

double fstDHudson(double p1, double p2);


#endif /* functions_hpp */
