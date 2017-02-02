#ifndef MM_H
#define MM_H


#include <vector>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>


class mm
{
    int verbosity;
    
    std::string modelName;
    unsigned int order;
    
    double alpha;
    
    std::vector<unsigned int> lowerOrderCounts,orderCounts;
    std::vector<double> pc; // pseudo count probabilities
    std::vector<double> p;
    
    void initArrays();
    int read(std::string modelFile);
    
    
    void computeModelProb();
    
    
    size_t hashKmer(std::string::iterator kmer, unsigned int k);
    static std::string mapToAlphabet(std::string seq);
    void countKmers(std::string genome);
    size_t lastNucl(size_t k);
    
    size_t head(size_t k);
    
    
public:
    mm(std::string modelFile, int verb);
    mm(unsigned int k,double a, int verb);
    ~mm();

    void printParameters();
    std::string getName();
    
    int write(std::string modelDir);
    int trainOn(std::string genomeFile);
    static std::vector<std::string> readGenome(std::string genomeFile);
    double evaluate(std::vector<std::string>);
};

#endif // MM_H
