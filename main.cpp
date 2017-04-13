#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <iostream>
#include <map>
#include <cmath>
#include <float.h>
#include <sstream>
#include "mm.h"
#include "main.h"

#ifdef OPENMP
    #include <omp.h>
#endif

#define VERBOSITY 3

#define BUILD_COMMAND "build"
#define PREDICT_COMMAND "predict"


void die(std::string text,std::string complement=std::string())
{
    std::cout<< "ERROR:" << text << complement << "." <<std::endl;
    exit(-1);
    
}

void build(std::string genomeDir, std::string modelDir, unsigned int order, double alpha, unsigned int threads = 1)
{
    DIR* dir;
    struct dirent *genomeFile;
    std::vector<std::string> genomeFiles;
    
    #ifdef OPENMP
        omp_set_num_threads(threads);
    #endif
    
    if((dir = opendir(modelDir.c_str())) == NULL)
        die("Cannot open model directory ", modelDir);
        
    if((dir = opendir(genomeDir.c_str())) == NULL)
        die("Cannot open genome directory ", genomeDir);

    while ((genomeFile = readdir (dir)) != NULL) {
        if (std::string(genomeFile->d_name) != "." && std::string(genomeFile->d_name) != "..")
            genomeFiles.push_back(genomeDir + "/" + std::string(genomeFile->d_name));
    }
    closedir (dir);
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < genomeFiles.size() ; i++)
    {
        mm model(order,alpha,VERBOSITY);
        if (VERBOSITY > 2)
            std::cout << "Processing "<<genomeFiles[i]<<std::endl;
        model.trainOn(genomeFiles[i]);
        model.printParameters();
        if (model.write(modelDir) < 0)
            die("Cannot write to model directory ",modelDir);
    }
    

}



double mean(std::vector<double> v)
{
    double mean = 0.0;
    
    for (size_t i = 0; i < v.size();i++)
    {
        mean += v[i];
    }
    return mean/v.size();
}


double sd(std::vector<double> v, double mean)
{
    double sd = 0.0;
    
    for (size_t i = 0; i < v.size();i++)
    {
        double x = (v[i]-mean);
        sd += x*x;
    }
    return (double)(v.size() -1 ) * sqrt(sd/v.size()) / v.size();
}


double getPval(std::pair<double,double> param, double ll)
{
    return 0.5 - 0.5*erf((ll-param.first)/(sqrt(2.0)*param.second));
}


void predict(std::string genomeDir, std::string modelDir,std::string resultDir, unsigned int threads = 1, bool writeLLMatrix=true,bool writeBestPred = false,std::string negFitsFile = std::string(),bool zScores = false)
{
    
    DIR* dir;
    struct dirent *genomeFile,*modelFile;
    std::vector<std::string> genomeFiles,modelFiles, genomeNames;
    std::vector< std::vector<double> > ll;
    
    
    #ifdef OPENMP
        omp_set_num_threads(threads);
    #endif
    
    
    // List files in genome and model dir
    if((dir = opendir(resultDir.c_str())) == NULL)
        die("Cannot open result directory ", resultDir);
    closedir(dir);
    
    if((dir = opendir(genomeDir.c_str())) == NULL)
        die("Cannot open genome directory ", genomeDir);
        
    
    while ((genomeFile = readdir (dir)) != NULL) {
        if (std::string(genomeFile->d_name) != "." && std::string(genomeFile->d_name) != "..")
        {
            genomeFiles.push_back(std::string(genomeFile->d_name));
            genomeNames.push_back(genomeFiles.back().substr(0, genomeFiles.back().find_last_of(".")));
        }
    }
    closedir (dir);
    
    if((dir = opendir(modelDir.c_str())) == NULL)
        die("Cannot open model directory ", modelDir);
        
    while ((modelFile = readdir (dir)) != NULL) {
        if (std::string(modelFile->d_name) != "." && std::string(modelFile->d_name) != "..")
        {
            modelFiles.push_back(std::string(modelFile->d_name));
            ll.push_back(std::vector<double>());
        }
    }
    closedir (dir);
    
    
    std::vector<std::vector<std::string> > bactGenomes;
    
    for (size_t j = 0 ; j < genomeFiles.size() ; j++) {
        bactGenomes.push_back(mm::readGenome(genomeDir + "/" + genomeFiles[j]));
    }
    
    std::vector<std::string> modelNames(modelFiles.size(),std::string());
    // Compute log-likelihoods
    size_t progressCount = 0;
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < modelFiles.size() ; i++)
    {
        mm model(modelDir + "/" + modelFiles[i],VERBOSITY);
    
        if (VERBOSITY>2)
        {
            progressCount++;
            std::cout<<"Processing "<<model.getName()<<" (approx. "<<progressCount<<"/"<<modelFiles.size()<<")"<<std::endl; //model.printParameters();
        }
        
        modelNames[i] = model.getName();
        
        for (size_t j = 0 ; j < genomeFiles.size() ; j++) {
            ll[i].push_back(model.evaluate(bactGenomes[j]));
        }
        
    }
    
    if (writeBestPred)
    {
        double maxLL;
        size_t host = 0;
        std::map<std::string,std::pair<double, double> > negFits;
        
        // If we need to compute a p-value, read the fits from file
        if (negFitsFile != std::string())
        {
            std::ifstream fin(negFitsFile.c_str(), std::ios::in);
    
            if (!fin.good())
                die("Cannot open negative fits file ",negFitsFile);
            
            
            std::string line;
            
            while(std::getline(fin,line))
            {
                if(!line.empty())
                {
                    int secondField = line.find("\t");
                    if (secondField+1<line.size())
                    {
                        int thirdField = line.find("\t",secondField + 1);
                        if (thirdField+1<line.size())
                        {
                            std::string bactName = line.substr(0,secondField);
                            double mu = std::stod(line.substr(secondField+1,thirdField));
                            double s = std::stod(line.substr(thirdField+1,line.size()));
                            negFits[bactName] = std::make_pair(mu,s);
                        } else {
                            std::cout << "Warning: "<< negFitsFile <<" is missformatted. Should be [bactName]\\t[mu]\\t[standard deviation]"<<std::endl;
                        }
                    } else {
                        std::cout << "Warning: "<< negFitsFile <<" is missformatted. Should be [bactName]\\t[mu]\\t[standard deviation]"<<std::endl;
                    }
                    
                }
            }
            fin.close();
        }
        
        
        
        std::ofstream fout((resultDir + "/prediction.list").c_str(), std::ios::out);
        
        if (!fout.good())
            die("Cannot open ",resultDir);
        fout << "\"Phage\"\t\"Best hit among provided hosts\"\t\"LogLikelihood\"\t\"p-value if null parameters provided\"\n";
        
        for (size_t j = 0 ; j < genomeNames.size() ; j++)
        {
            maxLL = -DBL_MAX;
            for (size_t i = 0; i < modelNames.size() ; i++)
            {
                if (ll[i][j] > maxLL)
                {
                    host = i;
                    maxLL = ll[i][j];
                }
            }
                
            fout << genomeNames[j] <<'\t'<<modelNames[host]<<'\t'<<ll[host][j];
            
            
            if ( negFits.find(modelNames[host]) != negFits.end() )
            {
                fout << '\t' << getPval(negFits[modelNames[host]],ll[host][j]);
            } else {
                fout << "\tNA";
            }
            
            fout << std::endl;
        }
        fout.close();
    }
    
    
    // Normalize the log-likelihood to z-scores
    if(zScores)
    {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < modelFiles.size() ; i++)
        {
            double mu = mean(ll[i]);
            double s = sd(ll[i],mu);
            
            for (size_t j = 0 ; j < genomeFiles.size() ; j++) {
                ll[i][j] -= mu;
                ll[i][j] /= s;
            }
            
        }
    }    
    
    
    // Output results according to user choice
    if (writeLLMatrix)  {
        std::ofstream fout((resultDir + "/llikelihood.matrix").c_str(), std::ios::out);
        
        if (!fout.good())
            die("Cannot open ",resultDir);

        if (genomeNames.size())
            fout<<genomeNames[0];
        for (size_t j = 1 ; j < genomeNames.size() ; j++)
        {
            fout << '\t'<< genomeNames[j];
        }
        fout << std::endl;
        
        for (size_t i = 0; i < modelNames.size() ; i++)
        {
            fout << modelNames[i];
            for (size_t j = 0 ; j < genomeNames.size() ; j++)
            {
                fout << '\t' << ll[i][j];
            }
            fout << std::endl;
            
        }
        fout.close();
    }
    
        

}



int main(int argc, char **argv)
{
    std:: string genomeDir,modelDir,resultDir,command,negFitFile;
    
    
    unsigned int order = 8;
    unsigned int threads = 1;
    double alpha = 16.0;
    bool bestPred = false;
    bool zScores = false;
    bool printHelp = false;
    
    int option;
    while ((option = getopt (argc, argv, "g:m:r:c:k:bn:zha:t:")) != -1)
    {
        switch(option)
        {
            case 'g':
                genomeDir = optarg;
                break;
            case 'm':
                modelDir = optarg;
                break;
            case 'r':
                resultDir = optarg;
                break;
            case 'c':
                command = optarg;
                break;
            case 'n':
                negFitFile = optarg;
                break;
            case 'b':
                bestPred = true;
                break;
            case 'k':
                order = atoi(optarg);
                break;
            case 't':
                threads = atoi(optarg);
                break;
            case 'a':
                alpha = std::stod(optarg);
                break;
            case 'z':
                zScores = true;
                break;
            case 'h':
                printHelp = true;
                break;
            default:
                die("Bad option");
                break;
        }
    }
    
    std::string helpText;
        helpText = "WIsH (v" + std::to_string(WIsH_VERSION_MAJOR) + "." + std::to_string(WIsH_VERSION_MINOR) + ") is a tool for predicting bacterial hosts from phage (meta)genomic data.\n\
© Clovis Galiez (clovis.galiez@mpibpc.mpg.de)\n\n\
Usage :" + std::string(argv[0]) + " [options] \n\
Options:\n\
\t-c\tCommand to be executed (build or predict)\n\
\t-k\tOrder for building the Markov chain (default is " + std::to_string(order) +")\n\
\t-a\tPseudo-count parameter (default is " + std::to_string(alpha) +")\n\
\t-t\tNumber of threads to be used (default is " + std::to_string(threads) +")\n\n\
Path specifications:\n\
\t-g\tSpecifies the genome directory (read access)\n\
\t-m\tSpecifies the model directory (read/write access)\n\
\t-r\tSpecifies the result directory (write access)\n\n\
Score options:\n\
\t-b\tOutputs a file containing for each viral sequence the host with highest likelihood\n\
\t-z\tNormalize the matrix of log-likelihood as z-scores\n\
\t-n\tSpecifies the parameters for the distribution of negative values of each model\n\
\t\t\tFormat should be: modelName<Tab>mean<Tab>standardDeviation\n\n\
Example for building models:\n\n\
\t\tWIsH -c build -g prokaryoteGenomesDir -m modelDir\n\n\
Example for predicting hosts:\n\n\
\t\tWIsH -c predict -g virusGenomesDir -m modelDir -r outputResultDir\n\n";

    if (printHelp)
    {   
        std::cout<<helpText;
        #ifdef OPENMP
            std::cout<<"OpenMP supported."<<std::endl;
        #else
            std::cout<<"Compiled *without* OpenMP support."<<std::endl;
        #endif
        exit(0);
    }
        
    
    
    if (command == BUILD_COMMAND)
    {
        build(genomeDir,modelDir,order,alpha, threads);
    } else if (command == PREDICT_COMMAND)
    {
        predict(genomeDir,modelDir,resultDir, threads, true, bestPred, negFitFile,zScores);
    } else {
        die(std::string("Bad options.\n") + helpText);
    }
        
	return 0;
}
