#include "mm.h"

mm::mm(std::string modelFile, int verb)
{
    verbosity = verb;
    modelName = modelFile.substr(0, modelFile.find_last_of("."));
    unsigned int slashPos;
    if ((slashPos = modelName.find_last_of('/'))!= modelName.size())
        modelName = modelName.substr(slashPos+1,modelName.size());
        
    read(modelFile);
}

mm::mm(unsigned int k,double a, int verb)
{
    verbosity = verb;
    order = k;
    alpha = a;
    initArrays();
}

void mm::printParameters()
{
    std::cout<<"Parameters for model "<<modelName<<std::endl;
    std::cout<<"Order: "<<order<<std::endl;
    std::cout<<"Alpha: "<<alpha<<std::endl;
    std::cout<<"Probability vector size: "<<p.size()<<std::endl;
    std::cout<<"lowerOrderCounts vector size: "<<lowerOrderCounts.size()<<std::endl;
    std::cout<<"Pseudo counts vector size: "<<pc.size()<<std::endl;
}

std::string mm::getName() {
    return modelName;
}

void mm::initArrays()
{
    lowerOrderCounts = std::vector<unsigned int> (1<<(2*(order)),0);
    orderCounts = std::vector<unsigned int>(1<<(2*(order+1)),0);
    p = std::vector<double>(1<<(2*(order+1)),0);
    pc = std::vector<double>(4,0.25);
}

int mm::trainOn(std::string genomeFile)
{
    modelName = genomeFile.substr(0, genomeFile.find_last_of("."));
    unsigned int slashPos;
    if ((slashPos = modelName.find_last_of('/'))!= modelName.size())
        modelName = modelName.substr(slashPos+1,modelName.size());
        
    
    std::vector<std::string> genomeChunks = readGenome(genomeFile);

    if (!genomeChunks.size())
        return -1;
        
    for (size_t i = 0 ; i < genomeChunks.size() ; i++)
    {
        countKmers(genomeChunks[i]);
    }
    
    computeModelProb();
    
    return 0;
}


std::string mm::mapToAlphabet(std::string seq)
{
    bool seenUnknown = false;
    std::string res;
    for (size_t i = 0 ; i < seq.size() ; i++)
    {
        switch (seq[i]) {
            case 'A':
                res.push_back(0);
                break;
            case 'T':
                res.push_back(1);
                break;
            case 'C':
                res.push_back(2);
                break;
            case 'G':
                res.push_back(3);
                break;
            default:
                seenUnknown=true;
                break;
        }
        
    }
    
    /*if (seenUnknown)
        std::cout<< "Warning: there are letters in genome not in {A,T,C,G}."<<std::endl;
    */
    return res;
}

size_t mm::hashKmer(std::string::iterator kmer, unsigned int k = -1)
{
    size_t val = 0;
    
    if (k==-1)
        k=order+1;
    
    for (size_t pos = 0; pos < k ; pos++)
    {
        val += (1<<(2*(k-1-pos))) * kmer[pos];
    }
    return val;
}


void mm::countKmers(std::string genome)
{
    size_t start = order;
    for (size_t pos = start; pos<genome.size();pos++)
    {
        size_t kmer = hashKmer(genome.begin() + pos-start);

        orderCounts[kmer]++;
        lowerOrderCounts[head(kmer)]++;
    }
    
    // count the last k-1 mer
    if (genome.size() >= order)
        lowerOrderCounts[hashKmer(genome.end() - order,order)]++;
}

size_t mm::head(size_t k)
{
    return k>>2;
}


size_t mm::lastNucl(size_t k)
{
    return k&3;
}



void mm::computeModelProb()
{
    for (size_t k=0; k<p.size();k++)
    {
        p[k] = log((double)(orderCounts[k] + alpha*pc[lastNucl(k)]) / (lowerOrderCounts[head(k)] + alpha));
        if (p[k] >0 && verbosity > 1)
            std::cout <<"OUPS ! Positive log-probability : "<< p[k] <<std::endl;
    }
}

mm::~mm()
{
}


int mm::write(std::string modelDir)
{
    std::string modelFile(modelDir);
    modelFile += "/" + modelName + ".mm";
    
    
    std::ofstream fout(modelFile.c_str(), std::ios::out|std::ios::binary);
    if (!fout.good())
        return -1;
    
    fout.write((char*)&order,sizeof(order));
    fout.write((char*)&alpha,sizeof(alpha));

    fout.write(reinterpret_cast<char*>(&p[0]),p.size() * sizeof(p[0]));
    fout.close();
    return 0;
}


int mm::read(std::string modelFile)
{
    std::ifstream fin(modelFile.c_str(), std::ios::in|std::ios::binary);
    if (!fin.good())
        return -1;
        
    fin.seekg(0,fin.end);
    size_t sizep = static_cast<unsigned long>(fin.tellg()) - sizeof(alpha) - sizeof(order);
    
    
    p = std::vector<double>(sizep / sizeof(double));
    
    fin.seekg(0,std::ios::beg);
    fin.read(reinterpret_cast<char*>(&order),sizeof(order));
    fin.read(reinterpret_cast<char*>(&alpha),sizeof(alpha));
    
    
    
    fin.read(reinterpret_cast<char*>(p.data()),p.size() * sizeof(double));
    fin.close();
    
    return 0;
}

std::vector<std::string> mm::readGenome(std::string genomeFile)
{
    std::ifstream fin(genomeFile.c_str(), std::ios::in);
    
    if (!fin.good())
        return std::vector<std::string>();
    
    std::vector<std::string> genomeChunks;
    std::string line;
    
    while(std::getline(fin,line))
    {
        if(!line.empty())
        {
            if (line[0]=='>')
            {
                genomeChunks.push_back(std::string());
            } else {
                if (!genomeChunks.size())
                {
                    std::cout<<"Corrupted Fasta file: "<<genomeFile<<std::endl;
                    genomeChunks.push_back(std::string());
                }
                genomeChunks.back() += mapToAlphabet(line);
            }
        }
    }
    fin.close();
    return genomeChunks;
}

double mm::evaluate(std::vector<std::string> genomeChunks)
{
    double ll = 0;
    size_t length = 0;
    for (size_t i = 0 ; i < genomeChunks.size() ; i++)
    {
        size_t start = order;
        for (size_t pos = start; pos<genomeChunks[i].size() ; pos++)
        {
            length++;
            ll += p[hashKmer(genomeChunks[i].begin() + pos - start)]; // TODO may be optimized by precomputing the hashing
        }
    }
    
    
    return ll/length;
}
