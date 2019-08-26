# README #
[ ![Build Status](https://travis-ci.org/soedinglab/WIsH.svg?branch=master)](https://travis-ci.org/soedinglab/WIsH)

### LICENSE ###
WIsH is licensed under the General Public License (see the LICENSE file). Copyright Clovis Galiez (clovis.galiez@mpibpc.mpg.de).

<p align="center"><img src="https://raw.githubusercontent.com/soedinglab/WIsH/master/WIsHAlpha0Bg.png" height="256" /></p>


### What is this repository for? ###

* WIsH can identify bacterial hosts from metagenomic data, keeping good accuracy even on smaller contigs.
* version 1.0
* [Availabity](git@github.com:soedinglab/WIsH.git)
 
### How do I get set up? ###

#### Installation: ####

##### Docker: #####
Build the docker container:
```
#!bash
cd /path/to/repository/of/WiSH
docker build -t wish .
```

To run WIsH from the container:
```
docker run -v /some/host/folder:/data wish <some WiSH commands>
```

##### Linux: #####
```
#!bash
git clone https://github.com/soedinglab/WIsH.git
cd WIsH
cmake .
make
```

##### MacOS #####
Compiling with OpenMP support on MacOS requires a recent gcc compiler. You can get gcc from [homebrew](https://brew.sh/).

```
#!bash
brew install gcc@6
export CC=gcc-6
export CXX=g++-6

git clone https://github.com/soedinglab/WIsH.git
cd WIsH
cmake .
make
```


#### Dependencies ####
If you want to enjoy the parallelization of WIsH, you should have the OpenMP library installed. WIsH uses C++11. The model construction and the interaction prediction are both parallelized. In both cases, it spreads one bacterial genome/model per thread as soon as you give WIsH the number N of threads to use (with the parameter "-t N").

#### Database configuration ####

You need two different directory containing only sequence data in FASTA format. One should contain your potential host genomes (one becteria per FATSA file), the other should contain the viral contigs/genomes (one virus per FATSA file).

#### Usage example ####
To run a prediction, you should proceed in two steps:

1 - Create the models from the bacterial genomes you stored in FASTA format in prokaryoteGenomesDir:
```
#!bash
mkdir modelDir
./WIsH -c build -g prokaryoteGenomesDir -m modelDir
```
This will create a model in modelDir for every bectrial genome.

2 - Run the prediction on the viral sequences you stored in FASTA format in phageContigsDir:

```
#!bash
mkdir outputResultDir
./WIsH -c predict -g phageContigsDir -m modelDir -r outputResultDir -b
```
This will output a file *llikelihood.matrix* containing a matrix of log-likelihood (rows are a bacteria, and columns are a viral contigs), and a "summary" file *prediction.list* containing for every viral sequence the host corresponding to highest log-likelihood (-b option).

The files can be further analyzed with any text editor or with R:

```
#!R

ll = read.table("outputResultDir/llikelihood.matrix")
predictions = read.table("outputResultDir/prediction.list")

# Show the number of viral contigs targeting every potential hosts:
table(predictions$V2)

# Show the histogram of the log-likelihoods for the best predictions:
hist(predictions$V3)

```

3 - Advanced analysis

Finding the lowest common ancestor of the top WIsH predictions provides a more robust estimate of the host taxonomy. First, you can call WIsH to generate a matrix of p-values for each of the predictions:
```
./WIsH -c predict -g phageContigsDir -m modelDir -r outputResultDir -b -p -n KeggGaussianFits.tsv
```

Then in R, get the 5 top best predictions:
```
library(dplyr)
library(tidyr)

pvals = read.table("outputResultDir/pvalues.matrix")
ll = read.table("outputResultDir/llikelihood.matrix")

# Mask the predictions having a bad p-value:
ll[pvals > 0.05] <- NA

# Set the model name as a column:
ll$model = rownames(ll)

# Extract the best 5 predictions:
topPredictions = as.matrix(ll %>% gather(phage,likelihood,-model) %>% group_by(phage) %>% top_n(5))

```

You can now use tools such as (https://github.com/pmenzel/taxonomy-tools) to perform your LCA analysis.

#### Cut-offs ####

The cut-off are highly dependent on your use-case.
The likelihood value is dependent on the model and the best is to use null dataset to get the parameters of a null-model in order to get p-values [see next section for getting parameters for custom models](#getting-the-null-paramters-for-new-bacterial-models). Then you can refer to the benchmark in the supplementary materials: Figure 1, page 5 gives the p-value cut-offs and their subsequently associated accuracy and recall values.

In case you only want to get an idea of the usual log-likelihood values, you can check the provided parameter file (computed for Kegg Genomes models) available [here](KeggGaussianFits.tsv). The first column is the mean the second is standard deviation. As to set a cut-off, *the more negative, the better the model fits the sequence*. 


#### Getting the null paramters for new bacterial models ####
If you want to get p-values for your predictions, you need to know the null parameters for a new bacterial model. To get them, for each bacterial model, you must run the predictions on a large set of phage genomes that are known *not* to infect your bacterial model (let's call it the null set of phages) and use the prediction likelihood to fit the null-model parameters. We achieved that in our benchmark by simply removing for every bacterial model B all the phages that were known to infect the same genus as B, then running the prediction for every genus G:
```
./WIsH -c predict -g setOfPhageContigsNotInfectingGenus_G -m modelsOfGenus_G -r outputNullModelResultDir -b
```

And then for each genus G, computing the parameters for the associated null-models using the script computeNullParameters.R.


However, a simpler but less accurate method can be used to get all the parameters at once or also when lack of the taxonomy of your bacterial models. In that case, you have to collect a set of phage contigs/genomes such that for every bacterial model B_i the size of the set of phages P_i that infects B_i can be neglected compared to the set of phages that does *not* infect B_i. Then you can get the parameters for each model at once by running the prediction:

```
./WIsH -c predict -g superDiverseSetOfPhageContigs -m newModelsDir -r outputNullModelResultDir -b
```
where newModelsDir contains all your bacterial model that are missing their null-model parameters. Then, you can use the script computeNullParameters.R that takes the predictions and create a file containing the null parameters for every bacterial model in newModelsDir.

Afterwards, to get the p-values while predicting interactions, please specify the options "-b -n nullParameters.tsv" in a prediction call.

#### Tricks ####

When building models, **only one genome should be stored per fasta file**. If you have a big fasta file containing all your genomes, you can split them using the following commands. Supposing the the header is in the format ">recordId|description|...", go to your genome directory and call:
```
mkdir -p splittedGenomes
awk '/^>/ {if(x>0) {close(outname); x=0} match($0, ">([^| ]*)", record);outname=sprintf("splittedGenomes/%s.fa",record[1]); if (x>0) {print >> outname} else {print > outname;} x++; next;} {if(x>0) print >> outname;}' *.fna
```

The directory *splittedGenomes* now contains one genome per fasta file (with the recordId as file names), and can be used to train your WIsH models:
```
./WIsH -c build -g splittedGenomes -m modelDir
```


### Getting the benchmark data ###
The accession numbers of the sequences used in the benchmark can be found in the benchmark directory of the project.

You can use the NCBI API to get the phage sequences, for instance:
NC_000896 through (https://www.ncbi.nlm.nih.gov/search/all/?term=NC_000896)

To get the bacterial genomes, you can use the KEGG API, for instance
aac through (https://www.genome.jp/kegg-bin/show_organism?org=aac)






### Troubleshooting - Bug reports ###

* Please open a github issue or contact clovis.galiez@mpibpc.mpg.de
