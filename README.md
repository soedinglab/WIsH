# README #
[ ![Build Status](https://travis-ci.org/soedinglab/WIsH.svg?branch=master)](https://travis-ci.org/soedinglab/WIsH)

### LICENSE ###
WIsH is licensed under the General Public License (see the LICENSE file). Copyright Clovis Galiez (clovis.galiez@univ-grenoble-alpes.fr).

<p align="center"><img src="https://raw.githubusercontent.com/soedinglab/WIsH/master/WIsHAlpha0Bg.png" height="256" /></p>


### What is this repository for? ###

* WIsH can identify bacterial hosts from metagenomic data, keeping good accuracy even on smaller contigs.
* version 1.1
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

##### Singularity: #####
Singularity is an alterative to Docker, which is more suited for
HPC environments.

Build the Singularity container:
```
#!bash
cd /path/to/repository/of/WiSH
singularity build -f --nv --force $PWD/WIsH.sif $PWD/WIsH.def 
```

To run the WIsH singularity container interactively: 
```
singularity shell WIsH.sif
```

To run the WIsH singularity container:
```
## This expects the [options] from WIsH
singularity run WIsH.sif [options]
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
./WIsH -c predict -g phageContigsDir -m modelDir -r outputResultDir -b 1
```
This will output a file *llikelihood.matrix* containing a matrix of log-likelihood (rows are a bacteria, and columns are a viral contigs), and a "summary" file *prediction.list* containing for every viral sequence the host corresponding to highest log-likelihood (-b k option output the k best prediction by likelihood).

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

#### Getting the null paramters for new bacterial models ####
If you want to get p-values for your predictions, you need to know the null parameters for a new bacterial model. To get them, for each bacterial model, you must run the predictions on a large set of phage genomes that are known *not* to infect your bacterial model (let's call it the null set of phages) and use the prediction likelihood to fit the null-model parameters. We achieved that in our benchmark by simply removing for every bacterial model B all the phages that were known to infect the same genus as B, then running the prediction for every genus G:
```
./WIsH -c predict -g setOfPhageContigsNotInfectingGenus_G -m modelsOfGenus_G -r outputNullModelResultDir -b 1
```

And then for each genus G, computing the parameters for the associated null-models using the script computeNullParameters.R.


However, a simpler but less accurate method can be used to get all the parameters at once or also when lack of the taxonomy of your bacterial models. In that case, you have to collect a set of phage contigs/genomes such that for every bacterial model B_i the size of the set of phages P_i that infects B_i can be neglected compared to the set of phages that does *not* infect B_i. Then you can get the parameters for each model at once by running the prediction:

```
./WIsH -c predict -g superDiverseSetOfPhageContigs -m newModelsDir -r outputNullModelResultDir -b 1
```
where newModelsDir contains all your bacterial model that are missing their null-model parameters. Then, you can use the script computeNullParameters.R that takes the predictions and create a file containing the null parameters for every bacterial model in newModelsDir.

Afterwards, to get the p-values while predicting interactions, please specify the options "-b 1 -n nullParameters.tsv" in a prediction call.

### Troubleshooting - Bug reports ###

* Please open a github issue or contact clovis.galiez@mpibpc.mpg.de
