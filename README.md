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

#### Getting the null paramters for new bacterial models ####
If you want to get p-values for your predictions, you need to know the null parameters for a new bacterial model. To get them, you must run the predictions on a large set of phage genomes that are know *not* to infect your bacterial model (let's call it the null set of phages) and use the prediction likelihood to fit the null-model parameters. You can use the script computeNullParameters.R that takes the predictions on this null set of phage and create a file containing the null parameters for every bacterial model.
To get the p-values while predicting interactions, please specify the options "-b -n nullParameters.tsv" in a prediction call.


### Troubleshooting - Bug reports ###

* Please open a github issue or contact clovis.galiez@mpibpc.mpg.de
