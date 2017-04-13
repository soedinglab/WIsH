# README #
[ ![Build Status](https://travis-ci.org/soedinglab/WIsH.svg?branch=master)](https://travis-ci.org/soedinglab/WIsH)

### LICENSE ###
WIsH is licensed under the General Public License (see the LICENSE file). Copyright Clovis Galiez (clovis.galiez@mpibpc.mpg.de).

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
If you want to enjoy the parallelization of WIsH, you should have the OpenMP library installed. WIsH uses C++11.

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



### Troubleshooting - Bug reports ###

* Please open a github issue or contact clovis.galiez@mpibpc.mpg.de
