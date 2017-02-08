# README #


### What is this repository for? ###

* WIsH can identify bacterial hosts from metagenomic data, keeping good accuracy even on smaller contigs.
* 1.0
* [Availabity](https://bitbucket.org/ClovisG/wish)

### How do I get set up? ###

#### Installation: ####

```
#!bash
git clone git@bitbucket.org:ClovisG/wish.git
cd wish
cmake .
make
```


#### Dependencies ####
If you want to enjoy the parallelization of WIsH, you should have the OpenMP library installed. WIsH uses C++11.

#### Database configuration ####

You need two different directory containing only sequence data in FASTA format. One should contain your potential host genomes (one becteria per FATSA file), the other should contain the viral contigs/genomes (one virus per FATSA file).

#### Usage example ####
To run a prediction, you should proceed in two steps:

1 - Create the models from the bacterial genomes stored in prokaryoteGenomesDir:
```
#!bash
mkdir modelDir
./WIsH -c build -g prokaryoteGenomesDir -m modelDir
```
This will create a model in modelDir for every bectrial genome.

2 - Run the prediction on your viral sequences stored in phageContigsDir:

```
#!bash
mkdir outputResultDir
./WIsH -c predict -g phageContigsDir -m modelDir -r outputResultDir -b
```
This will output a file containing a matrix of log-likelihood, and a "summary" file containing for every viral sequence the host corresponding to highest log-likelihood (-b option).



### Troubleshooting - Bug reports ###

* Please contact clovis.galiez@mpibpc.mpg.de