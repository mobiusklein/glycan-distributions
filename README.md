# glycan distributions

A Python library, webapp, and command-line script for quantifying differences in site-specific protein glycosylation

## Getting Started
To use this repo, clone it into a directory of your choosing and install the matplotlib, numpy, flask, and glycopeptidepy Python libraries.
```
# clone this repo
git init
git clone https://github.com/jturcino/glycan-distributions.git

# install libraries
pip install matplotlib
pip install numpy
pip install flask
pip install glycopeptidepy
```

## The Program
Whether run as a webapp or from the command line, the program takes a set of inputs for "expected" and "experimental" glycoproteins, models their site-specific glycan populations as n-dimensional probability vectors, and calculates similarity metrics for aligning glycosylation sites. 

The inputs, similarity metrics, and outputs are described in more detail below. If you'd like to skip to running the program, you can click [here](https://github.com/jturcino/glycan-distributions/README.md#running-the-webapp-on-the-scc) for the webapp (SCC access required) or [here](https://github.com/jturcino/glycan-distributions/README.md#running-the-script) for the command line.

### Inputs
The program compares glycosylation between "expected" and "experimental" glycoproteins. This could be a well characterized protein and a novel variant, or a wild type protein and mutant, etc.

For each protein, two inputs are required:
1. **[GlycReSoft](https://github.com/mobiusklein/glycresoft) CSV** created with the `glycopeptide-identification` step. It is advisable to provide more than one CSV per protein.
2. **[UniProt](http://www.uniprot.org/) accession ID**, which provides the program with the protein's sequence. If the protein's sequence is not available on UniProt, a fasta file can be submitted using the webapp.

Optional inputs:
* **Dataset names** labels expected and experimental datasets with something other than the accession ID
* **Replicate cutoff** (defaults to 1) removes glycans observed in fewer CSVs than the number given
* **Score cutoff** (defaults to 30) removes individual glycopeptide observations whose MS2 score is less than the number given
* **Sequon length** (defaults to 4) determines  the number of amino acids that must align for glycosylation sites to match between the expected and experimental proteins.

### Similarity metrics
Once each protein's glycan populations are modeled as n-dimensional probability vectors, their glycosylation sites are aligned based on amino acid sequence. The similarity between the vectors associated with that site is the quantified the following metrics
* **[Cosine similarity](https://en.wikipedia.org/wiki/Cosine_similarity)** calculates the cosine of the angle between the two vectors. A smaller angle means more similar vectors, which results in a cosine close to one.
* **[Normalized Kullback-Leibler divergence](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence)** measures how well the expected distribution describes the experimental distribution. A better description means more similar distributions, which results in a divergence value close to zero.

Note: If one or both of the aligned glycosylation sites do not have any coverage, the site will be ignored and similarity metrics will not be calculated.

### Outputs
Three sets of outputs are produced. When using the web app, all outputs are displayed in the webpage and can be downloaded at will. When using the script, outputs are saved to the output directory automatically.
* **Similarity CSV** contains information on aligned glycoyslation sites, cosine similarity values, and normalized Kullback-Leibler divergences.
* **Probability plots** show the glycan probability distributions of each protein's glycosylation sites as a bar chart.
* **Comparison plots** show the glycan probabiltiy distributions of aligned glycosylation sites side-by-side in a single bar chart.

## Running the Webapp on the [SCC](https://www.bu.edu/tech/support/research/computing-resources/scc/)
If working on the SCC, the webapp can easily be accessed in two steps:
1. **Launch** the app by executing the following command in the cloned directory: `./app.py`
2. **Open** the app in your web browser by vising the following address: `http://scc1.bu.edu:9898/`

This directs you to the home page, which explains how to use the app and from which you can access the submission form. If you prefer to go straight to the form, visit `http://scc1.bu.edu:9898/csv`.

## Running the Script
If you do not have access to the SCC, you can still run the program from the command line, provided the protein sequences in question have an entry in the UniProt database.

The protein sequence accession IDs, GlycReSoft CSVs, and path to an output directory are fed to the script via command-line arguments, as seen below. Addtional options can be viewed at any time using the script's `-h` flag.
```
./glycan_distributions.py -xc expected1.csv -xc expected2.csv -pc experiment1.csv -pc experiment2.csv -xa EXPECTID -pa EXPERIMENTID -o path/to/output/directory
```
The plots and CSV file produced by the script will be saved in the provided output directory.
