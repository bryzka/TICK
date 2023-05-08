# TICK
Repository for project for the "Architecture of large projects in bioinformatics".
We would like to create "TICK" - a tool for analysis sequences and structure of potential pathogenes and parasites of ticks. Our software will be using BLAST and DIAMOND to align sequences to records from our local database. Our aim is to create a user-friendly interface and provide some statistics from analysis. We also want to create models presenting structure of provided sequences.

# Demo
Demo is supposed to show that what we want to achieve is doable. Right now 
it's highly simplified version of our pipeline, currently using only NCBI BLAST, 
and doesn't have the possibility of using local blast instance and local database. 

## Work plan
 -[x] ~~Create demo~~
- [ ] Possibility of using local blast and database
  - [ ] Create tick's pathogens list
  - [ ] Create database out of tick's pathogens proteomes
- [ ] Elongation of the pipeline 
  - [ ] Decide what information we need about tick's pathogen sequences
    - Structure?
    - Amino acid composition?
    - Signal proteins?
    - Transmembrane helices?
    - ...?
  - [ ] Choose the tools needed for getting this information
  - [ ] Integrate them into the pipeline



# Few remarks on config and usage
Config file consists of two parts, required and optional parameters. 
Required parameters are path to blast (simply WWW, or empty string 
if NCBI blast is to be used), and path to file containing list of ticks pathogens. 
Optional are paths to tools that are yet to be added to function, or the ones already added, but optional. 

Usage of the demo is as follows( full version's basic usage should be similar):
```python
usage: demo.py [-h] [--config_file CONFIG_FILE] [--email EMAIL] input_file

positional arguments:
  input_file            Fasta file containing sequences that we want to check.

options:
  -h, --help            show this help message and exit
  --config_file CONFIG_FILE, -c CONFIG_FILE, -config CONFIG_FILE
                        Path to config file. DEFAULT: ./config
  --email EMAIL         email for blasting purposes


```


# Members
Julia Byrska, Stanisław Janik
