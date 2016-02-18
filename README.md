Metozoa Mitochondrial DNA Database
====================================

* Suitable for phylogenetic analysis of metazoa
* Phylogenetic representativeness of the sample confirmed

MeMitDB (Metazoan Mitchondrial DataBase) is an open database of curated mitochondrial genomes. 
This is a completely open source project, following the reproducible research 
principles; anyone can get copy or contribute.

Aim of the database is the collection of as much as possible different mitochondrial genomes in the same place, 
following standardized rules for annotations and retrieving. Therefore, publicly available sequences were downloaded 
from [GenBank](http://www.ncbi.nlm.nih.gov/), filtered, and parsed in order to eliminate redundancy.
Basic statistics of a given subset of sequences can be computed directly online, while for thorough analyses 
all sequences are easily available for download in FASTA format. Many filtering option are available in order 
to make selection quick and easy.

The pipeline of populating MetMitDB is entirely automated thanks to several Python scripts
(available in this repository), which is essential in order to efficiently update and keep 
the database alive. 

The date of the last update was 2016 February, 16.


### Running the workbook ###

All data collection was done with workbook.ipynb file, which can be run through ipython server. To run 
ipython server you can use script *start_ipython_notebook_server.sh* or just type 
```
ipython notebook --pylab=inline
``` 
to run it

### SQLite database ###

To get the database without taxonomical data, donwload it from [Yandex.Disk](https://yadi.sk/d/R4A-2OC5ogJJF)

### Running the database online ###

```
mkvirtualenv metmitdb
pip install -r requirements.txt
python metmitdb.py
```

### Example ###
[MetMitDb.info](http://metmitdb.info)