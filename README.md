### Metozoa Mitochondrial DNA Database ###

* Suitable for phylogenetic analysis of metazoa
* Phylogenetic representativeness of the sample confirmed


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