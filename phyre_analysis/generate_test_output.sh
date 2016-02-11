#!/bin/sh
python PhyRe.py data/carnivores_sample.txt data/carnivores_list.txt 40 100 -m y -p 100
python PhyRe.py data/bivalves_sample.txt data/bivalves_list.txt 5 15 -m y -p 100
python PhyRe.py data/coleoids_sample.txt data/coleoids_list.txt 20 40 -m y -p 100
python PhyRe.py data/termites_sample.txt data/termites_list.txt 20 60 -m y -p 100
