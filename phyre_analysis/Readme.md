# Phylogenetic Representativeness v.1.1

This repository contains partially rewritten and optimized code from
"Phylogenetic representativeness: a new method for evaluating taxon sampling in
evolutionary studies" by Ronald R. Ferrucci, Federico Plazzi, and Marco Passamonti.

[Full text] (http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-209)

Original version is available at [sourceforge](http://sourceforge.net/projects/phyrephylogenet/).


This version is capable of processing larger master list files. 

## Running software

Options reference is contained at [document] (http://www.mozoolab.net/downloads/manual.pdf)

## Running tests

Tests can be run from console

```
python PhyRe.py tests/coleoids_sample.txt tests/coleoids_list.txt 5 15 -m y -p 100
```

or

```
sh generate_test_output.sh
```


or

```
from PhyRe import phy_re_analysis

args = {'d2': 15, 'popfile': 'tests/coleoids_list.txt', 'samplefile': 'tests/coleoids_sample.txt', 'd1': 5}
options = {'c': None, 'b': None, 'm': 'y', 'l': None, 'o': None, 'p': 100, 's': True, 'parallel': False}

(phyre_results, funnel_plot) = phy_re_analysis(options,args)


args = {'d2': 100, 'popfile': 'tests/carnivores_list.txt', 'samplefile': 'tests/carnivores_sample.txt', 'd1': 40}
options = {'c': None, 'b': None, 'm': 'y', 'l': None, 'o': None, 'p': 100, 's': True, 'parallel': False}

(phyre_results, funnel_plot) = phy_re_analysis(options,args)


args = {'d2': 100, 'popfile': 'tests/termites_list.txt', 'samplefile': 'tests/termites_sample.txt', 'd1': 40}
options = {'c': None, 'b': None, 'm': 'y', 'l': None, 'o': None, 'p': 100, 's': True, 'parallel': False}

(phyre_results, funnel_plot) = phy_re_analysis(options,args)
```
