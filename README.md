# popgen_utils
Utilities for Lauren Sugden's lab to run simulations, statistics on the simulations, and analyze the output

### Installation
If you will be working on the utils repo, clone first and install in develop mode.
First, move to the directory in which you have code. The following line will create
a new directory in your code directory called popgen_utils, and everything
from the git repository will be loaded within.
```
git clone git@github.com:lasugden/popgen_utils popgen_utils
python setup.py develop --user
```

If instead you will be using the utilities, install with:
```
pip install --upgrade git+git://github.com/lasugden/popgen_utils.git
```

Next, set your main data directory and output directory. To do so, open python and type
```
import popgen_utils.config
popgen_utils.config.reconfigure()
```
And enter the path to your data directory and output file directory.
From now on, everything will be saved to one of those directories.

### Steps

1. Create inputs for a simulator ```write_slim_input.py``` (1 min/sim)
1. Run the simulator at the command line ???
1. Convert the slim output ```*_ms.txt``` into selscan input ```_map.txt```, ```*_p{0-9}.hap```


### SNP-based statistics:
_selscan_
* xp-ehh
* ihs
* ihh12

_vcftools_
* fst

_a guy who made software_
* sds

_a different guy who made software_
* isafe

### Region-based statistics:
Lauren wrote code to compute statistics
