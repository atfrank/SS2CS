# SS2CS: A Secondary Structure to Chemical Shift Prediction Tool

# Installation:

```
conda create -n ss2cs python=3.7
source activate ss2cs
pip install scikit-learn==0.22.2.post1 
pip install pandas
```

# Usage:
```
usage: ss2cs.py [-h] -i INPUT -o OUTPUT -s SS2CS_PATH

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to read input CT secondary structure file
  -o OUTPUT, --output OUTPUT
                        path to save output chemical shift file
  -s SS2CS_PATH, --ss2cs_path SS2CS_PATH
                        path to SS2CS repo
```

```shell
python ss2cs.py -i test/1HWQ_0.ct -o test/1HWQ_0.csv -s ~/Documents/GitHub/SS2CS
```

## Output format
residue-name, residue-number, nucleus, predicted-shifts, extra

```shell
GUA 1 C1' 92.11827999999988 .
GUA 2 C1' 93.75492200000001 .
URA 3 C1' 94.07885000000013 .
GUA 4 C1' 93.26077999999994 .
CYT 5 C1' 94.225774 .
GUA 6 C1' 93.28166333333328 .
ADE 7 C1' 93.81845999999994 .
ADE 8 C1' 92.4785283333333 .
GUA 9 C1' 92.60993366666662 .
GUA 10 C1' 93.25432316666662 .
```

# COMMERCIAL USE LICENSE:

If you are interested in commercial licensing of these applications (clinical, operational, etc.) please contact the University of Michigan Office of Technology Transfer for a quote and licensing options.

Drew Bennett - https://techtransfer.umich.edu/team/drew-bennett/

or

techtransfer@umich.edu
