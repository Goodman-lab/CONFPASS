# CONFPASS
Conformer Prioritizations &amp; Analysis for DFT reoptimisations


version: v4 25042023


## A.	Introduction 

There are two parts in this package:
PART 1: class conp (pre-reoptimisations: clustering and priority list generation) 
PART 2: class pas (post-reoptimisations: predict the completeness of the reoptimisation process) 

This package is suitable for processing organic molecules, including radical species (excluding organometallics and inorganic compounds). Explanations and theories behind the script can be found in the paper. 

The package is written entirely in python and developed under the following environment:

-	The Python (3.8.12) Standard Library (os, sys, itertools, collections, pickle, random, optparse, traceback, json)
-	Pandas (1.3.4)
-	Numpy (1.21.2)
-	Sklearn (1.0.1)
-	RDkit (2019.09.3)
-	Natsort (8.0.2)


A logistic regression model is used in PART 2, class pas.

We used the script from https://github.com/jensengroup/xyz2mol/blob/master/xyz2mol.py (which is based on the work of Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777) to convert xyz coordinates to a rdkit.Chem.mol object


## B.	Installation

* Option 1: via pypi (To be available soon)
* Option 2: Manually Cloning the repository and then adding the location of the compass directory to the PYTHONPATH environment variable.


## C.	Citing CONFPASS 

1.	C. C. Lam, J. M. Goodman, CONFPASS: fast DFT re-optimisations of structures from conformation searches; just accepted by Journal of Chemical Information and Modelling: 10.26434/chemrxiv-2023-vhlgg


## D.	Usage - Option 1: Execution as a command line


i.e. python -m confpass 

### Options:

* `-h`, --help: show the help message and exit


* `-m` M: priority list assembling method 
(choices=('pipe_x_as', 'pipe_x', 'pipe_as', 'nth', 'random', 'ascend')); pipe_x_as is equivalent to 'pipeline-mix', which is the notation used in the paper.


* `-x` X: hyperparameter for the pipe_x method, default = 0.8
 i.e. use the clustering result when n_clusters=total_conformer_no*x
 
 
* `--x_as` = X_AS:  hyperparameter for pipe_x_as method, default = 0.2 
i.e. 20% (pipe_x priority list) + 80%(pipe_as priority list); x_as is equivalent to Q, which is the notation used in the paper. 


* `-n` N: hyperparameter for nth method, default = 3 
i.e. prioritise every 3rd conformer


* `--csv`: get csv of the result (the priority list data frame or a summary of the pas test) 

* `--p2gjf`: execute priority2gjf() – create gjf files for DFT calculations according to the priority list   

* `-per`=PER: hyperparameter for priority2gjf() function; percentage of the conformers to be converted into gjf in the priority list

* `--togjf`: execute get_gjf(); create gjf files for DFT calculations based on a list of conformer indexes specified by the user in the keywords.py script

(NOTE: The keywords for conducting DFT calculations can be modified by changing lines in the keywords.py script. Please ensure that keywords.py is in the same directory when applying --p2gjf and --togjf. An example keyword script is located in the user_guide/ folder)

* `--path`=PATH: the path for the pas test, the path of the folder that contains the DFT output file from Gaussian16 calculations; the pas test can be activated by adding --pas or --pas_multi to the command line

* `--pas`: predict the completeness of the reoptimisation process

* `--pas_multi`: predict the completeness of the reoptimisation process; the pas test will be performed on all the folders in the --path specified

* `--radical`: this commend needs to be included when processing radical molecules (with a different number of atoms compared to the pseudo structure)

* `--rmatom`: list of atoms to be removed – need to be specified when processing radical molecules (with a different number of atoms compared to the pseudo structure)


## Examples

### PART 1: conp


You will need to go to execute the commands in the directory that contains your sdf files. 
E.g. 
`cd test_demo/`

1.	Generate the priority list for a sdf file from a conformational searching calculation with the default setting -- x=0.8, x_as=0.2, method = 'pipe_x_as'
```python
python -m confpass test_15.sdf 
```
2.	Generate the priority list for multiple sdf files with the default setting.
```python
python -m confpass *.sdf
```
3.	Generate the priority list for multiple sdf files with the default setting. Create a csv file for the priority list data frame. 
```python
python -m confpass *.sdf --csv
```
4.	Generate the priority list for multiple sdf files with the default setting. Create the gjf files for the top 5% of the conformers in the priority list 
```python
python -m confpass *.sdf --per 0.05 --p2gjf 
```
Please ensure that keywords.py is in the same directory as the sdf files. (keywords.py is in the /user_guide folder)

Note on the keywords.py: the keywords for DFT calculations are specified using keywords.py 

You can change:
* keyword – lines before the xyz coordinates
* space – lines after the xyz coordinates
* conf_idx_ls– required for the --togjf function; a list of conformer indexes 


5.	Generate the priority list for the sdf file of a pseudo radical structure with the default setting. Create the gjf files for top 5% of the conformers in the priority list.

If the radical structure and the pseudo structure from conformational searching are different, you will need to specify the removed atoms. This helps to remove the atom with index = 32 in the gjf file.

```python
python -m confpass radical_gu_liu1.sdf --rmatom '[32]' --per 0.05 --radical --p2gjf
```

6.	Generate the priority list for multiple sdf files with the every nth method, where n =5
```python
python -m confpass *.sdf -m nth -n 5
```

7.	Generate the priority list for multiple sdf files with the pipe_x method, where x = 0.7. Create a csv file for the priority list data frame. 
```python
python -m confpass *.sdf --csv -m pipe_x -x 0.7
```

### PART 2: pas 


8.	Perform the pas test on a result folder. The priority list is generated using the default setting.

You will need to go to execute the commands in the directory that contains your result folder. 
E.g. `cd test_demo/`
```python
python -m confpass --path ./test_15 --pas
```
Format of result folder – your files need to follow the below naming format for the script to work:  

```
mol/
├─ spe/
│   ├─ mol_1_spe.out
│   ├─ mol_2_spe.out
│   ├─ mol_3_spe.out
│   ├─ ...
├─ mol.sdf <- the conformational searching output file
├─ mol_1.out
├─ mol_2.out
├─ mol_3.out
├─ ...
```

(opt+)freq output file: Name_idx.out 

Single point energy output file: Name_idx_spe.out 


9.	Perform the pas test on multiple result folders under the same directory. The priority list is generated using the default setting. Generate a csv file to summarise the result. (Please do not include other irrelavent folders in the specified path)
```python
python -m confpass --path ./test_demo --pas_multi --csv
```

10.	Perform the pas test on a result folder for a radical species. The priority list is generated using the default setting.

If the radical structure and the pseudo structure from conformational searching are different, you will need to specify the removed atoms. ‘--pas_multi’ setting is currently not available for this type of radical structures.

```python
python -m confpass --path ./radical_gu_liu1 --rmatom '[32]' --radical --pas
 ```


## E.	Usage - Option 2: import and use in python scripts or Jupyter notebooks 

Please see CONFPASS_tutorial.ipynb for the user guide, tutorials and examples. Download and unzip the test_demo.zip file in user_guide/

Here, we have taken out a simple example from the CONFPASS_tutorial.ipynb as a taster sample:

Example 1:

```python

from confpass import confpass 

## 1. Generate the priority list for a sdf file from a conformational searching calculation with 
## the default setting -- x=0.8, x_as=0.2, method = 'pipe_x_as'

test1=confpass.conp(['test_15.sdf'])
test1.get_priority()

## for conp.priority_df, the index of the conformer starts from 0
test1.priority_df


```

Example 2:

```python

from confpass import confpass 
import os 

## 1. Generate the priority list for sdf files in a folder
## the default setting -- x=0.8, x_as=0.2, method = 'pipe_x_as'
path = os.getcwd()
sdf_ls = [f for f in os.listdir(path) if f.endswith('.sdf')]

test1 = confpass.conp(sdf_ls)
test1.get_priority()

## 2. generate the g16 input files for the first 20% of the conformers in the priority list

keywords = '''%nprocshared=32
%mem=4GB
# opt freq b3lyp/6-31g(d) int=ultrafine empiricaldispersion=gd3

Title Card Required

1 1
'''

space='''
'''

test1.priority2gjf(keywords, space, 0.2)



```

Example 3:

```python
## 8.	Perform the test using the LR model on a result folder. 
## The priority list is generated using the default setting.

path = os.getcwd()+'/test_15'

test6 = confpass.pas(path)
test6.preparation()
test2.make_prediction(repeat=3)

## generate the g16 input file if %conf is less than 85%
if test6.confidence < 85:
    print('NEED FURTHER CALCULATIONS!')
    test6.add_cal( keywords, space, per=0.05)
    
    
    

```

