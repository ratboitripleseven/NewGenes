# NewGenes
NewGenes is the pipeline created during my Master thesis

## Setup
Do this if you have conda installed

run once
```
conda create --name newgenes --file env_eda.txt -c conda-forge -c pytorch -c bioconda -c pypi
```

Whenever you want to use this

```
conda activate newgenes
```

## Setting up HGTDB

Download HGTDB
```
python data/utils/HGTDB/download_hgt_db.py
```

Prep HGTDB
```
python data/utils/HGTDB/preprocess_hgt_db.py
```

## Training and Eval a model
If done configuring

Example will create a test model
```
python main.py -c 000000_release_HGBC –m train
```

Eval
```
python main.py -c 000000_release_HGBC –m eval
```

## Cross-Eval
There is an example config that shows how to use cross-eval

file is 230524_HGBC_CV_maxbin_5_weighted_class.yaml

if you check the config file then you see that there are MULTIPLE partition file

```
python main.py -c 230524_HGBC_CV_maxbin_5_weighted_class -m cross_val
```
## Annotating samples
Add sample in sample_sequences
Preferably by creating a folder

Inside a folder should live a partition file

Example under the test folder is given

run
```
python set_partition_file.py –f sample_sequences/test_partition.csv
```

The script will then create a new .csv file for every sample_sequence
and will also produce a new partition file called test_partition_RTR.csv
For every partition given to the script a new partition file will be created with _RTR (Ready-to-run)

Then, to annotate with the previously trained model

```
python main.py –c 000000_release_HGBC –m eval –a sample_sequences/test/test_RTR.csv
```
NOTE: When annotating with a trained model, always set the mode to eval and add the flag -a with the RTR file!

The script will be found under annotated_file_HGT/MODELNAME
model name is whatever you input in -c

Then you can use the snakemake command to annotate the .fasta file found in the folder



