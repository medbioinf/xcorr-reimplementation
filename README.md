# Xcorr-reimplementation

## Peptide identification tool

### Installation

```
#With conda
conda env create -f environment.yml
```

### Usage 

```
#Activate Environment
conda activate xcorr-reimplementation

#Run
python xcorrcli.py [OPTIONS] SAMPLE_FILENAME PROTEIN_DATABASE
```

### Arguments

| Argument         | Type | Description                                       |
| ---              | ---  | ---                                               |                                                     
| sample_filename  | TEXT | The name of the sample file            [required] |
| protein_database | TEXT | The name of the protein database file  [required] |

### Options

| Option | Type    | Description                                                       |
| ---    | ---     | ---                                                               | 
| --p    | INTEGER | Amount of processes to use in parallel [default: cpu_count() - 2] |
| --s    | INTEGER | Amount of spectra loading at a time [default: 5000]               |
| --ps   |         | Predict the spectrum [default: no-ps]                             |