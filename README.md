# Xcorr-reimplementation

Peptide identification tool

Usage: cli.py [OPTIONS] SAMPLE_FILENAME PROTEIN_DATABASE

### Arguments

| Var              | Type | Description                                                      |
| ---              | ---  | ---                                                              |                                                     
| sample_filename  | TEXT | The name of the sample file           [default: None] [required] |
| protein_database | TEXT | The name of the protein database file [default: None] [required] |

### Options

| Var  | Type    | Description                                          |
| ---  | ---     | ---                                                  | 
| --p  | INTEGER | Amount of processes to use in parallel [default: 10] |
| --s  | INTEGER | Amount of spectra loading at a time [default: 5000]  |
| --ps |         | Predict the spectrum [default: no-ps]                |