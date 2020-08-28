All codes of this project are running under Python 3.7.0 version.
# Generate supplementary files
Run the script `gene_large_npy_files.py` under `code` directory to generate the necessary intermediate matrix files, including `a.npy`, `ECC.npy` and `PCC.npy`. These files will be saved in `tem_data` directory. Some functions may use these files.

# Generate figures
Run the script `gene_figures.py` under `code` directory to generate figures shown in paper. 

# Generate compressed prediction probability matrix files
Run the script `gene_probas.py` under `code`. This script will create three folders under the `tem_data` folder, including `hepatitis_probas`,`breast_probas` and `leukemia_probas`. Then the corresponding "npy" files will be saved under each directory, the naming form is `probas_nor_tau_x.x.npy` or `probas_ill_tau_x.x.npy`. The former is the probability matrix of healthy samples, the latter is the probability matrix of cancer samples, `x.x` is the value of the hyperparameter `tau`, with a value range of 0.1-2.0.

# Generate textual prediction probability matrix files
Run the script `gene_txt_results.py` under `code`, this script will convert the compressed probability matrix generated above into a text form, save them in the folder `tem_data/probas_txt`, and name them as `cancer_name_nor_probas_x.x .txt` or `cancer_name_ill_probas_x.x.txt`, `cancer_name` is the name of the corresponding disease. `x.x` has the same meaning as described above.
# Generate performance reports
Run the `gene_prediction_record.py` program under `code` to generate the corresponding performance report.
```
python3.7.0 gene_prediction_record.py cancer_name
```
The optional values of the `cancer_name` parameter are 'hepatitis', leukemia' and 'breast'.

# Generate probability difference matrix files and mislocation protein prediction report 
Run the `gene_location_change_record.py` program under `code`. This program will generate the result files in the `tem_data/cancer_name_records/` folder, where `cancer_name` is the name of the cancer. This program will generate files in three forms, namely `diff_x.x.txt`, `diff_x.x.npy` and `record_x.x.txt`. The first two kinds are difference matrix files, the files with suffix "npy" are used for calculation, and the files with suffix "txt"  are used for reading. The last kind sorts the difference matrix from largest to smallest, in which each row is a record. Every row has three columns, corresponding to the name of the protein, the location name and the probability difference.

# Quickly query the mislocalization score of the protein
Run the `query_alternative.py` program under `code`：
```
python3.7.0 query_alternative.py query_file
```
The protein to be queried is stored in the `query_file` file, both "tsv" and "csv" format are OK.

