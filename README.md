All codes of this project are running under Python 3.7.0 version.
# Download needed files and put them into proper folders
1. Download "BIOGRID-ORGANISM-Homo_sapiens-3.5.179.mitab" file from Biogrid and put it under `data` folder.
2. Download "uniprot\_sprot\_human.dat.gz" file from Uniprot and put it under `data` folder.
3. Download "GSE9476\_RAW.tar", "GSE27567\_RAW.tar" and "GSE121248\_RAW.tart" from GEO and put then under `R` folder. Then decompress them into right folders under `R` folder.
# Run R scripts
Before running R scripts, you need to modify the path variables according to your file system. Every scrip defined the `file_path` and `out_path` variables. You should use your path to replace them before you run the scipt.
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

# Generate compare report between our method and Hum-mploc 3.0
Run the `gene_performance_humploc.py` under `code` folder to generate performance report.

# Generate spreadsheets 
Run the `gene_excel.py` under `code` folder to generate files with extension '.xlsx'. The script will generate 3 files under folder `tem_data`, including "all.xlsx", "39.xlsx" and "max\_min.xlsx". The file "all.xlsx" contains 3 sheets, which correspond to Table s1-s3 in the supplementary materials. The file "39.xlsx" contains 36 sheets, which correspond to Table s4-s39 in the supplementary materials. The file "max\_min.xlsx" contains 6 sheets, which correspond to Table s40-s45 in the supplementary materials. Some files need to manually add headers, add notes and sort operations in Excel program.
