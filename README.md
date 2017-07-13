# nmd-processes
All processes for NMD index analysis, use with kn_tools
## Process A
BED File retrieval and preliminary RNA lists
### Design Files:
* **design_rnalist.csv:** Design file that includes all the names of the RNA after a few filters have been set:
	1. RNA has at least two isoforms
	2. For each RNA isoform, exonCount is the number of exons
		* To ensure that there is at least one pair of alternative exons:
			* max(exonCount) is at least 5
			* min(exonCount) is at least 3
			* max - min is at least 2
	3. Get rid of strange RNA
		* These have multiple strand directions, do not make sense, this is the case if get_rna_dfs (from kn_tools) returns None
* **design_rnalist_submit1:** Design file with the first 2999 RNA from **design_rnalist.csv**
* **design_rnalist_submit2:** Design file with the remaining 1224 RNA from **design_rnalist.csv**
* **design_rnalist_test:** Design file for testing, includes MBNL1 and a minus strand RNA

### SBATCH Files:
* **submit_part_fetch1.sbatch:** SBATCH that does **Process A** for the first 2999 RNA from **design_rnalist_submit1:**
	* Prerequisite: Have an empty directory named data_dir/
	* _Outputs_: For each RNA, in data_dir, example for MBNL1:
		* MBNL1-1_all_data:
			* id_MBNL1.txt - this has the directory ID for the bash to read, prevents overlap
			* MBNL1-1_analysis:
				* ss_MBNL1-1.csv
				* flag_MBNL1-1.csv
			* MBNL1bed1
				* MBNL1_bed_n.txt for n in number of paths, this is referred to ID in analysis CSV's
			* MBNL1bedinfo1
				* MBNL1_info_n.txt for n in number of paths
			* MBNL1made_beds1
				* MBNL1_bed_n.txt_made.txt for n in number of paths

* **submit_part_fetch2.sbatch:** SBATCH that does **Process A** for the remaining 1224 RNA from **design_rnalist_submit2**
* **submit_test_fetch.sbatch:** SBATCH that does **Process A** for the test RNA from **design_rnalist_test**

### Python Files:
* **make_bed.py:** Uses go_to_bed from kn_tools to get the BED files for an RNA for all paths
	* Options: --rna, --data (data filename)
	* Outputs: BED directories (such as MBNL1bed1/ and MBNL1bedinfo1/)
* **fetchFromBed.py:** Fetches sequences from BED files
	* Options: --fetch (requires BED file, directory of sequence files, and output file)
	* Outputs: Fetched directory (such as MBNL1made_beds1/)
* **startstop.py:** Gets information from sequences using RegEx, uses fetch_coords from kn_tools
	* Options: --rna, --id (directory ID)
	* Outputs: CSV File (such as ss_MBNL1-1.csv)
* **setflags.py:** Calculates NMD Indexes using set_flags from kn_tools
	* Options: --rna, --id, --dir (where the ss csv is located)
	* Outputs: CSV File (such as flag_MBNL1-1.csv)

## Process B
Narrowing down NMD candidates

### Python Files:
* **get_candidates.py:** Opens the flag csv created in **Process A**
	* Checks:
		1. If there is a UTR Intron
		2. If there exists at least one NMD index above, and at least one NMD below 50 bases
	* Outputs: TXT Files candidates.txt and strange_rna.txt if flag file does not exist

### SBATCH Files:
* **submit_candidates.sbatch:** Runs **get_candidates.py**

## Process C
Creating exon mod DataFrames

### Python Files:
* **get_mods.py:** Opens candidates.txt and creates the mod dataframe for each RNA
	* Outputs: CSV File (such as MBNL1-1_nmd_ind.csv) and TXT File (such as err_out.txt)

### SBATCH Files:
* **submit_mods.sbatch:** Runs **get_mods.py**

## Process D
Obtain NMD and Mutual Information analysis

### Python Files:
* **get_nmd.py:** Performs pairwise mutual information analysis
	* Outputs: CSV Files for all pairwise Mutual Information (mutual_inf.csv) and Splice Site info (splicy.csv)

### SBATCH Files:
* **submit_nmd.sbatch:** Runs **get_nmd.py**


