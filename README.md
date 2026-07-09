# Energy Correlators

**Investigating energy correlators**

This project investigates the differences in energy-energy correlators of monte carlo simulated pp collisions for data with and without hadronic final state interaction (hFSI), and for data with artificially injected elliptic anisotropy (v2). 

### Key files & analysis pipeline
The data analysis part of the project was run on Rice's NOTS cluster, using slurm scripts to call python files and analyse large .ROOT datasets. 

#### Pipeline
Note: the file paths in the slurm scripts will need to be adjusted to run properly

**Inclusive data**
1. submit code/slurm_scripts/0mb_inclusive/make_inclusive_EECs_0mb.txt as a .slurm script. This calls code/non_binned_inclusive_EEC_maker.py for the inclusive non-interacting dataset, and produces one output file
2. submit code/slurm_scripts/3mb_inclusive/make_inclusive_EECs_3mb.txt as a .slurm script. This calls code/non_binned_inclusive_EEC_maker.py for the inclusive 3.0mb hFSI dataset, and produces one output file
3. submit code/slurm_scripts/0mb_inclusive/make_v2_01_inclusive_EECs_0mb.txt. This calls injected_v2_inclusive_EEC_maker.py for the non-interacting inclusive dataset, and makes one ROOT file containing an EEC with an injected v2=0.1
4. submit code/slurm_scripts/0mb_inclusive/make_v2_03_inclusive_EECs_0mb.txt. This calls injected_v2_inclusive_EEC_maker.py for the non-interacting inclusive dataset, and makes one ROOT file containing an EEC using an injected v2=0.3

**High Nch data**
1. Submit code/slurm_scripts/0mb_nch60/make_raw_EECs_0mb.txt . This calls code/all_batch_EEC_maker.py for the non-interacting high multiplicity dataset, and produces one output file *per batch of data*. 
2. Submit code/slurm_scripts/0mb_nch60/merge_EECs_0mb.txt. This merges the per-batch results of the previous step, and normalises the combined EEC, producing one output file containing normalised EECs for each multiplicity bin in the non-interacting high multiplicity dataset.

3. similarly submit make_raw_EECs_3mb.txt and then merge_EECs_3mb.txt to make and then merge+normalise the per-batch results of calling all_batch_EEC_maker.py with the interacting high-nch dataset. 

4. call the injected v2-0.1 and v2=0.3 analyses, then merge those results


#### key files
- **all_batch_EEC_maker.py**: this file makes the raw (non-normalised) EECs, *binned by multiplicity* for a given batch of high multiplicity data. The batch should be a directory containing ROOT files of simulated collision data. The output is a ROOT file containing the raw EEC as well as info about number of jets used to make the EEC (this is used for normalisation of the merged EECs later on). 
- **non_binned_inclusive_EEC_maker.py**: this file makes the raw (non-normalised) EECs for a given batch of inclusive data. The output is ROOT file with the raw EEC and number of jets. 
- **code/all_batch_injected_v2_EEC_maker.py**: does the same as all_batch_EEC_maker.py but with injected v2 of a specified value. 
- **code/injected_v2_inclusive_EEC_maker.py**: does the same as non_binned_inclusive_EEC_maker.py but with injected v2
- **merge_energy_correlators.py**: merges the raw, per-batch EECs and then averages over the number of jets (i.e. scales the histogram by 1/num_jets). 



### Plots made
The plots were made locally using jupyter notebooks in code/local
1. Inclusive datasets
-  Interacting and Non-interacting EECs for inclusive datasets
-  Ratio of interacting over non-interacting EEC
-  EECs with injected v2 = 0.1 and v2 = 0.3
-  Ratios of injected v2 EECs over non-interacting dataset

2. High multiplicity datasets
- Non-interacting EEC binned by multiplicity (multiplicity bins are: [60,71], [71,78], [78,91], [91,97], and [97,1000])
- hFSI EEC binned by multiplicity
- Injected v2=0.1 EEC binned by multiplicity
- Injected v2=0.3 EEC binned by multiplicity
- ratio of interacting over non-interacting EECs for respective multiplicity bins
- ratio of v2=0.1 EECs over non-interacting EECs
- ratio of v2=0.3 EECs over non-interacting EECs


Key plots:
- *inclusive dataset EEC ratios*: EnergyCorrelators/output/ratios/inclusive/EEC_rebinned_inclusive_ratio_log_markers_v2.pdf
- *high Nch interacting vs non-interacting ratios*: EnergyCorrelators/output/ratios/high_Nch/rebinned/EEC_rebinned_plot_ratio_3_0_markers.pdf
- *high Nch v2=0.1 vs non-interacting ratios*: EnergyCorrelators/output/ratios/high_Nch/rebinned/EEC_v2_01_rebinned_plot_ratio_markers.pdf
- *high Nch v2=0.3 vs non-interacting ratios*: EnergyCorrelators/output/ratios/high_Nch/rebinned/EEC_v2_03_rebinned_plot_ratio_markers.pdf

### datasets
1. High multiplicity: Contains only jets with Nch > 60
    a. non interacting dataset (i.e. without hFSI)
    b. 3.0 mb hFSI (i.e. hadronic interaction included in simulation data)

2. Inclusive jets datasets
    a. non-interacting dataset
    b. 3.0 mb hFSI inclusive dataset

for v2 injection, the non-interacting datasets (inclusive and high multiplicity) datasets are used.