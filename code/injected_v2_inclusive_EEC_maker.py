# makes EEC for inclusive datasets
# modified from non_binned_inclusive_EEC_maker.py

# this file takes in 
# 1. a directory containing directories for each batch
# 2. an output directory where the results are saved (name gets modified based on v2 injection amount)
# 3. the batch number to analyse all files for
# 4. the v2 amount to inject (NEW!!)

# output: A root file containing (for both wta and standard, for all files in the batch) 
# - energy profile, injected EEC signal, total Nch, and number of jets

# skips invalid files


# IMPORTS ---
import argparse
import glob
import ROOT
import math
import fastjet
import os
import gc
import awkward as ak
import uproot
import numpy as np
import vector
vector.register_awkward()
# ---

def main():

    # 1. SET UP THE ARGUMENT PARSER AND OPEN THE FILES
    # This acts as the bridge allowing Slurm to feed file paths into Python
    parser = argparse.ArgumentParser(description="Run non-binned, inclusive Energy Correlators on NOTS Cluster")
    parser.add_argument("--input", required=True, type=str, help="Path to the directory of all batches to analyse")
    parser.add_argument("--outdir", required=True, type=str, help="Directory where output should be saved")
    parser.add_argument("--batchnum", required=True, type=int, help="Batch number to analyse")
    parser.add_argument("--v2", required=True, type=str, help="v2 value to inject (pass as a string)")
    args = parser.parse_args()

    file_list = sorted(glob.glob(os.path.join(args.input, f"batch{args.batchnum}", "*.root"))) # Get a list of all root files in that batch directory
    print(f"Found {len(file_list)} files to process in {os.path.join(args.input, f"batch{args.batchnum}")}")
    
    # 2. INITIALISING THINGS
    np.random.seed(42) # set reproducable seed for background
    ROOT.gROOT.SetBatch(True) #Force ROOT into headless batch mode - i.e. turn off pop up graphics

    # read v2 value to inject
    v2 = float(args.v2)

    # Jet definitions
    jet_radius = 1000 #max radius accepted by fastjet is 1000. radius value in prl paper is 0.8 in lab frame. Use 1000 because root data already has the clustered jets
    wta_def = fastjet.JetDefinition(fastjet.antikt_algorithm, jet_radius, fastjet.WTA_pt_scheme) #winner take all definition
    std_def = fastjet.JetDefinition(fastjet.antikt_algorithm, jet_radius) #standard E scheme definition

    # Initialize persistent histogram storage and tracking variables BEFORE the loop
    num_jets_dict = {"wta": 0.0, "std": 0.0} # number of jets per batch

    # Setup histogram bin specs
    deltarBW = 0.001

    # Book empty histograms once up front
    ROOT.gDirectory.Clear()
    wta_histograms = {
        'profile': ROOT.TH1D(f"WTA_profile", f"WTA Energy profile (inclusive);#Delta R;Profile" , int(1/deltarBW),0,1)
        }
    std_histograms = {
        'profile': ROOT.TH1D(f"STD_profile", f"Standard Energy profile (inclusive);#Delta R;Profile" ,int(1/deltarBW),0,1),
        'eec': ROOT.TH1D(f"STD_eec", f"Standard Energy-Energy Correlator (inclusive), v2={args.v2};#Delta R;EEC" ,int(1/deltarBW),0,1)
        }
        

    # -------------------------------------------------------------------------
    # DATA STREAMING LOOP OVER MULTIPLE FILES
    # -------------------------------------------------------------------------
    search_pattern = os.path.join(args.input, f"batch{args.batchnum}", "*.root")
    all_files = sorted(glob.glob(search_pattern))
    
    print(f"Found {len(all_files)} total files across batch {args.batchnum} in {args.input}")
    if len(all_files) == 0:
        raise FileNotFoundError(f"No root files found matching pattern: {search_pattern}")

    # --- NEW: PRE-SCREENING FILTER FOR INCOMPLETE/CORRUPTED FILES ---
    required_branches = {"genJetPt", "genJetEta", "genDau_pt", "genDau_eta", "genDau_phi", "genDau_chg"}
    valid_files = []
    
    print("Pre-screening files to verify required tree and branches...")
    for filepath in all_files:
        try:
            with uproot.open(filepath) as f:
                if "trackTree" in f:
                    tree = f["trackTree"]
                    # check if all required branches are subsets of the tree keys
                    if required_branches.issubset(tree.keys()):
                        valid_files.append(filepath)
                    else:
                        missing = required_branches - set(tree.keys())
                        print(f"Skipping {os.path.basename(filepath)}: missing branches {missing}")
                else:
                    print(f"Skipping {os.path.basename(filepath)}: 'trackTree' not found")
        except Exception as e:
            # Captures files that are partially written, zero-byte size, or corrupted
            print(f"Skipping {os.path.basename(filepath)}: error opening file ({e})")
            
    print(f"Retained {len(valid_files)} / {len(all_files)} valid files for processing.")
    if len(valid_files) == 0:
        raise FileNotFoundError("None of the discovered files contain a valid 'trackTree' with the necessary branches.")

    # 2. Append the TTree name ONLY to verified discovered file paths
    tree_paths = [f"{filepath}:trackTree" for filepath in valid_files]
    print("Beginning vectorized streaming iteration over valid files...")

    chunk_counter = 0
    # Pass the filtered list of paths directly to uproot.iterate
    for data in uproot.iterate(tree_paths, expressions=["genJetPt", "genJetEta", "genDau_pt", "genDau_eta", "genDau_phi", "genDau_chg"], step_size="150 MB"):
        chunk_counter += 1
        print(f"--- Processing Data Chunk #{chunk_counter} ---")

        # PRELIM ANALYSIS (Cuts) ---
        # lab frame jet cuts
        jetEtaCut = 1.6 # from PRL paper page 2
        jetPtCut = 550  # from PRL paper page 2
        data = data[(data.genJetPt > jetPtCut) & (abs(data.genJetEta) < jetEtaCut)] # applying cuts
        data = data[ak.num(data.genJetPt) > 0] # removing empty jets
        if len(data) == 0:
            continue

        # flattening to the relevant jets
        # Flatten the event/jet axes so that axis 0 is simply a list of all passing jets
        jet_dau_pt = ak.flatten(data.genDau_pt, axis=1)
        jet_dau_eta = ak.flatten(data.genDau_eta, axis=1)
        jet_dau_phi = ak.flatten(data.genDau_phi, axis=1)
        jet_dau_chg = ak.flatten(data.genDau_chg, axis=1)
        energies = jet_dau_pt * np.cosh(jet_dau_eta) # calculates total momentum, approximately equal to particle energy. p = pt * cosh(eta)

        # Build a jagged array of 4-vectors for all particles inside all jets
        particles = ak.zip({
            "pt": jet_dau_pt,
            "eta": jet_dau_eta,
            "phi": jet_dau_phi,
            "E": energies 
        }, with_name="Momentum4D")

        del data    # don't need jet level info anymore now that cuts have been made
        gc.collect()

        # 2. lab frame particle cuts
        particlePtCut = 0.3 #from prl paper page 2
        particleEtaCut = 2.4 #from prl paper page 2


        # particle selection
        particle_mask = (particles.pt > particlePtCut) & (abs(particles.eta) < particleEtaCut) & (jet_dau_chg != 0) #charged particles, meeting pt and eta cuts
        particles = particles[particle_mask] #now only valid particles left
        if len(particles) == 0 or ak.sum(ak.num(particles)) == 0:
            continue


        # 3. clustering the lab frame particles
        wta_clustered = fastjet.ClusterSequence(particles, wta_def)
        std_clustered = fastjet.ClusterSequence(particles, std_def)

        wta_jets = fastjet.sorted_by_pt(wta_clustered.inclusive_jets()) # the recombined wta jets
        std_jets = fastjet.sorted_by_pt(std_clustered.inclusive_jets()) # the recombined std jets

        wta_constituents = wta_clustered.constituents() #the constituent particles of the recombined jets
        std_constituents = std_clustered.constituents()

        # using leading jets.
        # fastjet sometimes recombines one jet into two or more, so in those cases i need to use the highest pt recombined one.
        wta_jets = wta_jets[:, -1] #slicing to be only the leading jet in the recombined lists. Other items in the list will only have a couple of soft particles
        std_jets = std_jets[:, -1]

        wta_constituents = wta_constituents[:, -1] # only need constituents of the relevant recombined jets
        std_constituents = std_constituents[:, -1] 

        # garbage collecting:
        del particle_mask, particles
        del wta_clustered, std_clustered
        gc.collect()
        # ---
    
        # SAVING Nch INFORMATION ----
        # 1. sorting by multiplicity bins
        # we get 1) Nch per jet, 2) number of jets in each bin
        # counting Nch for each jet
        wta_Nch = ak.num(wta_constituents) # number of particles in the recombined jets (before making jet frame cuts, after making lab frame cuts)
        std_Nch = ak.num(std_constituents)
        num_jets_dict['wta'] += len(wta_jets)
        num_jets_dict['std'] += len(std_jets)


        # MAKING THE SIGNAL (for this chunk) -----
        
        # Filter the jet-frame particles down to just the jets in this bin
        binned_wta_parts = wta_constituents
        binned_std_parts = std_constituents
        binned_wta_axes = wta_jets
        binned_std_axes = std_jets
        num_jets_in_chunk = len(binned_std_axes)
        psi_vals = np.random.uniform(-np.pi, np.pi, size=num_jets_in_chunk)

        # --- FILL SIGNAL ---
        # 1. WTA (only profile, no EEC)
        if num_jets_dict["wta"] > 0:
            # 1. Generate all unique pairs of particles within each individual jet (axis=1)
            wta_pairs = ak.combinations(binned_wta_parts, 2, axis=1) # gives tuple pairs
            if ak.sum(ak.num(wta_pairs)) > 0: # i.e. if more than 0 pairs
                # Fill in the energy profile:
                wta_deta_profile = ak.flatten((binned_wta_parts.eta - binned_wta_axes.eta), axis=None)
                wta_dphi_profile = ak.to_numpy(ak.flatten(np.arccos(np.cos(binned_wta_parts.phi - binned_wta_axes.phi)), axis=None)).astype(np.float64)
                
                wta_dRL_profile = np.sqrt(np.square(wta_deta_profile) + np.square(wta_dphi_profile))
                wta_dRL_profile_flat = ak.to_numpy(ak.flatten(wta_dRL_profile, axis=None)).astype(np.float64)
                wta_profile_energies = binned_wta_parts.pt / binned_wta_axes.pt
                wta_weights_profile = ak.to_numpy(ak.flatten(wta_profile_energies, axis=None)).astype(np.float64)
                
                if len(wta_dRL_profile_flat) > 0:
                    h_profile_wta = wta_histograms['profile'] #fill the right hist
                    h_profile_wta.FillN(len(wta_dRL_profile), wta_dRL_profile_flat, wta_weights_profile)


        # 2. STD signal
        std_n_sig_pairs = 0
        if num_jets_dict["std"] > 0:
            std_pairs = ak.combinations(binned_std_parts, 2, axis=1)

            # 3. Compute per-particle relative phi and flow weights for STD (for EEC & profile)
            std_parts_phi_wrt_jet = binned_std_parts.phi - binned_std_axes.phi
            std_parts_flow_weight = 1.0 + 2.0 * v2 * np.cos(2.0 * (std_parts_phi_wrt_jet - psi_vals))
            
            # Combine the single particle flow weights into unique pairs matching std_pairs
            std_weight_pairs = ak.combinations(std_parts_flow_weight, 2, axis=1)

            if ak.sum(ak.num(std_pairs)) > 0:
                std_p1, std_p2 = ak.unzip(std_pairs)
                
                std_w1, std_w2 = ak.unzip(std_weight_pairs) # Unzip particle 1 and particle 2 flow weights
                
                std_deta = ak.to_numpy(ak.flatten((std_p1.eta - std_p2.eta), axis=None)).astype(np.float64)
                std_dphi = ak.to_numpy(ak.flatten(np.arccos(np.cos(std_p1.phi - std_p2.phi)), axis=None)).astype(np.float64)
                std_dRL = np.sqrt( np.square(std_deta) + np.square(std_dphi) )
                
                std_jet_pt_sq = binned_std_axes.pt ** 2
                std_energies = std_p1.pt * std_p2.pt
                std_normalised_energies = std_energies / std_jet_pt_sq

                std_total_weights = std_normalised_energies * std_w1 * std_w2

                std_weights = ak.to_numpy(ak.flatten(std_total_weights, axis=None)).astype(np.float64)
                std_n_sig_pairs = len(std_deta)

                h_EEC_std = std_histograms['eec']
                h_EEC_std.FillN(std_n_sig_pairs,  std_dRL, std_weights)

                del std_pairs, std_p1, std_p2, std_deta, std_dphi, std_weights
                gc.collect()

                # Fill in the energy profile:
                std_deta_profile = ak.flatten(abs(binned_std_parts.eta - binned_std_axes.eta), axis=None)
                std_dphi_profile = ak.to_numpy(ak.flatten(np.arccos(np.cos(binned_std_parts.phi - binned_std_axes.phi)), axis=None)).astype(np.float64)
                
                std_dRL_profile = np.sqrt(np.square(std_deta_profile) + np.square(std_dphi_profile))
                std_dRL_profile_flat = ak.to_numpy(ak.flatten(std_dRL_profile, axis=None)).astype(np.float64)
                std_profile_energies = binned_std_parts.pt / binned_std_axes.pt
                std_weights_profile = ak.to_numpy(ak.flatten(std_profile_energies, axis=None)).astype(np.float64)
                
                if len(std_dRL_profile_flat) > 0:
                    h_profile_std = std_histograms['profile'] #fill the right hist
                    h_profile_std.FillN(len(std_dRL_profile_flat), std_dRL_profile_flat, std_weights_profile)


    print("All streaming chunks processed. Persistent histograms populated completely.")

    # making v2 tag for the directory name
    v2_tag = args.v2.replace(".", "_") #replace decimal with underscore

    # SAVING THE OUTPUT TO A SINGLE COMBINED COHESIVE ROOT FILE
    out_dir = os.path.join(args.outdir, f"injected_v2_{v2_tag}")
    os.makedirs(out_dir, exist_ok=True)

    full_output_path = os.path.join(args.outdir, f"EEC_non_binned_Output_Batch{args.batchnum}.root")
    
    print(f"Opening output file to save all results: {full_output_path}")
    out_file = ROOT.TFile(full_output_path, "RECREATE")
    

    # open the hists
    hEEC_std = std_histograms['eec']
    hProfile_wta = wta_histograms['profile']
    hProfile_std = std_histograms['profile']
    
    # Rename safely inside directory
    hEEC_std.SetName(f"hEEC_STD")
    hProfile_wta.SetName(f"hProfile_WTA")
    hProfile_std.SetName(f"hProfile_STD")

    # Write out items
    hEEC_std.Write()
    hProfile_wta.Write()
    hProfile_std.Write()

    param_num_jets_wta = ROOT.TParameter('double')(f"num_jets_WTA", float(num_jets_dict["wta"]))
    param_num_jets_std = ROOT.TParameter('double')(f"num_jets_STD", float(num_jets_dict["std"]))
    param_num_jets_wta.Write()
    param_num_jets_std.Write()

    out_file.Close()
    print(f"Results for batch {args.batchnum} written and saved")

    print("Job successfully completed")
    # END ----


# running:
if __name__ == "__main__":
    main()
