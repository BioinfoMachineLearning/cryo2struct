# Execute data preparation script
echo "Script started for data preparation"
relative_path=$1
density_map_dir=$(readlink -f "$relative_path")
for i in "${density_map_dir[@]}"; do
  echo "Running for density maps present in : $i"
  python3 get_resample_map.py "$i"
  python3 get_normalize_map.py "$i"
  python3 get_atoms_label.py "$i"
  python3 get_amino_labels.py "$i"
  python3 get_secondary_pdb.py "$i"
  python3 get_sec_stru_coil_label.py "$i"
  python3 get_sec_stru_helix_label.py "$i"
  python3 get_sec_stru_strand_label.py "$i"
  python3 merge_chains_fasta.py "$i"
  python3 get_pdb_seq.py "$i"
  python3 make_input_run_clustal.py "$i"
  python3 grid_division.py "$i"
  echo "Performing grid division now."
  wait
  echo "Done with maps in directory: $i"
  echo "ALL DONE!" 
done
