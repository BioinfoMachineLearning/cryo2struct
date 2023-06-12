# Execute script
echo "Script started for HMM"
a=(22898)
for i in "${a[@]}"; do
  echo "Running with density: $i"
  python3 ../cryo2struct.py --density_map_name "$i"
  wait
  echo "Done with density: $i"
done
