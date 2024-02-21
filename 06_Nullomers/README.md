# Optimise for maximum fragility: proof of concept
- Generate random sequences using combinations of only nullomers
- Extract feature matrices
- Predict DNA fragility at every base position
- Compare the fragility against a group of control sequences. These are randomly selected from the human genome, like how we do it for the proof of concept project.

# Optimise for maximum fragility: ROptimus
- Generate a random sequence
- Extract feature matrices
- Predict DNA fragility at every base position
- Update the original sequence using ROptimus or other general simulated annealing ways.

# Optimise for maximum fragility: Most flexible approach
- Generate a random sequence
- Extract feature matrices
- Predict DNA fragility at every base position
- Get SHAP values for each feature and see the biggest positive contributors to the fragility. Every feature needs to be mapped back to a k-mer for this to work. Optimise the sequence by changing a single base, but it should change those based on the k-mers of the features with the biggest positive contribution.