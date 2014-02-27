This python script reads MOL files in the current directory and groups them into clusters based on structural similarity.
It requires the pybel module for reading structures and computing structural fingerprints, and scikit-learn for clustering.
The distance between compounds is modeled as 1-T where T is the Tanimoto coeffecient between the corresponding chemical fingerprints.
The script prints the resulting clusters to stdout.
