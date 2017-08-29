# PDB file parser/data-extractor

A script which takes a **pdb** file as input and retreives and writes the following information to an output text file.

1. `Title/Name` of the protein.
2. Total `length` of the protein/number of residues in the given pdb file.
3. `Number of chains` present in the protein and their `names` (ascending order if the chains are named numerically followed by alphabetical order).
4. All `aminoacid ratios` present in the protein in alphabetical order.  
</t> eg: `Ratio(Leu) = no of leucine present/total length of protein`.
5. The `total count` of any `unknown aminoacids` present.
6. The names of any `ligand molecules` other than `water`.
7. Calculate all possible `phi`, `psi`, `omega` angles for the given pdb file.

#### To Run
`python bioscript.py <filename>`  
where filename is the name of a pdb file - eg: `2wsc.pdb` which is provided. Running the script gives the output in a file named `2wsc_output.txt`. The format matches the format of `Sample_output.txt` which is also provided.
