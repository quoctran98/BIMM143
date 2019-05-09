Class 12: Structural Bioinformatics II
================
Quoc Tran
5/9/2019

Loading 1HSG (HIV-1 protease + indinavir) from PDB

``` r
library(bio3d)
file.name <- get.pdb("1HSG")
```

    ## Warning in get.pdb("1HSG"): ./1HSG.pdb exists. Skipping download

``` r
hiv <- read.pdb(file.name)
summary(hiv)
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

**Q1: What is the name of the two non protein resid values in this
structure? What does resid correspond to and how would you get a listing
of all residue values in this structure?**

“HOH” and “MK1”. “resid” corresponds to the amino acid residue that an
atom belongs to.

Trimming PDB files to seperate protein and ligand PDB files

``` r
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```

all.pdbqt comes from AutoDock Vina (using config.txt)

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

Finding root mean square distance for each ligand in all.pdbqt vs the
original ligand in 1HSG

``` r
ori <- read.pdb("ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.649  4.206 11.110 10.529  4.840 10.932 10.993  3.655 10.996 11.222
    ## [11] 10.567 10.372 11.019 11.338  8.390  9.063  8.254  8.978

These are RMSD for each ligand in all.pdbqt vs the original ligan in
1HSG\!
