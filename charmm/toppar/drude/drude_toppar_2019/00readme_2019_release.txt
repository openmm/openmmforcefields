Overview of the 2019 polarizable Drude oscillator force field release.
May 2020

## master file and proteins

toppar_drude_master_protein_2019g.str
---------------
   * Includes updated Drude2019 protein parameters (see below)
   * Improved annotation of atom types, nonbond parameters etc. 
   * Delete extraneous parameters
   * comment out Hui-2015 CAL parameters until further
     validation. Remaining Hui et al parameters are maintained due to use
     in the 2019 protein FF tests.
   * Renaming of selected ions to be constent with additive FF.
   * Move RESI NMA to toppar_drude_model_2019g.str

toppar_drude_master_protein_2019f.str
---------------
   * Updated Drude2019 protein force field: 10.1021/acs.jctc.0c00057
     updated alpha and thole for all residues;
     Glu/Asp/Arg's charges are from updated molecular ions;
     HSD/HSE's charges are from new IMID, 4MIM, 4MIE; 
     chi1/chi2 are updated;
     backbone and CMAP are updated:
     NBFIX for anion/cation-pi interactions 
   * Update RESI LSN, based on neutral amine 
   * Update nonbond parameters/some charges for molecular ions;
   * Add nbfix for halogens with MAM1, EAM1, PAM1, IMID, 4MIM, 4MIE

Ref: Lin, F.-Y.; Huang, J.; Pandey, P.; Rupakheti, C.; Li, J.; Roux, B. t.;
MacKerell, A. D., Further Optimization and Validation of the Classical
Drude Polarizable Protein Force Field. Journal of Chemical Theory and
Computation 2020, 16 (5), 3221-3239.

toppar_drude_d_aminoacids_2019g.str:
---------------
   * Add lonepair to aromatic groups for cation-pi interaction; 
     this is based on toppar_drude_d_aminoacids_2019e.str 

## model compounds

toppar_drude_model_2019g.str
---------------
   * RESI DMET (Dimethylether) renamed to DMEE

toppar_drude_model_2019f.str
---------------
   * Added RESI NMA
   * Updated nonbond parameters for RESI NEA and NPA with new NMA parameters

toppar_drude_model_2018d.str:
---------------
   * Update molecular ions 
   * Optimized molecules: MAM1, EAM1, PAM1 (for neutral amine) and
     IMID, 4MIM, 4MIE (for HSD/HSE)
   * Halogen model compounds and parameters added.
   
## carbohydrates

toppar_drude_carbohydrate_2019a.str
---------------
   * Added new residues AGALA, BGALA, AGLCA, BGLCA, AIDOA,
     BIDOA, AGLCNA, BGLCNA, AGALNA, BGALNA, AMANNA, BMANNA
   * Added bonded parameters related to new residues.
   * File format edits by Ase.

Ref: Pandey, P.; Aytenfisu, A. H.; MacKerell, A. D. Jr.;
Mallajosyula, S. S., "Drude Polarizable Force Field
Parametrization of Carboxylate and N-Acetyl Amine Carbohydrate
Derivatives", Journal of Chemical Theory and Computation 2019,
15 (9), 4982-5000


## lipids

toppar_drude_lipid_2017c.str
---------------
   * No changes

## nucleic acids

toppar_drude_nucleic_acid_2017c.str
---------------
   * No changes
