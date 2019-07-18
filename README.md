# Simulation code of a flat lipid membrane with additional protein in LAMMPS


(LAMMPS version 30Jul16)

1. To generate the initial configuration of a flat membrane with two lipid species on an ideal triangular lattice run:
```bash
awk -f flat_membrane_generator_mod.awk /dev/null > membrane-2species.txt
```
Note: to change the proportion of the lipid species modify line 115 of `flat_membrane_generator_mod.awk` (rand()<0.875)
 
2. To generate the potential table with the desired `w_c` parameter run:
```bash
./lammps_table_potentialcos2_creator_v001.awk > table_potential_tail_tail_wC_1_75.txt
```
Note: change line 20 to tune w_c: w_c=1.75*sigma

3. To thermalise the membrane (heating up slowly up to a certain temperature) run the LAMMPS script: lipid_membrane_2species.lmp

4. To read the restart from heating simulation and introduce the protein run the LAMMPS script: lipid_membrane_2species_readrestart.lmp

