#!/bin/awk
# requires input from /dev/null

BEGIN {

  OFMT = "%10.8f";

  sqrt3 = 1.73205080757;
  sqrt3div2 = 0.8660254037844386;

  # number of 2-atom cells that compose the lattice, in each x and y directions
  xsize=87; # must be odd, multiple of 3 for good pbc
  ysize=44;

  latsize=1.12; #  triangle lattice size 

  # masses
  a1_mass = 1.0 # t1 head
  a2_mass = 1.0 # t2 tail
  a3_mass = 1.0 # t3 head
  a4_mass = 1.0 # t4 tail

  # fene bond, n1, head-tail1
  bond1_k    = 30.0
  bond1_r0   = 1.5
  bond1_e    = 1.0
  bond1_s    = 0.95

  # fene bond, n2, tail1-tail2
  bond2_k    = 30.0
  bond2_r0   = 1.5
  bond2_e    = 1.0
  bond2_s    = 1.0

  # harmonic bond, n3, head-tail2
  bond3_k    = 10.0
  bond3_r0   = 4.0

}

{

}

END {

# printing header
print("# Flat FENE coarse grained lipidic membrane; first version July 2014 by A.L.S. and G.C.");
print("# Lammps input data file generated with awk script v001.");
print("");

natoms = xsize * ysize * 2 * 2 * 3 ;
print("\t", natoms, "\t atoms");

nbonds = xsize * ysize * 2 * 2 * 3;
print("\t", nbonds, "\t bonds");

print("\t 0\t angles");
print("\t 0\t dihedrals");
print("\t 0\t impropers");
print("");

print("\t 4\t atom types");
print("\t 6\t bond types");
print("\t 0\t dihedral types");
print("\t 0\t improper types");
print("");


xlo = 0.0;
xhi = (xsize*latsize*1.0);
ylo = 0.0;
yhi = (ysize*latsize*sqrt3);
zlo = ((-2.5)*bond1_len) - vert_lipid_sep - 10.0 ;
zhi = ((+2.5)*bond1_len) + vert_lipid_sep + 10.0 ;


print("\t", xlo, xhi, " xlo xhi");
print("\t", ylo, yhi, " ylo yhi");
print("\t", zlo, zhi, " zlo zhi");
print("");

#initialize random function
seed=$RANDOM;
srand(seed);

# generating atom and bond data
# note: assumes atom style: molecular (atom-ID molecule-ID atom-type x y z)

currentatom=1;
currentbond=1;
currentmolecule=1;

for ( i=0 ; i<xsize; i++ ) # runs along x direction, must be less or equal
{
  for ( j=0 ; j<ysize; j++ ) # runs along y direction, must span odd lines
  {
    # 2 atom sets (for each triangle lattice cell)
    for ( as=0; as<=1; as++ )
    {
      if (as == 0)
      {
        cx = (i*latsize*1.0);
        cy = (j*latsize*sqrt3);
      }
      if (as == 1)
      {
        cx = (i*latsize*1.0 + latsize/2);
        cy = (j*latsize*sqrt3 + latsize*sqrt3div2);
      }


    #decide atom types | MH: I changed line 115 to get the ratio of 1:7 of raft and non-raft lipids

	if (rand()<0.875)
        {
	type=0;  #type 0 is the majoritarian
        }
	else
        {
        type=1; # type 1 is only ~10%
        }

      for ( k=0; k<6; k++ ) # runs along z
      {  

        if (type==0)
        {
		# atom coordinates
		atom_cx[currentatom] = cx;
		atom_cy[currentatom] = cy;

		if ( k < 3 )
		{         
		atom_cz[currentatom] = ((k-2.5)*bond2_s) ;
		}
		else
		{
		atom_cz[currentatom] = ((2.5-(k-3))*bond2_s) ;
		}        



		if ( k < 3 )
		{  
		  atom_molecule[currentatom] = currentmolecule;
		}
		else
		{
		  atom_molecule[currentatom] = currentmolecule+1;
		}
		
		# atom types
		if ( k==0 || k==3 )
		{
		  atom_type[currentatom] = 1; # head
		}
		else
		{
		  atom_type[currentatom] = 2; # tail
		}
		
		# atom bonds
		if ( k==0 ) # bottom head
		{
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+1;
		  bond_type[currentbond]=1;
		  currentbond++;
		  
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+2;
		  bond_type[currentbond]=3;
		  currentbond++;
		}
		if ( k==1 ) # bottom tail1
		{
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+1;
		  bond_type[currentbond]=2;
		  currentbond++;
		}
		if ( k==3 ) # top tail2
		{
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+1;
		  bond_type[currentbond]=1;
		  currentbond++;
		  
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+2;
		  bond_type[currentbond]=3;
		  currentbond++;
		}
		if ( k==4 ) # top tail1
		{
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+1;
		  bond_type[currentbond]=2;
		  currentbond++;
		}
        }
        else # if type==1
	{
		 # atom coordinates
		atom_cx[currentatom] = cx;
		atom_cy[currentatom] = cy;

		if ( k < 3 )
		{         
		atom_cz[currentatom] = ((k-2.5)*bond2_s) ;
		}
		else
		{
		atom_cz[currentatom] = ((2.5-(k-3))*bond2_s) ;
		}        



		if ( k < 3 )
		{  
		  atom_molecule[currentatom] = currentmolecule;
		}
		else
		{
		  atom_molecule[currentatom] = currentmolecule+1;
		}
		
		# atom types
		if ( k==0 || k==3 )
		{
		  atom_type[currentatom] = 3; # head
		}
		else
		{
		  atom_type[currentatom] = 4; # tail
		}
		
		# atom bonds
		if ( k==0 ) # bottom head
		{
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+1;
		  bond_type[currentbond]=4;
		  currentbond++;
		  
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+2;
		  bond_type[currentbond]=6;
		  currentbond++;
		}
		if ( k==1 ) # bottom tail1
		{
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+1;
		  bond_type[currentbond]=5;
		  currentbond++;
		}
		if ( k==3 ) # top tail2
		{
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+1;
		  bond_type[currentbond]=4;
		  currentbond++;
		  
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+2;
		  bond_type[currentbond]=6;
		  currentbond++;
		}
		if ( k==4 ) # top tail1
		{
		  bond_a1[currentbond]=currentatom;
		  bond_a2[currentbond]=currentatom+1;
		  bond_type[currentbond]=5;
		  currentbond++;
		}
	}
        
        currentatom+=1; # we can go to the next atom
        
      }
      currentmolecule+=2; # we have assigned 2 molecules in the previous step
      
    }
  }
}   

currentatom+=-1;
currentbond+=-1;
currentmolecule+=-1;


# checking
if ( currentatom != natoms ) print("FUCKED UP ATOM NUMBER: currentatom = ", currentatom, " ; natoms = ", natoms);   
if ( currentbond != nbonds ) print("FUCKED UP ATOM NUMBER: currentatom = ", currentbond, " ; natoms = ", nbonds);   


# now printing atom data

print("Masses");
print("");
print("\t 1 \t ", a1_mass);
print("\t 2 \t ", a2_mass);
print("\t 3 \t ", a3_mass);
print("\t 4 \t ", a4_mass);
print("");


print("Atoms # molecular"); # checks that molecular atom style is selected
print("");

for ( i=1; i<=currentatom; i++ )
{
  print("\t", i, atom_molecule[i], atom_type[i], atom_cx[i], atom_cy[i], atom_cz[i] );
}
print("");

print("Bonds # molecular");
print("");
for ( i=1; i<=currentbond; i++ )
{
  print("\t", i, bond_type[i], bond_a1[i], bond_a2[i] );
}


# end program
}
