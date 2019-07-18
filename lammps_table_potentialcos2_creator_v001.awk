#!/usr/bin/awk
# requires input from /dev/null
# 
# This script creates a lammps tabulated pair style potential for use in
# coarse-grained FENE lipid models.
# We use the FLAT-LJ model discussed in "Cooke, Deserno, J.Ch.Ph. 123, 224710, 2005"

BEGIN {

  OFMT = "%10.8f";

  pow_2_1o6 = 1.12246204831

  epsilon = 1.0

  sigma   = 1.0

  r_c     = pow_2_1o6 * sigma   
  #b       = 1.0*sigma
  w_c   = 1.75*sigma

  # cfr. Fig. 1

  minr    = 0.0
  maxr    = r_c + w_c

  pi = atan2(0, -1)

  steps   = 10000
  delta   = (maxr-minr)/steps
  
}

{

}

END{

  # generating potential and force data
  # NB start from nonzero epsilon, otherwise lammps complains!
  cx = 0.0
  
  for ( i=1; i<=steps; i++ )
  {
    cx += delta
    
    if ( cx <= r_c )
    {
      potv[i] = 4*epsilon*( (sigma/cx)**12 - (sigma/cx)**6 )
    }
    if ( cx > r_c && cx <= (r_c+w_c) )
    {
      potv[i] = -epsilon*cos(pi*(cx-r_c)/(2.*w_c))*cos(pi*(cx-r_c)/(2.*w_c))
    }
    if ( cx > (r_c+w_c) )
    {
      potv[i] = 0
    }
    
  }
  
  for ( i=0; i<steps; i++ )
  {
    potvp[i] = -(potv[i+1]-potv[i])/delta
  }
  
  # printing header
  print("# Tabulated Pair Style potential for inter-lipid interaction.");
  print("# Data file generated with awk scripts; original file 03/06/2014 by A.L.S.");
  print("# Model: FLAT-LJ model discussed in 'Cooke, Deserno, J.Ch.Ph. 123, 224710, 2005' ");
  print("");
  
  # printing tail-tail table
  print("PS_TAIL_TAIL_AB");
  print("N", steps);
  
  print("")
  
  
  for ( i=1; i<steps; i++ )
  {
    print( i, (i*delta), potv[i], potvp[i] );
  }
  print( 1000, (steps*delta), 0, 0 );
}
