**** NOTE: This repository was changed to the tpc_refitter repository after I migrated my physics work to the DESY BitBucket service for Belle II.  Therefore, this package is obsolete.  It is being kept for history and version control, but is no longer used ********

TPC Data Skimmer:

This is a program written to handle the BASF2 output ntuples produced by Igal
analysis/simulation suite and process them with only the minimal information
needed for doing TPC analysis.

Full variable list being used for analysis (will keep udpated):

nPoints

Notes on variables:
Energy (sum_e) is *not* linear with sum_tot!

Hitside array correspoinds to:
[Left (col<2), Right (col>78), Bottom (row<2), Top (row>334)]

*Need to verify that these values are correct*
