julia
=====

Julia scripts - 
NOTES: see julialang.org for details on the Julia language.

INVENTORY

conv_strat_partition.jl
A convective/stratiform/weak echo partitioner that analyzes radar reflectivity according to the criteria in 
Appendix A of Didlake & Houze (2009; MWR) based on Churchill & Houze (1984), Steiner et al (1995),
Yuter & Houze (1997). Writes a mask with the classification of size X x Y x Time to a new variable CSMASK
in the original NetCDF file. Requires Julia NetCDF package. 
