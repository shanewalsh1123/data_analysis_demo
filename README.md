# Sample Data Analysis Scripts

Research in turbulence consists of performing efficient high resolution computer simulations of turbulent systems. 
The large data sets output from these simulations need to be statistically analysed to understand the turbulent behaviour of the systems.
I had to create efficient C code to generate the data sets, and clever python post-processing code to extract the important information from these large data sets.
Included here are python scripts I used to analyse the flow of energy in the CHM equation as part of my PhD thesis. From the visualisations generated by these scripts, we were able to determine that at intermediate energy, the system favoured energy transfers to zonal scales. The HDF5 data sets on which the analysis was performed were generated by a separate MPI C code.
