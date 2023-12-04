# TargetedMetabolomics
R code for targeted metabolomics

The file "MultipleChrom.R" is an R file that can be run on the command line. 
You can use this program if you have a large amount of target ions that you want to find in your MS data. 
This program automates the process by iterating through a data frame of target ions, scanning the MS data for peaks of those ions, making an extracted ion chromatograph from the data, and creating m/z spectra from the identified peaks. 

## Inputs 
The inputs should be added after the MultipleChrom.R file in the command line in the following order.

- Directory containing the MS files in mzML format
- Data frame of target ion masses
- Retention time start (in seconds)
- Retention time end (in seconds)
- Data frame describing how the files are related to the sets in the experiment 

` practice`
