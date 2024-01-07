# TargetedMetabolomics
R code for targeted metabolomics

The file "MultipleChrom.R" is an R file that can be run on the command line. 
You can use this program if you have a large amount of target ions that you want to find in your MS data. 
This program automates the process by iterating through a data frame of target ions, scanning the MS data for peaks of those ions, making an extracted ion chromatogram from the data, and creating m/z spectra from the identified peaks. 

## Inputs 
The inputs should be added after the MultipleChrom.R file in the command line in the following order.

- Directory containing the MS files in mzML format
- Data frame of target ion masses
- Retention time start (in seconds)
- Retention time end (in seconds)

#### Example
` Rscript MultipleChrom.R  MSDataExample  TargetIonExample.csv  0  1000  `

## Outputs 
The program will produce folders for each target ion. The folders will be labeled from column 1 of the data frame of target masses and will include the target mass (ie pregnenolone_317.2475). In the folder, there will be files and folders named after the file names of the MS .mzML files. There will also be a file called "peaks.csv" file that identifies all the peaks in each MS file. 

- .pdf file - This is an extracted ion chromatogram (EIC)
- .csv file - This is an Excel readable file of intensity and retention times for recreating the chromatogram in your program of choice

There will also be a folder labeled "Compound_IonMass_Spectra". If there are any spectra under the peaks whose parent ion matches the target ion then the MS2 will be deposited in this folder. Two files will be created and labeled with the retention time of the spectra in seconds. 

- .pdf file - This is a labeled m/z spectra
- .csv file - this is an Excel readable file of the spectra

#### Example output
```
- WorkingDirectory
  - Compound1_319.2633_mz
    - MSFile1_EIC.csv
    - MSFile1.pdf
    - MSFile1_Spectra/
      - MSFile1_RT1_.csv
      - MSFile1_RT1_.pdf
      - MSFile1_RT2_.csv
      - MSFile1_RT2_.pdf
    - MSFile2_EIC.csv
    - MSFile2.pdf
    - MSFile2_Spectra/
    - peaks.csv
  - Compound2_315.2319_mz
```

## Explanations of Inputs 
#### Directory of MS files
The command line only needs the directory of the files. The program will then scan the folder for files ending in ".mzML". To convert the .raw files generated from mass specs you can use proteowizard's MSConvert program. You can also use the GPNS website to convert them as well. The following link has a lot of great resources for data conversion. 

[GNPS File Conversion Tool](https://gnps-quickstart.ucsd.edu/conversion)

[How to use MSConvert via GNPS](https://ccms-ucsd.github.io/GNPSDocumentation/fileconversion/)

Click the link and scroll down to the section labeled "Conversion with MSConvert"

#### Data frame of target ion masses
The table below gives you an idea of what the data frame needs to look like. 
The first column is the name of the compound of interest. 
All subsequent columns will be the ions that will be searched in the MS data. The ions don't need to be the same as displayed here. You can change them based on the needs of your project. The program will only look for ions in columns 2 and up. The program will use column 1 to make folders in the working directory that the chromatograms and spectra in which the data will be placed. Make sure the compound name has no special characters or spaces

| Name       | M+H        | M+H-H2O  |
| ------------- |:-------------:| -----:|
| Pregnenolone   | 317.2476 | 299.2375 |
| Progesterone  | 315.2319     |   297.2212 |
| 3B_hydroxy_5B_pregnane_20_one| 319.2633      |  301.2531 |

#### Retention Time 
Retention time that you would like to scan in seconds 

