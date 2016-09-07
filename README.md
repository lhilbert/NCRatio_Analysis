# NCRatio_Analysis
MatLab scripts to analyze nuclear-over-cytoplasmic intensity ratios, typically out of lightsheet microscopy datasets. Based on on bfmatlab of Open Microscopy Environment's BioFormats.

Created by Lennart Hilbert, Center for Systems Biology Dresden
Research group Nadine Vastenhouw, Max Planck Institute of Molecular Cell Biology and Genetics
Research group Vasily Zaburdaev, Max Planck Institute for the Physics of Complex Systems

These scripts were developed to extract ratios of nuclear fluorescence intensity over cytoplasmic fluorescence intensity in young zebrafish embryos. Our data were collected with a Zeiss Z.1 lightsheet microscope, and stored as .czi multi-color stacks. The analysis can deal with single color as well as multi-color data sets, and analyze single as well as multi-position images.

The data should first be ordered into a folder with subdirectories, where each subdirectory represents one experimental condition or experimental time point.

The code will then find all files in that source directory (recursive search), segment nuclei as well as cytoplasmic shells surrounding the nuclei, and save the intensities, coordinates, volumes ... for all nuclei.

All analysis results are saved within their respective subfolders. These results are stored in a per-stack basis, and can be reordered into new folders. After this, another script can be used to extract the main features out of these ordered folders and store them in excel sheet as well as in .mat format for further processing in MatLab.
