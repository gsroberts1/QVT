MATLAB-package
=========
![QVT Splash](files/splash2.png)

**Current Maintainer: Grant S. Roberts**

Previous Maintainer: Carson A. Hoffman


### Citations ### 
If you are using the QVT for cranial 4D flow analysis in your study, please cite the following papers:

- [Schrauben, E., Wahlin, A., Ambarki, K., Spaak, E., Malm, J., Wieben, O., & Eklund, A. (2015). Fast 4D flow MRI intracranial segmentation and quantification in tortuous arteries. J Magn Reson Imaging, 42(5), 1458-1464. doi:10.1002/jmri.24900](https://pubmed.ncbi.nlm.nih.gov/25847621/)

- Hoffman, C., Roberts, G. S., Berman, S. E., Eisenmenger, L. B., & Wieben, O. (2019). Towards Automated Cranial 4D Flow Cranial Analysis. Paper presented at the Society for Magnetic Resonance Angiography. p.80.

### License ###
BSD 2-Clause


## Installation ##
Requires MATLAB version > 2018

Download or clone the 'QVT' repository into the directory of your choice (e.g., 'C:\Users\username\Documents\MATALB'). From MATLAB, add the 'QVT' folder to your Matlab search path. This can be done several ways:
1. Move to the directory containing the 'QVT' folder, right click on 'QVT', select "Add to Path >> Selected Folders and Subfolders". 
2. `>> addpath(genpath('C:\Users\username\Documents\MATALB\QVT'))`

### Dependencies ###
Machine Learning Toolbox /
Optimization Toolbox \
Image Processing Toolbox


## Usage ##
(Optional) In Matlab, change to the directory where the 4D flow data exists. This is not necessary but is conventient for locating 4D flow data and for accessing saved data after processing.

From the command window, type the following command to open the GUI:
`>> paramMap`

Once opened, select 'Load Data'. From the pop-up window, select the folder which contains 4D flow data.

### IMPORTANT NOTE: Currently, this data must be in a format specific to our institution ###
From the PCVIPR reconstruction, data may be in .dat format (multiple .dat files of containing 3D volumes of magnitude, complex difference, and velocity data) or in HDF5 format (single file usually named 'Flow.h5'). Both formats can be loaded into the tool with the 'loadpcvipr.m' and 'loadHDF5.m' functions. 

In the near future, we plan to implement functions to load more universal 4D flow data formats (e.g., DICOM series, NIFTI?, or ISMRMRD?) from other institutions into our tool. However, this is currently not possible. If you have data from outside of UW-Madison and would like to use the QVT, please reach out and we can help develop functions to load in this data.


## Additional Resources ##

[Video Demo](https://mediaspace.wisc.edu/)

