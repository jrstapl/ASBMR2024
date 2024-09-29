# ASBMR 2024
## Joshua Stapleton, Soren Gam, Steven K. Boyd, Marjorie Howard, Daniel P. Beaver, Kristen M. Beavers, Ashley A. Weaver
## Validation Of Muscle Segmentation Program For High Resolution Peripheral Quantitative Computed Tomography

## Abstract 3064
## Poster SUN-LB 601


# Introduction
- HRpQCT can be used for soft tissue analysis
- This repository contains a similar style analysis program written in python with SimpleITK
- The environment can be set up using [anaconda](https://www.anaconda.com/)

# Methods
The original image is first downsampled in all dimensions by a factor of 3x, followed by isolation of the limb region of interest (ROI) through segmentation mask extraction. The skin is then segmented from the soft tissue mask and removed, this segmentation is also used to remove fat and muscle segmentations in the skin region. Bones are then segmented and removed from the image based on a binary threshold and mask closing step. Binary thresholding is then used to set fat and muscle mask starting seeds which are cleaned for small, unconnected islands which are likely to be noise. An iterative closing method is implemented to close both fat and muscle masks while preserving independent boundaries. These final masks are then cleaned by removing sections in the skin, and removing small areas of segmentation.

## Limb ROI Extraction
The limb area is extracted through the use of Otsu's algorithm [implemented](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1OtsuThresholdImageFilter.html) in SimpleITK. In brief, this image filter creates a binary mask which separates the foreground and background using a histogram separation method of the greyvalues. From this, the largest connected mask region is selected based on total volume and used to select the smallest bounding box which contains the full limb mask. The ROI is then used for the subsequent steps of mask creation.


## Skin segmentation and removal
Prior to the iterative skin mask removal, the image is smoothed using a smoothing recursive gaussian [filter](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1SmoothingRecursiveGaussianImageFilter.html) with a sigma [value](https://github.com/SpectraCollab/ORMIR_XCT/blob/main/ormir_xct/segmentation/ipl_seg.py) equal to 1/2 the voxel size. A contour is extracted from the outer edge of the limb mask which is used as a starting point for the skin segmentation. This mask is thickened by 15 voxels (kernel size of 15) initially, and the median intensity value of the label is extracted. A cutoff value is set at 1/3 of this number to establish a transition of fat to skin in the mask. This value was established through manual testing and examination of different discriminators for isolating the skin in a single mask. The kernel size is then decreased by one and the median label intensity is extracted at each kernel size and checked to see if the intensity has crossed the cutoff threshold. The final mask is then used to remove this section of the image from the ROI. 

## Bone segmentation and removal
Bones are initially segmented using a binary threshold for all voxels with HU > 300. An opening [filter](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1BinaryMorphologicalOpeningImageFilter.html) is then used to remove small voxel islands which may be noise in the image. The remaining masks are then relabled in order of descending volume (i.e. mask0 = background, mask1 = largest, ... maskN = smallest). The largest remaining nonzero mask is selected as the larger bone (in the distal leg this would be the tibia) and closed using a morphological dilation [filter](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1BinaryDilateImageFilter.html) followed by erosion [filter](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1BinaryErodeImageFilter.html). The same steps are then repeated with the second largest mask to segment the smaller bone. These segmentations are then combined into a single mask and inverted to be removed from the limb ROI. 


## Muscle and fat threshold seeds
Muscle and fat seeds are segmented using a binary threshold of [50 < HU < 600] and [-599 < HU < -150] respectively following a recursive Gaussian smoothing step. A connectivity [filter](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1ConnectedComponentImageFilter.html) is applied and used to remove an connected components containing less than 50 voxels.   

## Iterative closing proceedure
The iterative closing proceedure starts with the fat and muscle seed values established, and applies a single voxel dilation step to each mask independently. The resulting muscle mask is then subtracted from the fat mask, and the resulting fat mask is subtracted from the muscle mask. This step is then repeated 10 times, with a binary threshold at the end of each iteration to ensure each mask remains binary. Following the closing proceedure, the bone mask and skin mask are subtracted from the resulting masks. 

## Cleaning mask and finalizing results
The mask finalization proceedure starts by removing small islands which make up less than 1% of the total mask volume. The remaining mask then goes through an opening and closing proceedure to remove any remaining small components while preserving the overall connectivity of the mask. A binary holes [filter](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1BinaryFillholeImageFilter.html) is then used to patch up any remaining holes in the mask. The bones mask and skin mask are then inverted and used again to ensure no soft tissue mask outside of the ROI. These steps are performed for both the muscle and fat mask.

## Statistics calculation
Physical statistics were calculated using the label shape statistics [filter](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1LabelShapeStatisticsImageFilter) and the label statistics [filter](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1LabelStatisticsImageFilter.html). The total volume was calculated for the soft tissue masks and limb mask (without the bones), and average intensity was calculated for each soft tissue mask. A volume ratio was calculated by dividing the soft tissue mask volume by the total limb volume. The average cross sectional area was calculated by dividing the total volume by the physical height (# slices * image spacing) as previously [reported](https://www.sciencedirect.com/science/article/abs/pii/S1094695020301426). 

These final results are currently printed to the screen, however in the validate function used to batch images for this abstract, the results are compiled and saved to a csv file. 

# Usage
This script can be used by setting up the anaconda environment, and then providing the input file or directory to be segmented. 

## Cloning the repository
You can use git to clone the repository by running:
```bash
git clone https://github.com/jrstapl/HrpqctSoftTissueAnalysis.git
```

Alternatively, you can download the zip file to just download the project. This will not allow you to pull new updates from the project however, and each time you'd like to get the most recent version you will need to redownload it. 

If the repository is cloned, you can use the following command to recieve updates:
```bash
git pull
```

## Anaconda environment setup
If anaconda is installed on the system, the environment can be created by opening a terminal in the project directory and running the following command:
```bash
conda env create -f scanco_vtkbone.yml
```


## Running the program
The program can be run from the terminal command line by activating the conda environment and then running the python progam.
```bash
conda activate scanco_vtkbone
```

```bash
python3 HrpqctSoftTissueAnalysis/SoftTissueSegmentation.py -i <input> 
```

An example command would look like this:
```bash
python3 HrpqctSoftTissueAnalysis/SoftTissueSegmentation.py -i //medctr/path/to/img.nii.gz
python3 HrpqctSoftTissueAnalysis/SoftTissueSegmentation.py -i //medctr/path/to/img.isq
python3 HrpqctSoftTissueAnalysis/SoftTissueSegmentation.py -i //medctr/path/to/directory

```
Where the directory holds with .nii.gz or .isq files to be segmented.

# Future Work
In the future, I plan to continue refinement of the soft tissue analysis through further hardening of the program against edge cases. One current known issue is with any soft tissue out of field artifacts. Additionally, I'd like to work on extraction of a separate intramuscular fat (IMAT) region to add an additional measure of muscle quality. 

# Thank you!
If you'd like to use this tool, or have any issues, please reach out to my email:
```
jrstaple@wakehealth.edu
```
This is my first major project with XCT, and ITK so I welcome any feedback and suggestions for improvement! If you'd like to use this, I'd also be more than happy to help get everything setup and see how it works for you. If you have any features you'd like added, feel free to contact me and I can certainly do my best to try and help get what you need.
