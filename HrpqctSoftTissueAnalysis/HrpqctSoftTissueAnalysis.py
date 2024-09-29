from datetime import datetime
import itk
import itk.itkImageFileReaderPython
import matplotlib.pyplot as plt
import numpy as np
from operator import lt, gt
import pandas as pd
from pathlib import Path
import SimpleITK as sitk
import sys
from typing import NewType, Tuple, Iterable, Union
from math import ceil, floor

from imat_extraction import extract_imat

sitkBinaryMask = NewType("sitkBinaryMask",sitk.Image)

class Operator: 
    """Placeholder for operator typehints
    """
    pass

class FileNotSupportedError(Exception):
    """Raised when a filetype passed to a function is not supported
    """


class SitkBinaryMathFilter:
    """Class for performing math operations on binary images and maintaining
    binary labeling after the operation
    """
    def __init__(self):
        pass


    def SetImg1(self,_img: sitkBinaryMask):
        self.img1 = _img
    
    def setImg2(self, _img: sitkBinaryMask):
        self.img2 = _img

    def GetImg1(self):
        return self.img1
    
    def GetImg2(self):
        return self.img2
    

    def Add(self):
        addFilter = sitk.AddImageFilter()
        binFilter = sitk.BinaryThresholdImageFilter()
        binFilter.SetLowerThreshold(1)
        binFilter.SetUpperThreshold(100000)
        binFilter.SetInsideValue(1)
        binFilter.SetOutsideValue(0)

        tmp = addFilter.Execute(self.img1,self.img2)
        
        return binFilter.Execute(tmp)

    def Subtract(self):

        subtractFilter = sitk.SubtractImageFilter()
        binFilter = sitk.BinaryThresholdImageFilter()
        binFilter.SetLowerThreshold(1)
        binFilter.SetUpperThreshold(100000)
        binFilter.SetInsideValue(1)
        binFilter.SetOutsideValue(0)

        tmp = subtractFilter.Execute(self.img1,self.img2)
        
        return binFilter.Execute(tmp)
    
class SoftTissueOutput():
    """Class for holding output data from soft tissue analysis program
    """
    def __init__(self, _mv: float, _mp: float, _csa: float, _tv: float, _fv: float, _mv_tv: float):
        self.mv = _mv
        self.mp = _mp
        self.csa = _csa
        self.tv = _tv
        self.fv = _fv
        self.mv_tv = _mv_tv
    
        self.stats = [
            self.mv,
            self.mp,
            self.csa,
            self.tv,
            self.fv,
            self.mv_tv
        ]


# Doesn't actually update tags apparently
def write_image_with_date(_img: sitk.Image, _name: Union[str, Path], _CreationDate: datetime):
    """Set image metadata to include creation and modification date. Creation date passed
    to function, modified date determined during call. Image written to path specified

    :param _img: Input sitk image
    :type _img: sitk.Image
    :param _name: Filename of image
    :type _name: Union[str, Path]
    :param _CreationDate: Image creation date
    :type _CreationDate: datetime
    """

    name = str(_name) # in case it's path, doesn't affect if already string

    time_format = "%d-%b-%Y %H:%M:%S.%f"

    _img.SetMetaData("CreationDate",datetime.strftime(_CreationDate,time_format))
    _img.SetMetaData("ModifiedDate",datetime.now().strftime(time_format))

    sitk.WriteImage(
        _img,
        name
    )

def date_in_metadata(_img: sitk.Image) -> int:
    """Check if creation date is in metadata of an image and
    return value corresponding to tags present. 

    :param _img: Input sitk image
    :type _img: sitk.Image
    :return: return code for tag present
    :rtype: int
    """
    keys = _img.GetMetaDataKeys()

    if "CreationDate" in keys:
        return 0
    if "ModifiedDate" in keys:
        return 1
    
    else: return -1
    

def remove_skin_from_image(_img: sitk.Image, _mask: sitk.Image,
                            _v: bool = False, _start_kernel: int = 15 
                            ) -> Tuple[sitk.Image, sitk.Image]:
    """Program to iteratively check label statistics to isolate a skin mask
    in the input limb image and mask. Input mask should cover all of limb /
    soft tissue region of interest.

    :param _img: Input greyscale image of limb.
    :type _img: sitk.Image
    :param _mask: Mask which separates foreground (limb) and background of the greyvalue image
    :type _mask: sitk.Image
    :param _v: Verbose, print calculated statistics if True, defaults to False
    :type _v: bool, optional
    :param _start_kernel: Thickness of starting mask for skin, defaults to 15
    :type _start_kernel: int, optional
    :return: Returns the skin mask, and the limb mask without skin
    :rtype: Tuple[sitk.Image, sitk.Image]
    """
    erode = sitk.BinaryErodeImageFilter()
    stats = sitk.LabelStatisticsImageFilter()
    sub_filter = sitk.SubtractImageFilter()
    erode.SetBackgroundValue(0)

    kernel = _start_kernel
 

    while True:

        erode.SetKernelRadius(kernel)
        shrink_mask = erode.Execute(_mask)


        sub_img = sub_filter.Execute(_mask, shrink_mask)
        stats.Execute(_img,sub_img)

        if kernel == _start_kernel: 
            med = stats.GetMedian(1) 
            cut = ceil(med / 3) if med > 0 else floor(med / 3)
            _operator = lt if med > 0 else gt

            if _v:
                print(f"Cut: {cut}")

            

        if _v:
            print(f"Rad Size: {kernel}")
            print(f"Minimum: {stats.GetMinimum(1)}")
            print(f"Median: {stats.GetMedian(1)}")
            print(f"Mean: {stats.GetMean(1)}")
            print(f"std: {stats.GetSigma(1)}")
            print("\n")


                
        if _operator(stats.GetMedian(1),cut): 
            if (kernel == _start_kernel): 
                kernel -= 1
                continue
            break
        kernel -= 1
        if (kernel < 5) : break


    return shrink_mask, sub_img


def get_image_using_mask(_img: sitk.Image, _mask: sitk.Image) -> sitk.Image:
    """Extract sitk image which falls in the region of a specified binary mask.

    :param _img: Greyscale image
    :type _img: sitk.Image
    :param _mask: Mask which overlaps the desired region of the greyscale image
    :type _mask: sitk.Image
    :return: Greyvalue image with only mask ROI remaining.
    :rtype: sitk.Image
    """
    _img_nda = sitk.GetArrayFromImage(_img)
    _mask_nda = sitk.GetArrayFromImage(_mask)
    img = sitk.GetImageFromArray(_mask_nda * _img_nda)
    img.CopyInformation(_img)
    return img

def read_scanco(_ImagePath: Path, _preprocess_name: bool = True) -> Tuple[itk.image,dict]:
    """Read scanco image using ITKScancoIO and convert to sitk image. 
    Additionally extract and return image metadata dictionary

    :param _ImagePath: Filepath to input image
    :type _ImagePath: Path
    :return: Return sitk image and metadata dictionary
    :rtype: tuple[itk.image,str]
    """
    if _ImagePath.suffix in [".isq", ".aim"]:
        image_type = itk.Image[itk.ctype("signed short"), 3]
        reader = itk.ImageFileReader[image_type].New()
        io = itk.ScancoImageIO.New()
        reader.SetImageIO(io)
        reader.SetFileName(str(_ImagePath))
        reader.Update()

        img = reader.GetOutput()

        if _preprocess_name:
            img["PatientName"] = _ImagePath.stem

        return itk_to_sitk(img), dict(img)

    else: 
        print(f"Unknown file extension: {_ImagePath.suffix}")
        print(f"With file:\n{_ImagePath}")
        raise(FileNotSupportedError("Please use .isq or .aim file for scanco image type"))
    

def sitk_meta_to_dict(_reader: sitk.ImageFileReader) -> dict:
    """Gather sitk image metadata into a dictionary

    :param _reader: Reader used to read SimpleITK image
    :type _reader: sitk.ImageFileReader
    :return: Image metadata dictionary
    :rtype: dict
    """
    meta = dict()
    for key in _reader.GetMetaDataKeys():
        meta[key] = _reader.GetMetaData(key)

    return meta


def read_img_and_meta(_file: Path) -> Tuple[sitk.Image, dict]:
    """Read and return image and metadata with sitk

    :param _file: Filepath of image
    :type _file: Path
    :return: sitk Image and metadata dictionary
    :rtype: Tuple[sitk.Image, dict]
    """
    reader = sitk.ImageFileReader()
    reader.SetFileName(str(_file))
    reader.LoadPrivateTagsOn()
    reader.ReadImageInformation()
    img = reader.Execute()


    meta = sitk_meta_to_dict(reader)

    return img, meta
        


def itk_to_sitk(_imgITK: itk.image) -> sitk.Image:
    """Convert ITK image to SITK image. 

    :param _imgITK: ITK Image
    :type _imgITK: itk.image
    :return: SITK Image
    :rtype: sitk.Image
    """

    origin = _imgITK.GetOrigin()
    spacing = _imgITK.GetSpacing()
    direction = _imgITK.GetDirection()

    itk_nda = itk.GetArrayFromImage(_imgITK)
    sitk_img = sitk.GetImageFromArray(
        itk_nda,
        isVector = _imgITK.GetNumberOfComponentsPerPixel() > 1
        )
    sitk_img.SetOrigin(np.asarray(origin))
    sitk_img.SetSpacing(np.asarray(spacing))
    sitk_img.SetDirection(np.asarray(direction).ravel()) # needs to be flat

    return sitk_img


def sitk_cc_keep_largest(_img: sitk.Image) -> sitk.ConnectedComponent:
    """Extract and relabel connected components of binary image, return largest.

    :param _img: Binary mask
    :type _img: sitk.Image
    :return: Largest connected component of binary mask
    :rtype: sitk.ConnectedComponent
    """
    _cc = sitk.ConnectedComponent(_img)
    _sorted = sitk.RelabelComponent(_cc,sortByObjectSize=True)
    
    largest = _sorted == 1
    del _cc, _sorted

    return largest





def myshow(_img: sitk.Image):
    """Show single slice of 3D image.

    :param _img: Image to be plotted
    :type _img: sitk.Image
    """

    plt.imshow(
        sitk.GetArrayFromImage(_img)[84,:,:],
        cmap="Greys_r"
    )

def remove_small_islands(_binImg: sitk.Image, _volThresh: float = None, _pixNum: float = None) -> sitk.Image:
    """Remove isolated islants in a binary image using either a 
    volume percentage threshold or a number of pixels threshold.

    :param _binImg: Binary mask image.
    :type _binImg: sitk.Image
    :param _volThresh: Percentage of total volume to use as a minimum threshold
    :type _volThresh: float, optional
    :param _pixNum: Number of pixels to use as a minimum cutoff threshold
    :type _pixNum: float, optional
    :return: Binary mask image with small islands removed
    :rtype: sitk.Image
    """

    assert bool(_volThresh) ^ bool(_pixNum), "Please pass either _volThresh OR _pixNum"

    stats = sitk.LabelShapeStatisticsImageFilter()
    relabCCFilt = sitk.RelabelComponentImageFilter()
    agg = sitk.AggregateLabelMapFilter()
    labelToBin = sitk.LabelMapToBinaryImageFilter()

    cc = sitk.ConnectedComponent(_binImg)
    
    if _volThresh:
        stats.Execute(cc)
        sizes = dict()
        for label in stats.GetLabels(): 
            sizes[f"{label}"] = stats.GetNumberOfPixels(label)
        total = sum(sizes.values())
        frac = ceil(total * (_volThresh / 100) )
        relabCCFilt.SetMinimumObjectSize(frac)

    if _pixNum:
        relabCCFilt.SetMinimumObjectSize(_pixNum)

    cleaned = relabCCFilt.Execute(cc)
    aggregated = agg.Execute(sitk.Cast(cleaned,sitk.sitkLabelUInt8)) 
    bin_aggregated = labelToBin.Execute(aggregated)

    del cc, stats, relabCCFilt, agg, labelToBin, cleaned, aggregated

    return bin_aggregated



def finalize_mask(_mask: sitk.Image, _volThresh: float = 1.0, _rad: int = 8) -> sitk.Image:
    """Perform opening and closing steps with small island removal
    from masks. 

    :param _mask: Binary mask to be filtered.
    :type _mask: sitk.Image
    :param _volThresh: Percentage total volume to be used in small island removal, defaults to 1.0
    :type _volThresh: float, optional
    :param _rad: Kernel size used in opening and closing steps, defaults to 8
    :type _rad: int, optional
    :return: Finalized binary mask image
    :rtype: sitk.Image
    """
    opener = sitk.BinaryMorphologicalOpeningImageFilter()
    opener.SetKernelRadius(_rad)

    closer = sitk.BinaryMorphologicalClosingImageFilter()
    closer.SetKernelRadius(_rad)


    cleaned = remove_small_islands(_mask,_volThresh)


    o = opener.Execute(cleaned)
    c = closer.Execute(o)

    holes = sitk.BinaryFillholeImageFilter()
    
    filled = holes.Execute(c)

    del opener, closer, cleaned, o, c

    return filled


def verbose_out(_section: str):
    """Helper function to print section name.

    :param _section: Name of section
    :type _section: str
    """
    print(f"::::{_section}....\n")




def fill_diaphyseal_bone_area(_boneMask: sitkBinaryMask,
                   _limbMask: sitkBinaryMask) -> sitkBinaryMask:
    """Fill holes to create solid bone mask through inversion of the
    trabecular space and a binary hole filling.

    :param _boneMask: Binary mask of bone segmentation
    :type _boneMask: sitkBinaryMask
    :param _limbMask: Mask of full limb
    :type _limbMask: sitkBinaryMask
    :return: Solid binary bone mask
    :rtype: sitkBinaryMask
    """
    inv = sitk.InvertIntensityImageFilter()
    sub = sitk.SubtractImageFilter()
    add = sitk.AddImageFilter()
    mult = sitk.MultiplyImageFilter()
    inv.SetMaximum(1) # ensure binary image result


    invImg = inv.Execute(_boneMask)

    invLimb = mult.Execute(invImg,_limbMask)
    cc = sitk_cc_keep_largest(invLimb)

    holes = sub.Execute(invLimb,cc)

    filled = add.Execute(_boneMask,holes)
    return filled

def fill_distal_bones(_bone_mask: sitk.Image) -> sitkBinaryMask:
    """Perform large closing step to fill holes in bone

    :param _bone_mask: Binary mask of bone isolated from limb
    :type _bone_mask: sitk.Image
    :return: Solid mask of bone area.
    :rtype: sitkBinaryMask
    """
    invFilter = sitk.InvertIntensityImageFilter()
    dilateFilter = sitk.DilateObjectMorphologyImageFilter()
    multFilter = sitk.MultiplyImageFilter()
    openFilter = sitk.BinaryMorphologicalOpeningImageFilter()
    closeFilter = sitk.BinaryMorphologicalClosingImageFilter()
    relabFilter = sitk.RelabelComponentImageFilter()
    ccFilter = sitk.ConnectedComponentImageFilter()

    dilateFilter.SetKernelRadius(1)
    invFilter.SetMaximum(1)
    openFilter.SetKernelRadius(1)
    closeFilter.SetKernelRadius(5)
    relabFilter.SortByObjectSizeOn()

    # roi_mask = sitk.BinaryThreshold(_limb_roi,500,10000)




    mask_open = openFilter.Execute(_bone_mask)
    mask_close = closeFilter.Execute(mask_open)
    mask_inv = invFilter.Execute(mask_close)

    mask_cc = relabFilter.Execute(
        ccFilter.Execute(
            mask_inv
        )
    ) == 1

    hole_closed = closeFilter.Execute(mask_cc)
    hole_solid = invFilter.Execute(hole_closed)

    mask_solid_inv = multFilter.Execute(mask_inv, hole_solid)

    mask_solid = invFilter.Execute(mask_solid_inv)

    return mask_solid

    

def calculate_gauss_sigma(_sig: float, _voxDim: float) -> float:
    """Calculate sigma value for smoothing gaussian filter.

    Per ORMIR XCT

    :param _sig: Original sigma value
    :type _sig: float
    :param _voxDim: Voxel size of image (ideally isotropic)
    :type _voxDim: float
    :return: sigma to be used in gaussian smoothing
    :rtype: float
    """
    return _sig * _voxDim


def _check_minmax(_img: sitk.Image):
    """print minimum and maximum value of an image.

    Used in development.

    :param _img: simple itk
    :type _img: sitk.Image
    """
    print(sitk.GetArrayFromImage(_img).max())
    print(sitk.GetArrayFromImage(_img).min())

def _check_image_stats(_img: sitk.Image):
    """Print image metadata related to size and spatial characteristics

    Used in development.

    :param _img: SITK image
    :type _img: sitk.Image
    """
    print(r"////////////////////")
    print(f"Size: {_img.GetSize()}")
    print(f"Spacing: {_img.GetSpacing()}")
    print(f"Origin: {_img.GetOrigin()}")
    print(r"////////////////////")



def pad_image(_img: sitk.Image, _lower_dims: Tuple[int],_upper_dims: Tuple[float], _const: int = 0) -> sitk.Image:
    """Pad image with given constant to specified dimensions.

    :param _img: Image to be padded
    :type _img: sitk.Image
    :param _lower_dims: Number of padding voxels to be added to lower end of each dimension
    :type _lower_dims: Tuple[int]
    :param _upper_dims: Number of padding voxels to be added to upper end of each dimension
    :type _upper_dims: Tuple[float]
    :param _const: Value of padded pixels, defaults to 0
    :type _const: int, optional
    :return: Image with padded border
    :rtype: sitk.Image
    """

    padImageFilter = sitk.ConstantPadImageFilter()
    padImageFilter.SetConstant(_const)

    padImageFilter.SetPadLowerBound(_lower_dims)
    padImageFilter.SetPadUpperBound(_upper_dims)

    return padImageFilter.Execute(_img)


def create_limb_mask(_img: sitk.Image) -> sitkBinaryMask:
    """Separate limb and background using Otsu's method.

    :param _img: Full greyscale image.
    :type _img: sitk.Image
    :return: Binary mask containing the full limb as the foreground.
    :rtype: sitkBinaryMask
    """
    otsuFilter = sitk.OtsuThresholdImageFilter()
    closeFilter = sitk.BinaryMorphologicalClosingImageFilter()
    otsuFilter.SetInsideValue(0)
    otsuFilter.SetOutsideValue(1)

    seg = otsuFilter.Execute(
        sitk.Cast(_img,sitk.sitkFloat32) # required for otsu
    )

    seg_cc = sitk_cc_keep_largest(seg)
    
    closeFilter.SetKernelRadius(5)
    
    return closeFilter.Execute(seg_cc)


def SoftTissueSegmentation(_InputPath: Path, 
         _OutputPath: Path,
         _fat_thresh: list = [-600,-200],
         _muscle_thresh: list = [100, 600],
         _write_inter: bool = True,
         _preprocess_name: bool = False,
         _v: bool = False):
    """Perform segmentation of fat and muscle by isolating the limb
    from the rest of the image, extracting and removing the skin 
    and bones, and thresholding and closing the fat and muscle. 

    :param _InputPath: Path to input image or directory
    :type _InputPath: Path
    :param _OutputPath: Path to write output images to
    :type _OutputPath: Path
    :param _fat_thresh: Thresholds used for fat segmentation, defaults to [-600,-200]
    :type _fat_thresh: list, optional
    :param _muscle_thresh: Thresholds used for muscle segmentation, defaults to [100, 600]
    :type _muscle_thresh: list, optional
    :param _write_inter: Write intermediate images if True, defaults to True
    :type _write_inter: bool, optional
    :param _preprocess_name: Update patient metadata to be the name of the file, defaults to False
    :type _preprocess_name: bool, optional
    :param _v: Prints out script progress updates if True, defaults to False
    :type _v: bool, optional
    :return: Calculated physical characteristics of finalized masks
    :rtype: SoftTissueOutput
    """

    ##############################
    #####Initialize Filters#######
    ##############################
    andFilter = sitk.AndImageFilter()
    addFilter = sitk.AddImageFilter()
    binFilter = sitk.BinaryThresholdImageFilter()
    clampFilter = sitk.ClampImageFilter()
    clampFilter.SetUpperBound(1)
    clampFilter.SetLowerBound(0)
    ccFilter = sitk.ConnectedComponentImageFilter()
    closeFilter = sitk.BinaryMorphologicalClosingImageFilter()
    dilateFilter = sitk.BinaryDilateImageFilter()
    erodeFilter = sitk.ErodeObjectMorphologyImageFilter()
    smoothingGaussFilter = sitk.SmoothingRecursiveGaussianImageFilter()
    invFilter = sitk.InvertIntensityImageFilter()
    invFilter.SetMaximum(1)
    labelStatsFilter = sitk.LabelStatisticsImageFilter()
    labToBinFilter = sitk.LabelMapToBinaryImageFilter()
    labToBinFilter.SetBackgroundValue(0)
    labToBinFilter.SetForegroundValue(1)
    maskFilter = sitk.MaskImageFilter()
    multFilter = sitk.MultiplyImageFilter()
    negFilter = sitk.MaskNegatedImageFilter()
    openFilter = sitk.BinaryMorphologicalOpeningImageFilter()
    relabFilter = sitk.RelabelComponentImageFilter()
    relabFilter.SortByObjectSizeOn()
    roiFilter = sitk.RegionOfInterestImageFilter()
    shapeStatsFilter = sitk.LabelShapeStatisticsImageFilter()
    subFilter = sitk.SubtractImageFilter()

    ##############################
    #####Initialize Filters#######
    ##############################


    if _v: verbose_out("Reading Image")
    ## Read image to ITK

    filetype = ".".join(_InputPath.suffixes).lower()

    scanco_types = [".isq",".aim"]
    standard_types = [".nii",".nii.gz"]

    assert filetype in scanco_types or filetype in standard_types,(
            f"File type of {_InputPath.suffix} not supported")
    
    if filetype in scanco_types:
        original_image_sitk, img_meta = read_scanco(
                                _InputPath,
                                _preprocess_name = _preprocess_name
                                )
            
        name = f"{img_meta['PatientIndex']}_{img_meta['PatientName']}"
    
    else:
        
        name = _InputPath.stem

    CreationDate = datetime.strptime(img_meta["CreationDate"],"%d-%b-%Y %H:%M:%S.%f")


    output_dir = _OutputPath / name
    output_dir.mkdir(exist_ok=True,parents=True)

    height = original_image_sitk.GetSpacing()[2] * original_image_sitk.GetSize()[2]


    if _write_inter:
        # Write original image
        original_name = output_dir / f"{name}_Original_Image.nii.gz"
        write_image_with_date(
            original_image_sitk,
            str(original_name),
            CreationDate
        )


    # Resample Image
    if _v: verbose_out("Resampling Image")
    img_size = original_image_sitk.GetSize()

    new_size = [int(img_size[0] / 3), int(img_size[1] / 3), int(img_size[2] / 3)]

    ref_img = sitk.Image(new_size,original_image_sitk.GetPixelIDValue())
    ref_img.SetOrigin(original_image_sitk.GetOrigin())
    ref_img.SetDirection(original_image_sitk.GetDirection())
    ref_img.SetSpacing(
        [
            sz * spc / nsz for nsz, sz, spc in zip(new_size, original_image_sitk.GetSize(), original_image_sitk.GetSpacing())
        ]
    )

    resampled = sitk.Resample(original_image_sitk,ref_img)
    voxDim = resampled.GetSpacing()[2]


    if _write_inter:
        resampled_name = output_dir / f"{name}_Resampled_Image.nii.gz"
        write_image_with_date(
            resampled,
            str(resampled_name),
            CreationDate
        )


    ## Create SOLID limb mask for removal of skin,
    # Needs to be fully filled for subsequent stats/erosion steps to work
    if _v: verbose_out("Creating Limb Mask")


    filled_otsu = create_limb_mask(resampled)

    labelStatsFilter.Execute(resampled, filled_otsu)
    bb_limb = labelStatsFilter.GetRegion(1)


    roiFilter.SetRegionOfInterest(bb_limb)
    roi = roiFilter.Execute(resampled)
    roi_mask = roiFilter.Execute(filled_otsu)
    
    if _write_inter:
        limb_mask_name = output_dir / f"{name}_Limb_Mask.nii.gz"
        write_image_with_date(
            filled_otsu,
            str(limb_mask_name),
            CreationDate
            )
        

        
        write_image_with_date(
            roi,
            str(output_dir/f"{name}_ROI.nii.gz"),
            CreationDate
        )



    smoothingGaussFilter.SetSigma(
        calculate_gauss_sigma(1.0,voxDim)
    )
    smoothed = smoothingGaussFilter.Execute(roi)
    if _v: verbose_out("Removing Skin")

    op = gt # seems to work for both, leaving just in case something changes

    shrink_mask, skin_mask = remove_skin_from_image(smoothed,
                                                   roi_mask,
                                                   _v = _v,
                                                   _operator = op
                                                   )
    


    if _write_inter:
        skin_remove_name = output_dir / f"{name}_Skin_Removed_Mask.nii.gz"
        skin_mask_name = output_dir / f"{name}_Skin_Mask.nii.gz"
        write_image_with_date(
            shrink_mask,
            str(skin_remove_name),
            CreationDate
        )

        write_image_with_date(
            skin_mask,
            str(skin_mask_name),
            CreationDate
        )

    if _v: verbose_out("Extracting Bones")
    

    binFilter.SetOutsideValue(0)
    binFilter.SetInsideValue(1)


    no_skin = maskFilter.Execute(roi,shrink_mask)


    smoothingGaussFilter.SetSigma(
        calculate_gauss_sigma(0.8,voxDim)
    )
    smoothed = smoothingGaussFilter.Execute(no_skin)

    # Bone fill needs to be stronger for distal images to close holes

    binFilter.SetLowerThreshold(300)
    binFilter.SetUpperThreshold(100000)
    bones_bin = binFilter.Execute(smoothed)


    openFilter.SetKernelRadius(2)
    bones_open = openFilter.Execute(bones_bin)

    bones_cc = ccFilter.Execute(bones_open)

    bones_relabel = relabFilter.Execute(bones_cc)

    big_bone = bones_relabel == 1

    bone_inv = invFilter.Execute(big_bone)
    
    rest_of_img = sitk_cc_keep_largest(bone_inv)

    bone_center = bone_inv - rest_of_img

    full_bone = bone_center + big_bone

    small_bone = bones_relabel == 2
    small_bone_bin = clampFilter.Execute(small_bone)

    small_bone_inv = invFilter.Execute(small_bone_bin)
    rest_of_img = sitk_cc_keep_largest(small_bone_inv)
    small_bone_center = small_bone_inv - rest_of_img
    full_small_bone = small_bone_center + small_bone_bin

    bones = full_bone + full_small_bone



    bones_name = output_dir / f"{name}_Bones.nii.gz"
    write_image_with_date(
        bones,
        str(bones_name),
        CreationDate
    )


    if _v: verbose_out("Creating Seeds")


    no_bone = negFilter.Execute(no_skin,
                                sitk.Cast(bones,no_skin.GetPixelID()
                                          ))
    
    no_bone_mask = subFilter.Execute(shrink_mask,bones)

    smoothingGaussFilter.SetSigma(
        calculate_gauss_sigma(1.0,voxDim)
    )
    smoothed = smoothingGaussFilter.Execute(no_bone)

    binFilter.SetLowerThreshold(_muscle_thresh[0])
    binFilter.SetUpperThreshold(_muscle_thresh[1])
    musc_init = binFilter.Execute(smoothed)

    binFilter.SetLowerThreshold(_fat_thresh[0])
    binFilter.SetUpperThreshold(_fat_thresh[1])
    fat_init = binFilter.Execute(smoothed)



    musc_seed = remove_small_islands(musc_init,_pixNum = 50)
    fat_seed = remove_small_islands(fat_init,_pixNum = 50)


    if _write_inter: 
        musc_seed_name = output_dir / f"{name}_Muscle_Seed.nii.gz"
        fat_seed_name = output_dir / f"{name}_Fat_Seed.nii.gz"

        write_image_with_date(
            musc_seed,
            str(musc_seed_name),
            CreationDate
        )

        write_image_with_date(
            fat_seed,
            str(fat_seed_name),
            CreationDate
        )


    # Check if array is empty
    if not np.any(sitk.GetArrayFromImage(musc_seed)):
        sys.exit("No muscle seed formed")

    if not np.any(sitk.GetArrayFromImage(fat_seed)):
        sys.exit("No fat seed formed")
    


    if _v: verbose_out("Iterative Mask Generation")
    max_iter = 10
    dilateFilter.SetKernelRadius(1)
    erodeFilter.SetKernelRadius(1)
    binFilter.SetUpperThreshold(255)
    binFilter.SetLowerThreshold(254)

    fat = fat_seed
    musc = musc_seed

    

    for i in range(max_iter):
        fat_dil = dilateFilter.Execute(fat)
        mus_dil = dilateFilter.Execute(musc)


        tmp_fat = subFilter.Execute(fat_dil,mus_dil)
        tmp_mus = subFilter.Execute(mus_dil,fat_dil)


        fat = binFilter.Execute(tmp_fat)
        musc = binFilter.Execute(tmp_mus)

    # Used for removal of these areas in final mask

    bone_inv = invFilter.Execute(bones)
    skin_inv = invFilter.Execute(skin_mask)



    

    if _v: verbose_out("Finalizing Masks")
    musc_clean = finalize_mask(musc, _rad = 10)
    fat_clean = finalize_mask(fat, _rad = 8)


    musc_tmp = multFilter.Execute(musc_clean,bone_inv)
    musc_final = multFilter.Execute(musc_tmp,skin_inv)

    fat_tmp = multFilter.Execute(fat_clean,bone_inv)
    fat_tmp2 = multFilter.Execute(fat_tmp,skin_inv)
    fat_select = andFilter.Execute(fat_tmp2,sitk.Cast(no_skin,sitk.sitkUInt8))
    closeFilter.SetKernelRadius(8)
    fat_final = closeFilter.Execute(fat_select)
    

    write_image_with_date(
        musc_final,
        str(output_dir / f"{name}_Muscle.nii.gz"),
        CreationDate
    )

    write_image_with_date(
        fat_final,
        str(output_dir / f"{name}_Fat.nii.gz"),
        CreationDate
    )




    if _v: verbose_out("Calculating Statistics")


    shapeStatsFilter.Execute(fat_final)
    fv = shapeStatsFilter.GetPhysicalSize(1)


    shapeStatsFilter.Execute(no_bone_mask)
    tv = shapeStatsFilter.GetPhysicalSize(1)
    print(tv)

    shapeStatsFilter.Execute(musc_final)
    mv = shapeStatsFilter.GetPhysicalSize(1)
    csa = mv / height

    # Ratio of Muscle Volume to total limb volume
    mv_tv = mv / tv 

    # Density
    labelStatsFilter.Execute(roi, musc_final)

    musc_p = labelStatsFilter.GetMean(1)
        
    ## TODO:: save to output csv file or something
    print("::::::::::::::::::::::::::::::::::::")
    print(":::::::::::::::Stats::::::::::::::::")
    print("::::::::::::::::::::::::::::::::::::")
    print(f"Muscle Volume: {mv:.2f} [mm\N{SUPERSCRIPT THREE}]")
    print(f"Muscle Cross Sectional Area: {csa:.2f} [mm\N{SUPERSCRIPT TWO}]")
    print(f"Fat Volume: {fv:.2f} [mm\N{SUPERSCRIPT THREE}]")
    # print(f"Subcutaneous Fat Volume: {satv:.2f} [mm\N{SUPERSCRIPT THREE}]")
    # print(f"Intramuscular Fat Volume: {imatv:.2f} [mm\N{SUPERSCRIPT THREE}]")
    print(f"MV/TV: {mv_tv:.2f} [%]")
    print(f"Muscle \N{GREEK SMALL LETTER RHO}: {musc_p:.2f} [mg HA/cm\N{SUPERSCRIPT THREE}]")
    print("::::::::::::::::::::::::::::::::::::")
    



    if _write_inter: 
        if _v: verbose_out("Writing Combined Labels")
        combined_labels = [
            skin_mask,
            fat_final,
            musc_final,
            bones
        ]

        clampFilter.SetUpperBound(
            len(combined_labels) + 1
        )


        img_tmp = sitk.Image(combined_labels[0].GetSize(),
                    combined_labels[0].GetPixelID())
        img_tmp.CopyInformation(combined_labels[0])

        for idx, label in enumerate(combined_labels, start=1):
            img_tmp += label * idx

        combined_masks = clampFilter.Execute(img_tmp)

        combined_labels_name = output_dir / f"{name}_Combined_Labels.nii.gz"
        write_image_with_date(
            combined_masks,
            str(combined_labels_name),
            CreationDate
        )


        return SoftTissueOutput(mv,musc_p,csa,tv,fv,mv_tv)
    
def validate(input_directory: str):
    """Run soft tissue segmentation batch for validation images.
    Write results to excel file. 

    :param input_directory: Input directory to segment.
    :type input_directory: str
    """
    input_directory = Path(input_directory)
    fat_thresh = [-599,-150]
    muscle_thresh = [50,600]

    output_directory = input_directory / "Output"
    output_directory.mkdir(exist_ok = True)
    results = list()

    for file in input_directory.iterdir():
        if file.suffix != ".isq": continue

        name = file.stem
        print(name)
        results_folder = output_directory / name
        output = SoftTissueSegmentation(
            file,
            results_folder, 
            fat_thresh,
            muscle_thresh,
            _write_inter = True,
            _preprocess_name = True
        )

        row = list()
        row.extend([name])
        row.extend(output.stats)

        results.append(row)
    
    df = pd.DataFrame(results, columns = ["Subject","Muscle.Volume",
                                 "Muscle.Density","CSA","Total.Volume","Fat.Volume",
                                 "MV.TV"])
    
    results_file = output_directory / "Results.xlsx"

    try:
        df.to_excel(str(results_file),"Results",index=False)
    
    except PermissionError:
        results_file = output_directory / "Results1.xlsx"

        try:
            df.to_excel(str(results_file),"Results",index=False)

        except PermissionError:
            csv = df.to_csv(index=False)
            retry_file = results_file.stem
            retry_file = retry_file.with_suffix(".csv")
            with open(retry_file,'w') as f:
                f.write(csv)


def parse_args():
    parser = argparse.ArgumentParser("Input for Soft Tissue Segmentation Program")

    parser.add_argument("--input","-i")
    parser.add_argument("--output","-o")
    parser.add_argument("--muscle_lower","-mL",type=int,default=50)
    parser.add_argument("--muscle_upper","-mU",type=int,default=600)
    parser.add_argument("--fat_lower","fL",type=int, default=-599)
    parser.add_argument("--fat_upper","fU",type=int, default=-150)
    parser.add_argument("--write_intermediate","-I",action="store_true")
    parser.add_argument("--preprocess_name","-p", action="store_true")
    parser.add_argument("--verbose","-v",action="store_true")

    return parser.parse_args()



if __name__ == "__main__":
    import argparse

    args = parse_args()
    SoftTissueSegmentation(Path(args.input),
                           Path(args.output),
                           _fat_thresh = [args.fat_lower, args.fat_upper],
                           _muscle_thresh = [args.muscle_lower, args.muscle_upper],
                           _write_inter = args.write_intermediate,
                           _preprocess_name = args.preprocess_name,
                           _v = args.verbose)