
# An implementation of the method ALTIS for the automated segmentation of each lung and trachea in CT scans

An implementation of the method ALTIS for the automated segmentation of each lung and trachea in CT scans is provided, following its description presented in the paper **"ALTIS: A fast and automatic lung and trachea CT-image segmentation method"** published in the journal of Medical Physics **[1]**, with some improvements as described in the next paragraph.

The original ALTIS algorithm sometimes presents problems with images when the entire stretcher appears in the image. In some of these cases, the stretcher is segmented instead of the respiratory system. In particular, the solution proposed in **[1]** was not working properly in the LCTSC database (Lung CT Segmentation Challenge) **[2]**. Our ALTIS implementation is adapted to automatically filter the stretcher whenever it is contained in the image.

We also provide an alternative version, that computes improved segmentation of the lungs with seeds by ALTIS **[1]** and final delineation based on ROIFT (Relaxed Oriented Image Foresting Transform) **[3]**, herein referred to as ALTIS + ROIFT.

**[1]** Sousa AM, Martins SB, Falc&atilde;o AX, Reis F, Bagatin E, Irion K,
Med Phys. 2019 Nov; 46(11):4970-4982; doi: [10.1002/mp.13773](https://doi.org/10.1002/mp.13773).

**[2]** Yang, J., Sharp, G., Veeraraghavan, H., Van Elmpt, W., Dekker, A., Lustberg, T., & Gooding, M. (2017). Data from Lung CT Segmentation Challenge (LCTSC) (Version 3) [Data set]. The Cancer Imaging Archive. doi: [10.7937/K9/TCIA.2017.3R3FVZ08](https://doi.org/10.7937/K9/TCIA.2017.3R3FVZ08).

**[3]** Caio L. Demario, Paulo A.V. Miranda, RELAXED ORIENTED IMAGE FORESTING TRANSFORM FOR SEEDED IMAGE SEGMENTATION, 26th IEEE International Conference on Image Processing (ICIP). Sep 2019; Taipei, Taiwan, pp. 1520-1524; doi: [10.1109/ICIP.2019.8803080](http://dx.doi.org/10.1109/ICIP.2019.8803080).


## Method

ALTIS consists of a sequence of image foresting transforms (IFTs) organized in three main steps:

1. lung-and-trachea extraction,
2. seed estimation inside background, trachea, left lung, and right lung, and
3. their delineation such that each object is defined by an optimum-path forest rooted at its internal seeds.

ALTIS + ROIFT changes the last delineation step to a method based on reference **[3]**.

## Referencing and citing

If you are working with this code to your project, in addition to referencing **[1]**, please also cite:

> Choi J., Condori M. A. T., Miranda P. A. V. and Tsuzuki M. S. G. Lung automatic seeding and segmentation: a robust method based on relaxed oriented image foresting transform.

## Authors

This particular code was implemented by:

- Jungeui Choi
- Marcos A. T. Condori
- Paulo A.V. Miranda
- Marcos S. G. Tsuzuki

## Source code

The source code was implemented in C/C++ language, compiled with gcc 9.4.0, and tested on a Linux operating system (Ubuntu 20.04.5 LTS 64-bit), running on an Intel® Core™ i5-10210U CPU @ 1.60GHz × 8 machine. 
The code natively supports volumes in the NIfTI format.

To compile the program, enter the folder and type **"make"**.
If you get the error **"fatal error: zlib.h: No such file or directory"**, then you have to install the zlib package: zlib1g-dev.

To segment a volume, you must run the **"altis"** executable inside the folder.
This implementation assumes that the patient orientation in the CT scan is from inferior to superior along the axial slices (z-axis), from right to left along the sagittal slices (x-axis), and from posterior to anterior along the coronal slices (y-axis).

### usage of ALTIS:

```
altis <volume> <output_type> [T=value] [left=file] [right=file]
	output_type......... 0 - segmentation labels,
	                     1 - seeds only.
Optional parameters:
	T................... threshold integer value
	                     (if not specified Otsu is used).
	left................ ground truth for left lung.
	right............... ground truth for right lung.
```


As output, the program generates the label image of the resulting segmentation in file **"segm_altis.nii.gz"** in the **"out"** subfolder, when **output_type** is zero.

### Program execution examples of ALTIS:

#### To execute the ALTIS method:

The following command computes the segmentation by ALTIS for the volume in the hypothetical file **"example01.nii.gz"**.

```
./altis example01.nii.gz 0
```

The following command computes the segmentation by ALTIS for the volume in the hypothetical file **"example01.nii.gz"**, using a fixed threshold of 200 on the residual image, instead of the default threshold defined as a percentage above Otsu's threshold.

```
./altis example01.nii.gz 0 T=200
```


### usage of ALTIS + ROIFT:

```
usage:
altis_roift <volume> [T=value] [left=file] [right=file]
Optional parameters:
	T................... threshold integer value
	                     (if not specified Otsu is used).
	left................ ground truth for left lung.
	right............... ground truth for right lung.
```

As output, the program generates improved segmentations of the left and right lungs, respectively, in the files **"segm_left_lung.nii.gz"** and **"segm_right_lung.nii.gz"** in the **"out"** subfolder.

### Program execution examples of ALTIS + ROIFT:

#### To execute seed generation by ALTIS and delineation by ROIFT:


The following command computes the segmentation of the lungs with seeds by ALTIS and delineation by ROIFT for the volume in the hypothetical file **"example01.nii.gz"**.

```
./altis_roift example01.nii.gz
```


The following command computes the segmentation of the lungs with seeds by ALTIS and delineation by ROIFT for the volume in the hypothetical file **"example01.nii.gz"**, using a fixed threshold of 200 on the residual image, instead of the default threshold defined as a percentage above Otsu's threshold.

```
./altis_roift example01.nii.gz T=200
```


## Contact

If you have any doubts, questions or suggestions to improve this code, please contact me at:
**pmiranda@ime.usp.br**

