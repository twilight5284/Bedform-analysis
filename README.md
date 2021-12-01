# DPA: A Bedform Analysis Procedure
DPA (Dune Parameters Analysis) is presented as a free MATLAB software for analyzing bedform geometry parameters. 

## Method framework
This is an automated method, combining 2D Fourier analysis, wavelet transform, zero-crossing analysis and a variety of filters.
Firstly, the wavelength of interest can be automatically determined by a series of 2D Fourier analyses. Secondly, the dominant regional orientation of the surfaces is calculated by the 2D DFT (Cazenave et al., 2013). Thirdly, the matrix is rotated, re-gridded, and split into numerous profiles. Target bedform profiles are separated by wavelet transform and spline filters (Gutierrez et al., 2013). Then, dune crests and troughs are extracted by zero-crossing analysis. Finally, the individual geometric parameters such as wavelength, height, asymmetry, lee-slope angle and so on, are calculated. 

## Dune parameters definition 
The wavelength is defined as the distance between subsequent trough locations and height as the perpendicular distance between the crest and the straight line connecting the adjacent troughs. The asymmetry is defined as the difference of the distance between the crest and the trough west of the crest and the distance between the crest and the trough east of the crest divided by the wavelength. The average and maximum lee-slope angles are detected by the derivative of the lee-slope of the dune.
## Important
This program is a free software in the hope that it can be useful. And you can redistribute it and/or modify it. 

If you use it and want to quote it, here is the paper.
1.	Wang, L., Yu, Q., Zhang, Y., Flemming, B.W., Wang, Y., Gao, S., 2020. An automated procedure to calculate the morphological parameters of superimposed rhythmic bedforms. Earth Surface Processes and Landforms. 45, 3496–3509. https://doi.org/10.1002/esp.4983
2.	Wang, L., Yu, Q., Gao, S., 2019. A combined method to calculate superimposed 2-D dune morphological parameters. In: Lefebvre, A.; Garlan, T., and Winter, C. (eds), Proceedings of the Sixth International Conference on Marine and River Dune Dynamics (MARID VI) (Bremen, Germany), pp. 243-248. https://www.marum.de/Binaries/Binary18546/MARIDVI-Wang-Li.pdf

## Instructions
### 1.	Bathymetry data preprocessing

The input data must be a rectangular surface saved in a mat file. It can be a matrix named ‘data’ with three columns of x, y, and depth(z) (inputtype: 1, such as the example data ‘testdata_type_1.mat’) or three meshed matrixes named ‘x’, ‘y’, and ‘z’ (inputtype: 2, such as the example data ‘testdata_type_2.mat’). 

### 2.	Measurement of the wavelength of interest

The function LT can be used to calculate the wavelength of interest, and display the bathymetry map and a typical bedform profile perpendicular to the dune crest. After the calculation, the results will be displayed. 

Here is the example:

[LT, T, PHID] = LT('D:\DPA\testdata_type_2.mat', 2，1);
1
2
.
.
13 

There are 2 wavelength(s) of interest, they are 10 m, 208 m.

### 3.	Dune geometries analyses

The wavelength of interest calculated by the LT function is the input of the dune parameters calculation. Besides, the surface matrix, inputtype, data resolutions, project storage path, and folder name are also the inputs of the DPA function. After the calculation, the results will be saved in the project folder. 

There are two wavelengths of interest in the example data, so we need to run it twice.

Here is the example: 

First run: 10 m as the wavelength of interest.

DPA('D:\ DPA\', 'test1', 'D:\DPA\testdata_type_2.mat', 2, 10, 1);

The subset 1_1 is finished. (1/55)
.
.
.

The subset 5_11 is finished. (55/55)

Second run: 208 m as the wavelength of interest.

DPA('D:\ DPA\', 'test2', 'D:\DPA\testdata_type_2.mat', 2, 208, 1);

The subset 1_1 is finished. (1/1)

### 4.	Outputs

Outputs are saved in the output.mat file.
phiM is the subset dominant dune crest orientation.
BPALL is the individual dune parameters. They are x of dune crest, y of dune crest, wavelength, left wavelength, right wavelength, wave height, asymmetry, stepness, average lee slope angle, average lee slope angle between knick points, maximum lee slope angle, and depth of dune crest.

![image](https://user-images.githubusercontent.com/58336082/129526180-4c5c4820-f0cc-4c58-9e1b-e244dfcb688a.png)

## Example data
The example data is a 500 m × 1000 m rectangular domain located on a sandbank in the north of the Dover Strait, UK, called South Falls. The bathymetric data was downloaded from a website called UKHO INSPIRE Portal & Bathymetry DAC (http://aws2.caris.com/ukho/mapViewer/map.action), which is provided by the UK Hydrographic Office (www.ukho.gov.uk) under an open government license. 

## If you have any questions, please send an email to twilight528400@hotmail.com
