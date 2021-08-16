# DPA: A Bedform Analysis Procedure
DPA (Dune Parameters Analysis) is presented as a free MATLAB software for analyzing bedform geometry parameters. 

Method framework
This is an automated method, combining 2D Fourier analysis, wavelet transform, zero-crossing analysis and a variety of filters.
Firstly, the wavelength of interest can be automatically determined by a series of 2D Fourier analyses. Secondly, the dominant regional orientation of the surfaces is calculated by the 2D DFT (Cazenave et al., 2013). Thirdly, the matrix is rotated, re-gridded, and split into numerous profiles. Target bedform profiles are separated by wavelet transform and spline filters (Gutierrez et al., 2013). Then, dune crests and troughs are extracted by zero-crossing analysis. Finally, the individual geometric parameters such as wavelength, height, asymmetry, lee-slope angle and so on, are calculated. 

Dune parameters definition 
The wavelength is defined as the distance between subsequent trough locations and height as the perpendicular distance between the crest and the straight line connecting the adjacent troughs. The asymmetry is defined as the difference of the distance between the crest and the trough west of the crest and the distance between the crest and the trough east of the crest divided by the wavelength. The average and maximum lee-slope angles are detected by the derivative of the lee-slope of the dune.
Important
This program is a free software in the hope that it will be useful. And you can redistribute it and/or modify it. 

If you use it and want to quote it, here is the paper.
Wang, L.; Yu, Q.; Gao, S.; Zhang, Y., and Flemming BW. Submitted. An automated procedure to calculate the morphological parameters of superimposed two-dimensional dunes.

Instructions
1.	Bathymetry data preprocessing
The input data must be a rectangular surface saved in a mat file. It can be a matrix named ‘data’ with three columns of x, y, and depth(z) (inputtype: 1, such as the example data ‘testdata_type_1.mat’) or three meshed matrixes named ‘x’, ‘y’, and ‘z’ (inputtype: 2, such as the example data ‘testdata_type_2.mat’). 

2.	Measurement of the wavelength of interest
The function LT can be used to calculate the wavelength of interest, and display the bathymetry map and a typical bedform profile perpendicular to the dune crest. After the calculation, the results will be displayed. 
Here is the example:
[LT, T, PHID] = LT('D:\DPA\testdata_type_2.mat', 2);
1
2
.
  .
  13 
There are 2 wavelength(s) of interest, they are 9 m, 208 m.

3.	Dune geometries analyses
The wavelength of interest calculated by the LT function is the input of the dune parameters calculation. Besides, the surface matrix, inputtype, data resolutions, project storage path, and folder name are also the inputs of the DPA function. After the calculation, the results will be saved in the project folder. 

There are two wavelengths of interest in the example data, so we need to run it twice.
Here is the example: 
First run: 9 m as the wavelength of interest.
DPA('D:\ DPA\', 'test1', 'D:\DPA\testdata_type_2.mat', 2, 9, 1);
The subset 1_1 is finished. (1/55)
.
.
  .
The subset 5_11 is finished. (55/55)
Second run: 208 m as the wavelength of interest.
DPA('D:\ DPA\', 'test2', 'D:\DPA\testdata_type_2.mat', 2, 208, 1);
The subset 1_1 is finished. (1/1)

4.	Outputs
 
Outputs are saved in the output.mat file.
  
phiM is the subset dominant dune crest orientation.
 
BPALL is the individual dune parameters. They are x of dune crest, y of dune crest, wavelength, left wavelength, right wavelength, wave height, asymmetry, stepness, average lee slope angle, average lee slope angle between knick points, maximum lee slope angle, and depth of dune crest.

Example data
The example data is a 500 m × 1000 m rectangular domain located on a sandbank in the north of the Dover Strait, UK, called South Falls. The bathymetric data was downloaded from a website called UKHO INSPIRE Portal & Bathymetry DAC (http://aws2.caris.com/ukho/mapViewer/map.action), which is provided by the UK Hydrographic Office (www.ukho.gov.uk) under an open government license. 

Author
Li Wang - State Key Laboratory for Estuarine and Coastal Research, East China Normal University, Shanghai, China 
Qian Yu - Ministry of Education Key Laboratory for Coast and Island Development, Nanjing University, Nanjing, China
Shu Gao - State Key Laboratory for Estuarine and Coastal Research, East China Normal University, Shanghai, China

If you have any questions about code, please send an email to twilight528400@hotmail.com

References
Cazenave PW, Dix JK, Lambkin DO, McNeill LC. 2013. A method for semi-automated objective quantification of linear bedforms from multi-scale digital elevation models. Earth Surface Processes and Landforms 38(3): 221–236.
Gutierrez RR, Abad JD, Parsons DR, Best JL. 2013. Discrimination of bed form scales using robust spline filters and wavelet transforms: Methods and application to synthetic signals and bed forms of the Río Paraná, Argentina. Journal of Geophysical Research: Earth Surface 118(3): 1400–1418.
