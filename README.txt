Optimization_4x_FLIM_Fig4.m is the main code for upsampling a low resolution FLIM data. The code first upsamples the photoncount datacube which is followed by lifetime evaluation based on least squares fitting. 
SlideWinSum.m, DownSampleBlur_2.m, TempoIntegration.m and decaymodelSingle.m are functions required therein.

CMOS_Raw.mat is the CMOS intensity from experiment. Spad_raw.mat is the low resolution Spad measurement.

One would need to change the path of the raw data files prior to running the code.

(written at University of Glasgow,UK)
