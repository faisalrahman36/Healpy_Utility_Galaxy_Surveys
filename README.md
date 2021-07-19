# Healpy_Utility_Galaxy_Surveys

This repo contains some utility code to help with healpy based galaxy survey analysis.

The files have code for:

1-Function to convert right ascensions given in hour,minutes and seconds format into degrees.
2-Function to convert declination given in degrees,arcminutes and arcseconds format into degrees only format.
3-Function (pix2ang_radec) to do the reverse of ang2pix_radec.
4-Function to perform pixf2ang operation but returns values in theta and phi format not in ra and dec format.
5-Function to give counts of galaxies on each healpix pixel. 
6-Function to create density fluctuation maps with mask.
7-Function to convert Cls to CCF_theta or ACF_theta.
8-Function to calculate error bars on CCF_theta or ACF_theta.
9-Function to calculate shot noise from average source count and nside values.
10-Function to calculate fsky ratio using the survey mask. Survey mask should have value=0 for masked pixels and value=1 for valid survey pixels. 

You can read :
a)C. Blake,P . G .Ferreira,J . Borrill,Mon. Not. R. Astron. Soc.,351,923(2004) 
b)S.F. Rahman, CJP, 93(4): 384-394, 10.1139/cjp-2014-0339 (2015)
c)S.F. Rahman and M. J. Iqbal Number Counts, Confusion, Mapping Issues, and Sky Coverage Analysis for Radio Continuum Surveys through Emu Early Science, EMU-ASKAP, and WODAN Especially for Cosmology Science Goals. Astron. Rep. 63, 515–526 (2019). https://doi.org/10.1134/S1063772919070072 https://arxiv.org/abs/1612.08226 (2019)

If you find the code useful in your analysis the please cite the repo and papers (b) and (c) which used earlier version of this code and so were part of the development process.

Please feel free to send suggestions or queries at: faisalrahman36@hotmail.com

    
    
