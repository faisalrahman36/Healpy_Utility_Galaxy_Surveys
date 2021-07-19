# -*- coding: utf-8 -*-
__author__ = 'Syed Faisal ur Rahman '
__version__ = "1.0.1"
__email__ = "faisalrahman36@hotmail.com"
__status__ = "Production"

#please feel free to email me if you find bugs or if you have suggestions for new functionality or improvement in methods


import numpy as np
import sympy.mpmath as smp
import healpy as hp


RadPerDeg = np.pi / 180.0 #Radians per degree multiplier factor
nside=32 #Default


'''
function to convert right ascensions given in hour,minutes and seconds format into degrees.
The function take three inputs h(hour), m(minute) and s(seconds).

To get correct results use hours between 0 and 24, minutes between 0 and 60 and seconds between 0 and 60.

'''
def ra_hms2deg(h,m,s):
    
    ra_deg=0.
    ra_deg= h*15. + m/4. + s/240.
    return ra_deg
    
    
'''
function to convert declination given in degrees,arcminutes and arcseconds format into degrees only format.
The function take three inputs d(degrees), am(arcminute) and asec(arcseconds).

To get correct results use degrees between -90 and +90, arcminutes between 0 and 60 and arcseconds between 0 and 60.

'''
def dec_damas2deg(d,am,asec):
    
    dec_deg=0.
    dec_deg= d + am/60. + asec/3600.
    return dec_deg
    
           


#It takes ra and dec information and converts into healpix or healpy grid. 
#ra and dec should be in degrees format
def ang2pix_radec(nside,ra, dec):
            phi = RadPerDeg * ra
            theta = (np.pi / 2.0 )- (RadPerDeg * dec)
    #For coordinate transformation

            theta_gal, phi_gal = hp.Rotator(coord=['C','G'])(theta, phi)
            ipix = hp.ang2pix(nside, theta_gal, phi_gal,nest=False)
            return ipix
            
#pix2ang_radec does the reverse of ang2pix_radec
#this will not return exact ra and dec but will return the centre of pixel ra and dec values for the given healpix pixel
def pix2ang_radec(nside, ipix, ra_rotate=0):
         (theta_gal, phi_gal) = hp.pix2ang(nside, ipix, nest=False)

    #For coordinate transformation
         theta, phi = hp.Rotator(coord=['G','C'])(theta_gal, phi_gal)
         ra = phi / RadPerDeg
         dec = ((np.pi / 2.0 ) - theta)/RadPerDeg
         if (ra_rotate):
                 if (ra >= 360 + ra_rotate):
                     ra = -360 + ra
                 elif (ra < 0 + ra_rotate):
                     ra = 360 + ra
         return (ra, dec)

#performs pixf2ang operation but returns values in theta and phi format not in ra and dec format
def pix2ang_theta_phi(nside,ipix):
             (theta_gal, phi_gal) = hp.pix2ang(nside, ipix, nest=False)

             return (theta_gal, phi_gal)

#The python function below, counts galaxies on each healpix pixel. 
#It is different from galaxy density maps . we will just need to subtract avg density and divide by 
#avg density after that. You can calculate avg density using the count map by something like total count/total pix in survey.

def ang2pix_count_radec(nside,ra,dec):
        ra=ra
        dec=dec

        map=np.zeros(hp.nside2npix(nside))

        for x in range(0,len(ra),1):
           tpix=ang2pix_radec(nside,ra[x],dec[x])
           if map[tpix]==0:
               map[tpix]=1
           else:
               map[tpix]=map[tpix]+1
        return map




#Code to create density fluctuation maps with mask
#mask file format:shoult contain value=1 for  surevey area and value=0 for masked regions


def masked_density_map(inp_map,mask,nside):
    #size of inp_map, mask and density_fluctuation_map should be same
    #inp_map: the galaxy count map generated from the ang2pix_count function
    #mask: the survey mask map: it should be a combined mask of all masks
    #combined mask can be obtained by multiplying two or more masks
    #nside : healpy.pixelfunc.nside2npix(nside)  function can be used to check if the length
    #of inp_map, mask and density_fluctution_map are same.

    masked_pix=0
    survey_pix=0
    galaxy_count=0
    average_count=0.
    density_fluctuation_map=np.zeros(hp.nside2npix(nside),dtype=float)



    for x in range(0,len(density_fluctuation_map),1):
       if mask[x]==0:
           masked_pix=masked_pix+1
       else:
           survey_pix=survey_pix+1
           #adding all the galaxy counts of survey pixels
           galaxy_count=galaxy_count + inp_map[x]



    average_count=float(galaxy_count)/survey_pix

    for y in range(0,len(density_fluctuation_map),1):
       if mask[y]==0:
           density_fluctuation_map[y]= -1.6375e+30
       else:
           density_fluctuation_map[y]=(inp_map[y]-average_count)/average_count


    return density_fluctuation_map,survey_pix,galaxy_count,masked_pix,average_count

'''
Healpy anafast function only produces Cls and not correlation function values. 
We can use the function below to convert Cls to CCF_theta or ACF_theta
You can read :
a)C. Blake,P . G .Ferreira,J . Borrill,Mon. Not. R. Astron. Soc.,351,923(2004) 
b)S.F. Rahman, CJP, 93(4): 384-394, 10.1139/cjp-2014-0339 (2015)
c)A. Raccanelli, et al., Mon. Not. R. Astron. Soc.,386, 2161(2008)
d)Hernández-Monteagudo, CAstronomy and Astrophysics, Volume 520, id.A101, 16 pp
e)Giannantonio, Tommaso et al. , Physical Review D, Volume 89, Issue 2, id.023511 

To understand the mathematics behind the routines.


'''


def CF_theta(cln,l_start=3,end_theta=10.,theta_samples=50):
    #print cln
    #theta_samples: number of samples between 0 and end_theta
    #more samples will require more computations
    #you can use the Shot_Noise function to calculate the shotnoise for a galaxy survey helpix map
    #CF_theta can be used for both theoretical and observed correlation functions 
      

    cl_theta=np.zeros(theta_samples, dtype =float)
    ctheta=np.linspace(0.,end_theta,num=theta_samples)

    for x in range(0,len(cl_theta),1):
        #cl_arb is used in the loop to keep the summation values
        cl_arb=0.
        #Will run through the starting multipole value to the last
        for y in range(l_start,len(cln),1):
          #Legendre Polynomial
            Pl=smp.legenp(y,0,np.cos(np.radians(ctheta[x])))
            arb=((2.*y+1.)/(4*np.pi))*cln[y]*Pl
            cl_arb+= arb

        cl_theta[x]=cl_arb
    #return theta values used and respective CCF_theta values
    return ctheta,cl_theta

#For error bars on CCF_theta or ACF_theta
def e_CF_theta(e_cln,l_start=3,end_theta=10.,theta_samples=50):
    #e_cln is variance in the cl values
    #theta_samples: number of samples between 0 and end_theta
    #more samples will require more computations

    e_cl_theta=np.zeros(theta_samples, dtype =float)
    ctheta=np.linspace(0.,end_theta,num=theta_samples)
    for x in range(0,len(e_cl_theta),1):
        #cl_arb is used in the loop to keep the summation values
        cl_arb=0.
        #Will run through the starting multipole value to the last
        for y in range(l_start,len(e_cln),1):
            #Legendre Polynomial
            Pl=smp.legenp(y,0,np.cos(np.radians(ctheta[x])))
            arb=(((2.*y+1.)/(4*np.pi))*e_cln[y]*Pl)**2
            cl_arb+= arb
        e_cl_theta[x]=cl_arb**0.5 #square root
    #return theta values used and respective error  bar values for CCF_theta
    return ctheta,e_cl_theta
    
    
'''
Function to calculate shot noise from average source count and nside values

'''

def Shot_Noise(average_count,nside=32):
    
    #average_count: number of sources per pixel
    pixelarea=0.
    conversionfactor=0.
    count_persteradian=0.
    ShotNoise=0.
    
    #number of squre degrees per steradian = 3282.810
    sq_deg_persteradian=3282.810
    #First step is to convert nside into square degrees pixel area
    
    pixelarea=hp.nside2pixarea(nside, degrees = True)
    
    #now calculate the multiplying factor to conver square degrees source count into steradian source count
    conversionfactor=sq_deg_persteradian/pixelarea 
    
    count_persteradian=average_count * conversionfactor
    
    #Shot Noise=1/Number of sources per steradian
    
    ShotNoise=1./count_persteradian
    
    return ShotNoise
    
'''
Function to calculate fsky ratio using the survey mask. Survey mask should have value=0 for masked pixels and value=1 for valid survey pixels. 

You can also use survey_pix and masked_pix from masked_density_map function to directly calculate the ratio

fsky=survey_pix/(survey_pix+masked_pix)
'''
def fsky(mask_map):
    
    survey_pix=0.
    masked_pix=0.
    
    for x in range(0,len(mask_map),1):
       if mask_map[x]==0:
           masked_pix=masked_pix+1
       else:
           survey_pix=survey_pix+1
           
    fsky=survey_pix/(survey_pix+masked_pix)
    
    return fsky
    
    