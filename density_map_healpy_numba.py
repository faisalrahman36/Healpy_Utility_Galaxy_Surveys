import numpy as np
import numba as nb
import datetime
import healpy as hp
import pylab as pl
#This implementation may be useful for large density maps generation. This implementation uses numba jit utility.
#Suitable for analysis involving Cls from galaxy density fluctuations healpix maps.
#for Healpix resolution parameter 
nside=8
#We do not have return parameters in this implementation so we will need to define  and initialize output parameters here:
survey_pix=0.
galaxy_count=0.
masked_pix=0.
average_count=0.
density_fluctuation_map=np.zeros(hp.nside2npix(nside),dtype=float)
#Provide mask and Healpix number count map
'''
:param inp_map: Healpix Number count map
:param mask: Healpix map mask
'''
map1=np.random.randint(10,size=hp.nside2npix(nside))*1. #inp_map:assign  actual count map
mask=np.random.randint(2,size=hp.nside2npix(nside))*1. #assign  actual mask

#In jit case we also have to write output parameters with input
@nb.jit('void(float64[:],float64[:],float64[:],float64,float64,float64,float64)', nopython=True)
def masked_density_map_numba(inp_map,mask,density_fluctuation_map,survey_pix,galaxy_count,masked_pix,average_count):
    '''
    Create density fluctuation maps and combine with a mask file

    :param inp_map: Healpix Number count map
    :param mask: Healpix map mask
    :param nside: Nside used in creating both the healpix map and mask
    :parameters for output: Masked galaxy density map
    
    The function will not return anything but will modify the input parameters (arrays or scalars)
    which should be initialized as global.
    '''

    #print('inp_map',inp_map)
    masked_pix=0
    survey_pix=0
    galaxy_count=0
    #density_fluctuation_map=np.zeros(n,dtype=float)



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


    
    
 
    



'''
Numba jit function call and computation time calculation
'''
   
c = datetime.datetime.now()
masked_density_map_numba(map1,mask,density_fluctuation_map,survey_pix,galaxy_count,masked_pix,average_count)

d = datetime.datetime.now()

print('numba')
delta = d - c
print(delta)
print(float(delta.total_seconds() * 1000),'ms') # milliseconds


print(c,'c')
print(d,'d')
print('map',density_fluctuation_map)
#Testing maps with mollview
hp.mollview(density_fluctuation_map)
#Testing maps with anafast
print(hp.sphtfunc.anafast(density_fluctuation_map,lmax=5))
pl.show()
