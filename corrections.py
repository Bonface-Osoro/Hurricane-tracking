# -*- coding: utf-8 -*-


def resolutioncorr(satdatach2lath, satdatach2lonh, satdatach2radh):
    import numpy as np
    satdatach2lat=np.zeros((700,1000)) #<-- Equal the size of the arrays for channel 1 and 3
    satdatach2lon=np.zeros((700,1000))
    satdatach2rad=np.zeros((700,1000))        
    for i in range(0,699):
        for j in range(0,999):
           satdatach2lat[i,j]=np.mean(satdatach2lath[2*i:2*(i+1),2*j:2*(j+1)])
           satdatach2lon[i,j]=np.mean(satdatach2lonh[2*i:2*(i+1),2*j:2*(j+1)])
           satdatach2rad[i,j]=np.mean(satdatach2radh[2*i:2*(i+1),2*j:2*(j+1)])
    return satdatach2lat, satdatach2lon, satdatach2rad
    

def rayleighcorr(k1,satdatach1lat,satdatach1lon,satdatach2lat,satdatach2lon,satdatach1rad,satdatach2rad,satdatach3rad):

    from datetime import datetime, timedelta, timezone
    from pyspectral.rayleigh import Rayleigh
    from pysolar.solar import get_azimuth_fast
    import numpy as np
    
    # Atmospheric corrections for channel 1 and 2
    abi=Rayleigh('GOES-16', 'abi') #Loading specific satellite instrument characteristics for rayleigh correction
    doy=float(k1[31:34])-1.
    hh=float(k1[34:36])
    dt = datetime(2017,1,1,0,0,0,363238,tzinfo=timezone.utc)
    dtdelta = timedelta(days=doy,hours= hh)
    date=dtdelta+dt
    
    sunz=get_azimuth_fast(satdatach1lat,satdatach1lon,date) #Calculate local sun zenith angle for rayleigh correction for channel 1
    
    #Calculate local satellite zenith angle for rayleigh correction for channel 1
    teil1=np.cos(satdatach1lat[:,:])
    teil2=np.cos(-75.2+satdatach1lon[:,:])
    teil3=teil1*teil2
    satz=np.arccos(teil3)

    
    #Calculate angle difference for both channels
    ssadiff=satz-sunz
    
    #Calculate Rayleigh Correction for each pixel for both channels
    refl_cor_m1 = abi.get_reflectance(sunz, satz, ssadiff, 'ch1')
    refl_cor_m2 = abi.get_reflectance(sunz, satz, ssadiff, 'ch2')
    
    #Substract the correction from the corresponding channel data
    RRC=satdatach2rad-refl_cor_m2 #Channel 2 wavelength corresponds to Red in the optic spectrum
    BRC=satdatach1rad-refl_cor_m1 #Channel 1 wavelength corresponds to Blue in the optic spectrum
    
    #RRC=satdatach2rad
    #BRC=satdatach1rad
    
    #There is no Channel that corresponds to the green portion of the optic spectrum, but
    #we can calculate it from the Channels 1, 2, and 3 as follows:
    GRC=BRC*0.45+RRC*0.45+0.1*(satdatach3rad) 
   
    
    #Sorting out values below 0
    for i in range(0,699):
        for j in range(0,999):
            if RRC[i,j]<0:
                RRC[i,j]=0.
            if BRC[i,j]<0:
                BRC[i,j]=0.
            if GRC[i,j]<0:
                GRC[i,j]=0.
    
    return RRC, BRC, GRC    
