'''
Routine to check for particle detection bias.  Given coreas simulation directory, and event direction and core information, this will return 1 for a bias free event and 0 for a biased event.  
For each simulation, the particle desity at each detector is determined from the geant .lora files.  The detectors are filled with a particle density nTries times (100 default) using a Poisson distribution, and then the probability that the event triggers is calculated.  If this probability is over 'bias_condition' (here 95% as a default), then the simulation is considered bias free.  If all the simulations for an event meet this criteria, then the event is considered bias free.  (Note: as it is now, I am considering (N simulations -1)  passing the bias criteria for the event to be good.  This is because I have seen a few cases where one simulation is flawed for an otherwise good event.  To be studied further.)
    
'''


import os
import numpy as np
import matplotlib.pyplot as plt
import re

em_peak=6.7  # this converts energy deposited into equivalent (all sky) muons, which corresponds to the database partice density values.  The trigger conditions are roughly around this number, depending on background noise levels, so this says a detector is triggered if the particle density > 1 equivalent muon.
#nTries=500  # number of times to assign particle density for a given simulation
nDet=20 # number of detectors
nLasas=5 # number of stations

#trigger_condition=13    # number of LORA detectors required for trigger
#bias_condition=95.0    # percent detection rate to clear bias cut (for one event)
detectorpath='/vol/astro3/lofar/sim/kmulrey/Detector_Cord.dat'

def theta_phi(theta,phi,psi,x0,y0,z0):  # transform X0,y0,z0 to shower plane
    
    x1=x0*np.cos(phi)+y0*np.sin(phi)
    y1=-1*x0*np.sin(phi)+y0*np.cos(phi)
    z1=z0
    
    x2=x1*np.cos(theta)-z1*np.sin(theta)
    y2=y1
    z2=x1*np.sin(theta)+z1*np.cos(theta)
    
    x=x2*np.cos(psi)+y2*np.sin(psi)
    y=-1*x2*np.sin(psi)+y2*np.cos(psi)
    z=z2
    
    return x,y,z


def find_bias(directory,theta,phi,xcore,ycore, ntrials, working_detectors,trigger_type,local_trigger_condition=np.asarray([3.0,3.0,3.0,3.0,3.0]),min_percentage_trigger=95.0, trigger_condition_nof_detectors=13,verbose=False):
    '''
    - directory is the path to the coreas files, for example '/vol/astro3/lofar/sim/pipeline/events/175494398/0/coreas'
    - theta is the zenith angle of the event in radians
    - phi is the aziuthal angle of the event in radians
    - xcore amd ycore are shower core positions in meters
    '''
    

    psi=2*np.pi-phi
    iron_directory=directory+'iron/'
    proton_directory=directory+'proton/'
    Aeff=0.9*np.cos(theta)
    pass_cut=0  # return value, 1 = no bias, 0 = bias

    # collect geant files for both iron and proton
    
    files=[]
  
    for filename in os.listdir(iron_directory):
        if filename.endswith('.lora'):
            pathname=iron_directory+filename
            files.append(pathname)
        
    for filename in os.listdir(proton_directory):
        if filename.endswith('.lora'):
            pathname=proton_directory+filename
            files.append(pathname)



    # get LORA detector positions and put them in shower plane, relative to shower core

    detFile=open(detectorpath)
    #detFile=open('/Users/kmulrey/LOFAR/LORA/sims/Detector_Cord.dat')
    LORAdetectors=np.genfromtxt(detFile,comments="//",usecols=(1,2,3))
    #LORAdetectorsUse=LORAdetectors[working_detectors==1]

    #nUse=len(LORAdetectors)

    xdet, ydet, zdet = theta_phi(theta,phi,psi,LORAdetectors.T[0]-xcore,LORAdetectors.T[1]-ycore,LORAdetectors.T[2])

    #find radius of detectors from shower core and find bin associated with geant file binning
    rad=np.sqrt(xdet*xdet+ydet*ydet)
    radBin=(rad/5.0).astype(int)

    percent_hits_all=np.zeros(len(files)) # pecrentage of events (ntrials) that pass trigger condition for each .lora file
 
    for i in np.arange(len(files)):
        data=np.genfromtxt(files[i])
        runnr = files[i].split('on/')[1].split('T')[1].split('.')[0]
        
        EM=np.zeros([ntrials, nDet])  # equivalent muons (all sky) in each detector
        hits=np.zeros([ntrials])  # number of dectors hit
        LASA_triggered=np.zeros([nLasas])

        #EM = np.zeros(nDet)
        #data.T[1][radBin[j]]/em_peak is the geant deposit at the positition of the detector in 'all sky equivalent muon' units
        # Check for too high distances to detectors
        if max(radBin) > len(data.T[1]):
            print 'WARNING: one or more detectors too far from shower core. Possibly due to outlying core pos. or outlying LOFAR stations in detection.'
        
        for k in np.arange(ntrials):
            for j in np.arange(nDet):
                if working_detectors[j]==1:  # now including this to account for detectors that are not working
                    if radBin[j] < len(data.T[1]):
                        EM[k][j] = np.random.poisson(Aeff*data.T[1][radBin[j]]/em_peak)  # select EM deposit from poisson distribution
                    else:
                        EM[k][j] = 0.0
    
                            #poisson_lambda = Aeff*data.T[1][radBin[j]]/em_peak
            #EM[j] = 1 - np.exp( - poisson_lambda) # probability of >= 1 particle hits


            ### account for working detectors:

            EM[k]=EM[k]*working_detectors

            ### divide into LASAs
            LASAs=np.zeros([nLasas,4])

            LASAs[0]=EM[k][0:4]
            LASAs[1]=EM[k][4:8]
            LASAs[2]=EM[k][8:12]
            LASAs[3]=EM[k][12:16]
            LASAs[4]=EM[k][16:20]

            ########## check station level trigger
            
            for l in np.arange(nLasas):
                if len(LASAs[l][LASAs[l]>=1])>=local_trigger_condition[l]:
                    LASA_triggered[l]=1
                        #print LASA_triggered
            ########## if station trigger
            if trigger_type=='s':
                hits[k]=len(LASA_triggered[LASA_triggered==1]) # number of detectors triggered for each event (trigger condition > 1 EM)-> also what is used for LORA file threshold
    
    
            ########## if detector trigger

            if trigger_type=='d':

                hits[k]=len(EM[k][EM[k]>=1]) # number of detectors triggered for each event (trigger condition > 1 EM)-> also what is used for LORA file threshold

        percent_hits = 100.0 * float(len(hits[hits>=trigger_condition_nof_detectors])) / ntrials # chance of triggering for this simulation
        percent_hits_all[i] = percent_hits # save trigger info for this simulation

    #import pdb; pdb.set_trace()
    pass_bias=percent_hits_all[percent_hits_all > min_percentage_trigger]  # number of simulations that pass bias
    min_percentage_hit = np.min(percent_hits_all)
    print trigger_condition_nof_detectors
    if verbose:
        print 'number of total events: {0}'.format(len(percent_hits_all))
        print 'number of events that pass bias condition: {0}'.format(len(pass_bias))
        print 'Minimum percentage of hit: %2.4f' % min_percentage_hit
        if min_percentage_hit < 100.0:
            print 'File with minimum percentage hits: %s' % files[np.argmin(percent_hits_all)]
    '''
    - to assign bias to event, check that the all the simulations (-1) pass the bias cut, for example, that 95% of the time a given simulation will trigger
    - there is a -1 there in len(percent_hits_all-1) because sometimes one simulation would not trigger for an otherwise clearly good event-> check further
    '''


    if len(pass_bias)<len(percent_hits_all-1):
        pass_cut = False

    if (len(pass_bias)==len(percent_hits_all) or len(pass_bias)==(len(percent_hits_all)-1)):
        pass_cut = True


    return (pass_cut, min_percentage_hit, percent_hits_all)  # 1 = passed bias cut, 0 = biased event





