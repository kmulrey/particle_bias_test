import os
import numpy as np
import cPickle as pickle
from optparse import OptionParser
import scipy.fftpack as fftp
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import signal
import ROOT
import time
import glob
import re
from pycrtools import crdatabase as crdb
import pycrtools as cr

nTrace=4000
nDet=20
dt=2.5
times=np.arange(0,4000*dt,dt)

nLasa=5

utc_day=60*60*24
utc_hour=60*60

jan2011=1293840000
jan2012=1325376000
jan2013=1356998400
jan2014=1388534400
jan2015=1420070400
jan2016=1451606400
jan2017=1483228800
jan2018=1514764800
jan2019=1546300800

data_dir='/vol/astro3/lofar/vhecr/lora_triggered/LORAraw/'

on_off_filename='/vol/astro3/lofar/sim/kmulrey/lora/detectors_on_off.txt'


def find_trigger(event_id):

    #print 'running event: {0}'.format(event_id)
    
    detector=np.zeros([nDet])
    ymd=np.zeros([nDet])
    gps=np.zeros([nDet])
    ctd=np.zeros([nDet])
    nsec=np.zeros([nDet])
    trigg_condition=-1*np.ones([nDet])
    trigg_pattern=np.zeros([nDet])
    total_counts=np.zeros([nDet])
    pulse_height=np.zeros([nDet])
    pulse_width=np.zeros([nDet])
    counts=np.zeros([nDet,nTrace])
    on_off=np.zeros([nDet])
    local_trigger=-1*np.ones([nLasa])
    trigger_setting=-1
    dbManager = crdb.CRDatabase("crdb", host="coma00.science.ru.nl",user="crdb", password="crdb", dbname="crdb")
    db = dbManager.db

    event = crdb.Event(db = db, id = event_id)
    lora_nsec=event["lora_nsecs"]
    lora_utc=event["lora_utc_time_secs"]
    #print 'lora utc: {0}'.format(lora_utc)
    
    
    
    data=np.genfromtxt(open(on_off_filename,'r'))
    logUTC=data.T[0]
    on_off_log=data.T[1:21].T
    
    ## find closest day to lora utc
    idx = (np.abs(logUTC - lora_utc)).argmin()
    if logUTC[idx]>lora_utc:
        idx=idx-1
    on_off= on_off_log[idx]

    ## find detectors that are on/off daily with
    
    
    
    
    
    
    
    
    
    


    year=time.gmtime(lora_utc).tm_year
    month=time.gmtime(lora_utc).tm_mon
    day=time.gmtime(lora_utc).tm_mday

    # find correct log file for event (daily, sometimes not with a standard name)
    log_list= glob.glob(data_dir+'{0:04d}{1:02d}{2:02d}*.log'.format(year,month,day))
    log_list.extend(glob.glob(data_dir+'{0:04d}{1:02d}{2:02d}*.log'.format(year,month,day-1)))

    if day==1:
        log_list.extend(glob.glob(data_dir+'{0:04d}{1:02d}{2:02d}*.log'.format(year,month-1,30)))
        log_list.extend(glob.glob(data_dir+'{0:04d}{1:02d}{2:02d}*.log'.format(year,month-1,31)))
    if day==1 and month==1:
        log_list.extend(glob.glob(data_dir+'{0:04d}{1:02d}{2:02d}*.log'.format(year-1,12,31)))
        log_list.extend(glob.glob(data_dir+'{0:04d}{1:02d}{2:02d}*.log'.format(year-1,12,30)))


    found_utc=0
    for l in np.arange(len(log_list)):
        if str(int(lora_utc)) in open(log_list[l]).read():
            log_file_name=log_list[l]
            root_file_name=log_file_name.split('.')[0]+'.root'
            found_utc=1

    if found_utc==1:

        log_file=open(log_file_name,'r')
        root_file=ROOT.TFile.Open(root_file_name)

        tree_sec = root_file.Get("Tree_sec")
        tree_event = root_file.Get("Tree_event")
        tree_log = root_file.Get("Tree_log")
        tree_noise = root_file.Get("Tree_noise")

        active_lasas,trigger_setting=find_active_stations(log_file)
        
        event_index=find_entry_number(lora_utc,lora_nsec,tree_event)


        for d in np.arange(nDet):
            detname='Det'+str(d+1)
            det=tree_event.GetBranch(detname)

            detector[d],ymd[d],gps[d],ctd[d],nsec[d],trigg_condition[d],trigg_pattern[d],total_counts[d],pulse_height[d],pulse_width[d],counts[d]=getData(det,event_index)
        
            if np.max(counts[d])>0.0:
                on_off[d]=1

        # to get the lasa local trigger condition
        #katie: question- how can only 2 detectors have condition 2/4??
        lasa1=np.min(trigg_condition[0:4])
        lasa2=np.min(trigg_condition[4:8])
        lasa3=np.min(trigg_condition[8:12])
        lasa4=np.min(trigg_condition[12:16])
        lasa5=np.min(trigg_condition[16:20])
        local_trigger=np.asarray([lasa1,lasa2,lasa3,lasa4,lasa5])

        ## change to stricter condition in case of broken files
        for t in np.arange(len(local_trigger)):
            if local_trigger[t]<0.5:
                local_trigger[t]=3.0
        print 'event okay'


    else:
        print 'didn\'t find matching file'




    ## if there is a problem reading local trigger setting, change everything to strictest condition (3/4)
    if (-1 in local_trigger)==True:
        local_trigger=3*np.ones([nLasa])

    ## check if station or detector trigger
    ## catch weird settings

    trigger_type='d'
    if (int(trigger_setting)<6.0 and int(trigger_setting)>0.0):
        print 'station trigger'
        trigger_type='s'
    elif (int(trigger_setting)>5.0 and int(trigger_setting)<21.0):
        print 'detector trigger'
        trigger_type='d'
    else:
        if lora_utc<(jan2012+6*30*utc_day):
            trigger_setting=5.0
            trigger_type='s'
        elif lora_utc>(jan2012+6*30*utc_day) and lora_utc<jan2013:
            trigger_setting=4.0
            trigger_type='s'
        else:
            trigger_setting=13
            trigger_type='d'


    return on_off,int(trigger_setting),local_trigger,trigger_type


def find_active_stations(file):
    
    trigger_string='LOFAR trigger settings'
    stations=['CS003','CS004','CS005','CS006','CS007']
    active_stations=np.zeros([nLasa])
    trigger_setting=0
    
    for line in file:
        for n in np.arange(nLasa):
            if re.match(stations[n], line):
                #print line.split()[1]
                active_stations[n]=int(line.split()[1])
    
        if re.match(trigger_string, line):
            trigger_setting= line.split()[len(line.split())-1]
    return active_stations,trigger_setting

def getTime(det, entry):
    det.GetEntry(entry)
    ymd=det.GetLeaf('YMD').GetValue()
    gps=det.GetLeaf('GPS_time_stamp').GetValue()
    ctd=det.GetLeaf('CTD').GetValue()
    nsec=det.GetLeaf('nsec').GetValue()
    return ymd,gps,ctd,nsec


def find_entry_number(lora_utc,lora_nsec,tree_event):
    event=-1
    
    det1=tree_event.GetBranch('Det1')
    det5=tree_event.GetBranch('Det5')
    det9=tree_event.GetBranch('Det9')
    det13=tree_event.GetBranch('Det13')
    det17=tree_event.GetBranch('Det17')
    
    det1.GetLeaf('GPS_time_stamp')
    nE= det1.GetEntries()
    #times=np.zeros([nLasa,nE])
    trigger_check=0
    diff_best=1e10
    
    for e in np.arange(nE):
        
        ymd1,gps1,ctd1,nsec1= getTime(det1,e)
        ymd2,gps2,ctd2,nsec2= getTime(det5,e)
        ymd3,gps3,ctd3,nsec3= getTime(det9,e)
        ymd4,gps4,ctd4,nsec4= getTime(det13,e)
        ymd5,gps5,ctd5,nsec5= getTime(det17,e)
        times=[gps1,gps2,gps3,gps4,gps5]
        times_ns=[nsec1,nsec2,nsec3,nsec4,nsec5]
        
        if lora_utc in times:
            diff=np.max(np.abs(lora_nsec-np.asarray(times_ns)[np.asarray(times_ns)>1]))
            if diff<10000:  # w/in 10 us to avoid mis-triggers
                if diff<diff_best:
                    diff_best=diff
                    trigger_check=1
                    event=e



    return event

def getData(det, entry):
    
    det.GetEntry(entry)
    
    detector=det.GetLeaf('detector').GetValue()
    ymd=det.GetLeaf('YMD').GetValue()
    gps=det.GetLeaf('GPS_time_stamp').GetValue()
    ctd=det.GetLeaf('CTD').GetValue()
    nsec=det.GetLeaf('nsec').GetValue()
    trigg_condition=det.GetLeaf('Trigg_condition').GetValue()
    try:
        trigg_pattern=det.GetLeaf('Trigg_pattern').GetValue()
    except:
        trigg_pattern=-1
    total_counts=det.GetLeaf('Total_counts').GetValue()
    pulse_height=det.GetLeaf('Pulse_height').GetValue()
    pulse_width=det.GetLeaf('Pulse_width').GetValue()
    counts=det.GetLeaf('counts')
    hold=np.zeros([nTrace])
    for i in np.arange(nTrace):
        hold[i]=counts.GetValue(i)
    
    
    return detector,ymd,gps,ctd,nsec,trigg_condition,trigg_pattern,total_counts,pulse_height,pulse_width,hold



def find_trigger_from_file(event, dir):


    file=infile=open(dir+str(event)+'.dat')
    nLa=np.genfromtxt(infile,max_rows=1)
    nLa2=np.genfromtxt(infile,max_rows=1)
    nLa3=np.genfromtxt(infile,)



    if len(nLa)>0:
        trigger_cond=nLa2
        nActive=nLa3

    else:
        trigger_cond=0
        nActive=np.zeros([nDet])
    return nLa3,trigger_cond







