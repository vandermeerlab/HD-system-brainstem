# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:47:16 2024

@author: plachanc
"""

import h5py
import scipy
import numpy as np
import bisect


def load_spike_timestamps(t_file):
    ''' load the .t file and return timestamps '''

    f = open(t_file, 'rb')
    #how long is the header
    header_string = f.read()
    header_length = header_string.index(b'%%ENDHEADER') + 10
    
    
    nvt_dtype = np.dtype([ 
        ('timestamps'  , '>u8')]) 
    #memmap the file
    mmap = np.memmap(t_file, dtype=nvt_dtype, mode='r+', 
        offset=(header_length+2))
    
    return np.array(mmap['timestamps']).flatten()
    

def load_video_timestamps(smi_file):
    ''' load timestamps from the SMI file '''
    
    with open(smi_file, 'r') as f:
        timestamps = []
        start = False
        for line in f:
            if line.startswith("<SYNC Start=0>"):
                start = True
            elif line.startswith("</BODY>"):
                start = False
            if start:
                line = line.strip()
                timestamps.append(int((line.split('ENUSCC>'))[1].split('</SYNC>')[0]))
                
    return timestamps


def get_hds(orientation_file):
    ''' compute head directions based on orientation data --
    modified from https://github.com/vandermeerlab/HD-system-brainstem/blob/master/__in%20process/GetOrientationValues.m '''
    
    try:
        with h5py.File(orientation_file, 'r') as file:
            hd_voltages = np.array(file['csc_tsd']['data']).flatten()
            hd_timestamps = np.array(file['csc_tsd']['tvec']).flatten()
            
    except:
        with h5py.File(orientation_file, 'r') as file:
            hd_voltages = np.array(file['orientation_tsd']['data']).flatten()
            hd_timestamps = np.array(file['orientation_tsd']['tvec']).flatten() 
    
    #borrowed from GetOrientationValues.m
    rangetouse = 360
    baseline = hd_voltages[0]
    maxR = np.max(hd_voltages)
    maxL = np.min(hd_voltages)
    Rrange = np.abs(maxR-baseline)
    Lrange = np.abs(maxL-baseline)
    rangediff = Lrange/Rrange
    fullrange = np.abs(maxL - maxR)
    
    subtracted_voltage = hd_voltages - baseline
    divisionconstant = fullrange/rangetouse
    hds = -subtracted_voltage/divisionconstant
    
    return hds[::10], hd_timestamps[::10]


def grab_nev_data(filename): 
    ''' get timestamps, ttl ids, and ttl messages from nev file '''
    
    #read file
    f = open(filename, 'rb')
    #skip past header
    f.seek(2 ** 14)
    #specity data types
    dt = np.dtype([('filler1', '<h', 3), ('time', '<Q'), ('id', '<h'),
                   ('nttl', '<h'), ('filler2', '<h', 3), ('extra', '<i', 8),
                   ('estr', np.dtype('a128'))])
    #grab the data
    temp = np.fromfile(f, dt) 
    #return it
    return temp['time'], temp['nttl'], temp['estr']


def match_to_video(hds,hd_timestamps,ahvs,ahv_timestamps,speeds,speed_timestamps,spike_timestamps,video_timestamps):
    ''' match HDs, AHVs, and wheel speeds to the closest video frame '''
    
    downsampled_hds = np.zeros_like(video_timestamps)
    downsampled_ahvs = np.zeros_like(video_timestamps)
    downsampled_speeds = np.zeros_like(video_timestamps)
    spike_train = np.zeros_like(video_timestamps)
    
    for i in range(len(video_timestamps)):
        
        downsampled_hds[i] = hds[bisect.bisect_left(hd_timestamps,video_timestamps[i],hi=len(video_timestamps)-1)]
        downsampled_ahvs[i] = ahvs[bisect.bisect_left(ahv_timestamps,video_timestamps[i],hi=len(video_timestamps)-1)]
        downsampled_speeds[i] = speeds[bisect.bisect_left(speed_timestamps,video_timestamps[i],hi=len(video_timestamps)-1)]
        
    for i in range(len(spike_timestamps)):
        
        spike_train[bisect.bisect_left(video_timestamps,spike_timestamps[i],hi=len(video_timestamps)-1)] += 1
        
    return downsampled_hds, downsampled_ahvs, downsampled_speeds, spike_train


def process_data(t_filename,smi_filename,orientation_filename,saccade_filename,speed_filename,event_filename):
    ''' process the data to get HD, AHV, wheel speed, pupil pos, and spike train matched to video timestamps '''
    
    event_ts, event_ttl, event_str = grab_nev_data(event_filename)
    recording_start_ts = event_ts[0]

    spike_timestamps = (100*load_spike_timestamps(t_filename) - recording_start_ts) / 1000000.
    video_timestamps = (load_video_timestamps(smi_filename) - recording_start_ts) / 1000000.
    
    sac_data = scipy.io.loadmat(saccade_filename)
    
    raw_ahvs = sac_data['AHV_tsd']['data'][0,0].flatten()
    ahv_timestamps = sac_data['AHV_tsd']['tvec'][0,0].flatten()
    pupil_pos = sac_data['tsdHdeg']['data'][0,0].flatten()
    
    raw_hds, hd_timestamps = get_hds(orientation_filename)
    
    with h5py.File(speed_filename, 'r') as file:
        raw_speeds = np.array(file['wheel_speed_tsd']['data']).flatten()
        speed_timestamps = np.array(file['wheel_speed_tsd']['tvec']).flatten()
        
    downdampled_hds, downsampled_ahvs, downsampled_speeds, spike_train = match_to_video(raw_hds,hd_timestamps,raw_ahvs,ahv_timestamps,raw_speeds,speed_timestamps,spike_timestamps,video_timestamps)

    return downdampled_hds, downsampled_ahvs, downsampled_speeds, pupil_pos, spike_train
