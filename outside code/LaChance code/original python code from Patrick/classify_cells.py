# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:47:16 2024

@author: plachanc
"""

import os

import data_processing
import glm


if __name__ == '__main__':

    example_data_dir = os.getcwd() + '/example_data'

    t_filename = example_data_dir + '/M247-2022-01-14-2/M247-2022-01-14-2-TT08_1.t'
    smi_filename = example_data_dir + '/M247-2022-01-14-2/M247-2022-01-14-2-VT1.smi'
    orientation_filename = example_data_dir + '/M247-2022-01-14-2/M247-2022-01-14-2_orientation_tsd.mat'
    saccade_filename = example_data_dir + '/M247-2022-01-14-2/M247-2022-01-14-2-saccades-edited.mat'
    speed_filename = example_data_dir + '/M247-2022-01-14-2/M247-2022-01-14-2_wheel_speed_tsd.mat'
    event_filename = example_data_dir + '/M247-2022-01-14-2/M247-2022-01-14-2-Events.Nev'

    hds, ahvs, speeds, pupil_pos, spike_train = data_processing.process_data(t_filename,smi_filename,orientation_filename,saccade_filename,speed_filename,event_filename)
    
    best_model = glm.run_classifier(hds, ahvs, speeds, pupil_pos, spike_train)

    print(best_model)