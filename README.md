# HD-system-brainstem
Code for analyzing data from recordings in mouse brainstem head direction circuit. 

Makes extensive use of the vandermeerlab codebase: https://github.com/vandermeerlab/vandermeerlab

Experiments were carried out with headfixed mice on a running wheel. The platform situated below the mice could be rotated like a lazy susan so that the mice experienced periods of sinusoidal rotation at varying speeds. The platform speed is referred to as angular head velocity (AHV). 

A camera and infrared light source were situated near the left eye to record eye movement data. 

16 channel silicon probes from Neuronexus were advanced through a craniotomy to reach brainstem targets in the earliest part of the HD circuit; namely, nucleus prepositus (NPH), Gigantocellular nucleus (Gi) -- just below NPH, and supragenual nucleus (SGN). 

In some sessions, transgenic mice with viral injections into the dorsal tegmental nucleus (DTN) recieved pulses of blue laser light with combined optical silicon probes. The purpose of the injection + laser was to identify photoresponsive cells projecting to DTN (i.e., optotagging). Because the shutter that operates the laser makes an audible sound, a "dummy" sound stimulus was added in later sessions to control for neurons that might be response to sound or a surprising stimulus. 

Below is a list of important files/file types for each session. 
-----------------------------------------------------------------------------------------------------------------------------------------------------------
All files start with the mouse SSN. The SSN is the mouse identity and the date of the recording. M___-YEAR-MONTH-DAY-recording number(if more than one for the same day)		
	
**.mp4**		Raw video file from the eye-tracking camera. The eye camera uses infrared light.

**.ntt**		Raw neuralynx tetrode data file. Contains the event times and waveforms for all events that cross threshold. Each session has 8 files (1 per tetrode), but only a few will actually contain data with good cells.

**.t** 		Processed file that corresponds to one unique neuron that has been identified manually after the cluster cutting proccess. Contains time stamps and waveforms for the identified neuron only. 

**.smi** 		This file contains information about the time of each camera frame. It is used to align data in the mp4 file.

**proc.mat**		This file contains the eye tracking (plus some extras) information that is extracted from the mp4 the Facemap piece of software. Some manual input from the experimenter is required.  https://github.com/MouseLand/facemap

**-saccades-edited.mat**		After the mostly automated processing with Facemap, the x,y pupil position data is then thresholded to identify saccade times and amplitudes. There is a manual curation step for removing false positives and adding false negatives. The amount of manual tuning depends on the quality of the eyetracking. 
		variables are ...nasalSaccades (timestamps)		nasalAmplitudes (amplitude, in pixels)			temporalSaccades (timestamps)		                      temporalAmplitudes (amplitude, in pixels)		
		
**CSC33.Ncs**		CSC = continuously sampled channel. The neuralynx file type contains data from the platform encoder, which tracks the angle of the recording platform. Since the mouse is head-fixed, this corresponds to HD. This data is used to calculate AHV. 

**CSC34.Ncs**		Channel 1: Position data from the wheel that the mouse runs (or sits) on during the session. Quadrature encoding is used to encode the wheel position. 	

**CSC35.Ncs**		Channel 2 of wheel encoder

**-AHV_StationaryTimes.mat**		These are manually identified start and end times for periods when the platform is stationary.

**.keys**	This file contains important metadata about the recording session. Too many entries to list here. It does contain a line for which hemisphere each TT is likely in. It also may contain the position in the brain (A/P, M/L, D/V) and an estimate of confidence in that location (based on histology).

**Event.nev**		Labels and timestamps for events that occur during a session, usually signalled by a TTL pulse to the Lynx input board. The relevent events here are 'Laser On' and 'Laser Off' timestamps. Also, 'Starting Recording' and 'Stopping Recording'

Below is a list of important functions for extracting/plotting information about AHV and eye movements. 
-----------------------------------------------------------------------------------------------------------------------------------------------------------

All of the relevant files for a recording session are in one folder. Once inside that folder, use **sd = LoadSessionData([])** to calculate and gather all of the relevant information into one structure **[sd]** in the matlab workspace. sd will be an input for most functions. 

To calculate an AHV tuning curve, use **getAHV_TC.m**. Set 'doPlot'=1 to plot each TC. The default setting will smooth the data. One parameter that you may want to change is minOCC. This value determines the minimum number of samples needed in a given AHV bin to include that data. minOcc of 100 requires only 0.5 sec of occupancy to include those values. I usually set this to 200 to include at least 1 sec of data. But one may want to set it higher. 

Due to eye movement tuning (and natural variability) the tuning curve variance can be high. We prefer to plot the raw firing rate data (FR x AHV) and overlay the binned (average) tuning curve. Use **plotAHVscatter** to visualize the raw firing rate data and the tuning curve together. 


























