# HD-system-brainstem
Code for analyzing data from recordings in mouse brainstem head direction circuit. 

Makes extensive use of the vandermeerlab codebase. Use [this commit](https://github.com/vandermeerlab/vandermeerlab/tree/e4ff8327ee8b7b65a856499f25bad2c7d57524dc).

Execute these lines when you begin analysis to set up the correct MATLAB path for analysis. 

```
%% set up
restoredefaultpath;
clear classes;
addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared')); % replace with your path
addpath(genpath('D:\My_Documents\GitHub\HD-system-brainstem')); % replace with your path
```

Several example sessions are available for download in [this google drive folder](https://drive.google.com/drive/folders/11RaD-QtlHRowaEffT6OR2wNw98RGl02i?usp=sharing). One session has prototypical AHV cells. One session has prototypical eye movement modulated cells. One sessions has an opto-tagged cell. 

![exp_setup](https://github.com/vandermeerlab/HD-system-brainstem/assets/1922878/25ac3863-7612-477f-a943-fb0cd27e1097)


Experiments were carried out with head-fixed mice on a running wheel. The platform situated below the mice could be rotated like a lazy susan so that the mice experienced periods of sinusoidal rotation at varying speeds. The platform speed is referred to as angular head velocity (AHV). 

A camera and infrared light source were situated near the left eye to record eye movement data. The video data were analyzed with [Facemap](https://github.com/MouseLand/facemap) to extract a time series of pupil position. Velocity was calculated (from position) and thresholded to get candidate saccade event times. These candidate events were screened by the experimenter to exclude events that appeared to be false positives and to include missed events. Pupil position coordinates are in units of camera pixels (x,y). The technically difficult process of calibrating and measurement needed to obtain degrees of visual angle was not conducted. 

![sample eye trace](https://user-images.githubusercontent.com/16581827/235464246-c8276ff3-332f-431e-b532-2ef623d9a3a3.JPG)


16 channel silicon probes from Neuronexus were advanced through a craniotomy to reach brainstem targets in the earliest part of the HD circuit; namely, nucleus prepositus (NPH), Gigantocellular nucleus (Gi) -- just below NPH, and supragenual nucleus (SGN). 

In some sessions, transgenic mice with viral injections into the dorsal tegmental nucleus (DTN) recieved pulses of blue laser light with combined optical silicon probes. The purpose of the injection + laser was to identify photoresponsive cells projecting to DTN (i.e., optotagging). Because the shutter that operates the laser makes an audible sound, a "dummy" sound stimulus was added in later sessions to control for neurons that might be responsive to sound or a surprising stimulus. 

Below is a list of important files/file types for each session. 
-----------------------------------------------------------------------------------------------------------------------------------------------------------
All files start with the mouse SSN. The SSN is the mouse identity and the date of the recording. M___-YEAR-MONTH-DAY-recording number(if more than one for the same day)		
	
**.mp4** &nbsp;-  Raw video file from the eye-tracking camera. The eye camera uses infrared light.

**.ntt**		-  Raw neuralynx tetrode data file. Contains the event times and waveforms for all events that cross threshold. Each session has 8 files (1 per tetrode), but only a few will actually contain data with good cells.

**.t** 		-  Processed file that corresponds to one unique neuron that has been identified manually after the cluster cutting proccess. Contains time stamps and waveforms for the identified neuron only. 

**.smi** 		-  This file contains information about the time of each camera frame. It is used to align data in the mp4 file.

**proc.mat**		-  This file contains the eye tracking (plus some extras) information that is extracted from the mp4 by the [Facemap](https://github.com/MouseLand/facemap) software. Some manual input from the experimenter is required. 

**-saccades-edited.mat**		-  After the mostly automated processing with Facemap, the x,y pupil position data is then thresholded to identify saccade times and amplitudes. There is a manual curation step for removing false positives and adding false negatives. The amount of manual tuning depends on the quality of the eyetracking. 
		variables are ...nasalSaccades (timestamps)		nasalAmplitudes (amplitude, in pixels)			temporalSaccades (timestamps)		                      temporalAmplitudes (amplitude, in pixels)		
		
**CSC33.Ncs**		-  CSC = continuously sampled channel. The neuralynx file type contains data from the platform encoder, which tracks the angle of the recording platform. Since the mouse is head-fixed, this corresponds to HD. This data is used to calculate AHV. 

**CSC34.Ncs**		-  Channel 1: Position data from the wheel that the mouse runs (or sits) on during the session. Quadrature encoding is used to encode the wheel position. 	

**CSC35.Ncs**		-  Channel 2 of wheel encoder

**-AHV_StationaryTimes.mat**		-  These are manually identified start and end times for periods when the platform is stationary.

**.keys**	-  This file contains important metadata about the recording session. Too many entries to list here. It does contain a line for which hemisphere each TT is likely in. It also may contain the position in the brain (A/P, M/L, D/V) and an estimate of confidence in that location (based on histology).

**Event.nev**		-  Labels and timestamps for events that occur during a session, usually signalled by a TTL pulse to the Lynx input board. The relevent events here are 'Laser On' and 'Laser Off' timestamps. Also, 'Starting Recording' and 'Stopping Recording'

Below is a list of important functions for extracting/plotting information about AHV and eye movements. 
-----------------------------------------------------------------------------------------------------------------------------------------------------------

All of the relevant files for a recording session are in one folder. Once inside that folder, use **sd = LoadSessionData([])** to calculate and gather all of the relevant information into one structure **[sd]** in the matlab workspace. sd will be an input for most functions. 

To show AHV tuning curves for a single session, do **getAHV_TC([],sd)**. See the function help for options other than the defaults. 

Due to eye movement tuning (and natural variability) the tuning curve variance can be high. We prefer to plot the raw firing rate data (FR x AHV) and overlay the binned (average) tuning curve. Use **plotAHVscatter([],sd)** to visualize the raw firing rate data and the tuning curve together. See the function help for options other than the defaults.

![AHV scatter example](https://user-images.githubusercontent.com/16581827/235242261-32805e02-7141-437a-86eb-c0daeedbe0b4.jpg)

Use **saccadePETH([],sd)** to collect the spike times (relative to zero) and Firing Rate PETHs from a user-chosen window around saccade event times. See the function help for options other than the defaults ([]). Saccade times are calculated from a semi-automated process using the estimated pupil position from Facemap and from a velocity threshold and manual editing. saccadePETH will return the FR for each bin and can plot the raster plot for that time window with the FR overlaid. PETH information is generated separately for nasal and temporal saccades. A +- 200 ms window is good for seeing phasic changes. A +- 2 s window will show the saccade response on top of the underlying AHV oscillation (if present). Use **FR_PETH([],S,t)** if you want to generate a PETH with any general set of times t that you input. 
![saccadePETH example](https://user-images.githubusercontent.com/16581827/235328676-9724619a-ddb0-46bb-a84c-da6385d9097f.jpg)




Use **getWheel_TC([], sd)** to calculate and plot the linear velocity (wheel speed) tuning curve. See the function help for options other than the defaults ([]).
![wheel speed TC](https://user-images.githubusercontent.com/16581827/235240827-e979ac97-9f75-4538-8ef7-955baac154fb.JPG)



Use **getPupil_TC([], sd)** to calculate and plot the eye position tuning curve. See the function help for options other than the defaults ([]).
Efforts were made during surgery to cement the headbar in place so that it was level with the horizontal. There are sessions in which the eye is not perfectly on the horizontal. DLC could be used to correct.

![pupil position TC example](https://user-images.githubusercontent.com/16581827/235252994-2643811d-d877-4a60-8244-62af6bd075c4.JPG)



Use [blank] function to plot all of the relevant information in one figure with subplots. 




