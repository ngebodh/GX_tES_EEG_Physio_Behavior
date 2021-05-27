# Dataset of Concurrent EEG, ECG, and Behavior with Multiple Doses of transcranial Electrical Stimulation

This repository contains supporting code for the [GX dataset](https://zenodo.org/record/4456079#.YK8ak6hKiF5).
For an in-depth description of the dataset please see the accompanying dataset publication [TO DO : ADD NATURE DATA INFO]


## Description:

A dataset combining high-density electroencephalography (EEG) with physiological and continuous behavioral metrics during transcranial electrical stimulation (tES; including tDCS and tACS). Data includes within subject application of nine High-Definition tES (HD-tES) types targeted three brain regions (frontal, motor, parietal) with three waveforms (DC, 5Hz, 30Hz), with more than 783 total stimulation trials over 62 sessions with EEG, physiological (ECG or EKG, EOG), and continuous behavioral vigilance/alertness metrics.


### Folder set up
Generally the project folder is set up as below:

```
+GX
|--Analysis
    +--GX_Exp1_CTT_GeneralAnalysis.m
    +--GX_Exp2_CTT_GeneralAnalysis.m
    ...
|--Data
    +--0101
        +--0101
            +--ptracker-0101.csv
            +-- ptracker-summary-0101..txt
        +--GX_01_2019-09-24_15-45-53.cnt
        +--GX_01_2019-09-24_15-45-53.evt
        +--MATLABfilestream0101924.mat
        +--MATLABfilestream0101924.txt
    +--0102
    +--0103
    +--0104
    ... 
|--Documents
|--Results

```

### What do I need to get started?
You can get started right away by using the downsampled `.mat` files linked to in the main data repository [see 'Extras' here](https://zenodo.org/record/3837213#.YK9ThahKjZQ). 
The `.mat` files are compatible with MATLAB and Python (other platforms have not been tested). Each of the files contains the ~70 min recording combined with the CTT data and information on the stimulation trials. 

##### What's in the dowsampled ``.mat`` files?
The ``.mat`` files contain a matlab structre that contains the following:
* DSamp
    * triggers <- These are all the labeled EEG/Stimulation start/stop triggers
    * EEGdata <- Contains the downsampled EEG/ECG/EOG voltage data dims: 35 channelss X ~4E6 samples
    * fs <- The downsampled sampling frequency of the data : 1000 Hz
    * fsOld <- The original sampling frequency of the data
    * time <- Time vector for the data. Should be 1 X ~4E6
    * label <- Contains the channel label information. BIP1= ECG, BIP2=EOG, RESP1= N/A
    * nchan <- The number of channels in the data
    * rate <- Redundant to fs, sampling rate of data
    * npt <- Number of data points ~4E6
    * Subj <- Subject and session that data belong to. I.e. 0302 - Subject 03 session 03
    * ptrackerPerf <- The CTT data deviation/ the behavioral data
    * ptrackerTime <- Time vector for the CTT data
    * ptrackerfs <- The sampling frequency for the CTT data 100 Hz.

NOTE: The CTT and EEG data can be time alligned by using the 1st trigger in the EEG. The CTT starts at the same time the 1st trigger want sent to the EEG. 


TO DO: Insert pictures

### Going through the code


#### Experiment 1- CTT
Once all the raw data is downloaded and all supporting toolboxes are obtained, the ```GX_Exp1_CTT_GeneralAnalysis.m``` script can be run to analyze and extract data and outcomes for Experiment 1. This script mainly looks at the outcomes of the CTT data. Please note that directories were hard coded and would have to be recoded for the appropriate directories.

TO DO: Insert pictures


#### Experiment 2- CTT
To explore the outcomes for Experiment 2 the ```GX_Exp2_CTT_GeneralAnalysis.m``` script can be used. This script mainly looks at the outcomes of the CTT data. 

TO DO: Insert pictures

#### I want to look at the stimulation trials
To pull out and look at all the stimulation trials for the whole study you can do so by running the ``GX_PullingDataIntoTrials_PlottingTopoplots.m`` script. This script runs through all the data and pull out the EEG and CTT data into 30 sec Pre During Post stimulation trials. Please pay attention to the code flags which allow for the plotting and saving for different things. 

TO DO: Insert pictures


#### The files are too large how do I downsample them?!
Since the raw data are sampled at 2k Hz moving and loading files may become difficult on some machines. If you would like to downsample the data please use the ```GX_DataDownSample.m``` script. The script features a GUI that allows you to paste in file names, locations, and downsample factor according to how much you want to downsample the data. 

TO DO: Insert pictures

#### I want to explore the demographic data of the study
To pull in and look at the demographic data of the study you can use the ``GX_PlottingDemographicInfo.m`` script. This pulls in the tabulated demogrpahic data that were compiled over the course of the study and plots them. 

TO DO: Insert pictures

## Acknowledgments

Portions of this study were funded by X (formerly Google X), the Moonshot Factory. The funding source had no influence on study conduction or result evaluation. MB is further supported by grants from the National Institutes of Health: R01NS101362, R01NS095123, R01NS112996, R01MH111896, R01MH109289, and (to NG) NIH-G-RISE T32GM136499.
