## 2024_Murphy_Belief-updating-psychosis-proneness
Analysis and task code for

**Murphy PR, Krkovic K, Monov G, Kudlek N, Lincoln, T & Donner TH (2024).** [**Individual differences in belief updating and phasic arousal are related to psychosis proneness.**](https://www.biorxiv.org/content/10.1101/2024.01.14.575567v1) ***bioRxiv*.**

Raw behavioural and eye-tracking data to accompany this code will be made available upon publication.

Code shared here was developed and tested using Matlab R2019b and the [FieldTrip Toolbox](https://www.fieldtriptoolbox.org/) version 20160221.

Further detail on analysis and task scripts is provided below. For questions, contact the first author at murphyp7@tcd.ie.

### behaviour:
Scripts for analysing behavioural data from change-point decision-making task and delayed match-to-sample working memory task. Main script for analysis of decision-making is `analyse_dm.m`, which computes choice accuracies, psychophysical kernels and other key measures, assesses their relationships to scores from the CAPE questionnaire, and creates a variety of figures including many that feature as figure panels in the manuscript. `analyse_wm.m` computes average response accuracy on the working memory task and assesses its relationship to scores from the CAPE questionnaire.
### modelling:
Scripts for fitting a variety of model variants to choice data from the decision-making task. Code for each distinct model variant is contained in individual subdirectories.
-	Models are fit via particle swarm optimization ([Birge, 2003, *IEEE Swarm Intelligence Symposium*](https://ieeexplore.ieee.org/document/1202265); [code here](https://www.mathworks.com/matlabcentral/fileexchange/7506-particle-swarm-optimization-toolbox)), with very slightly adapted version of code included here in the `particle_swarm_Glaze` subdirectory (changes are commented with `%%% PM` and made code work with global variables passed from higher-level scripts).
-	`model_comp.m`: compares goodness of fit of different model variants.
### param_recovery:
Scripts for running parameter recovery analyses for different model variants, referred to as Model 1 (only H and noise as free parameters; `Glaze_parameter_recovery.m`) and Model 3 (H, noise and IU as free parameters; `Glaze_parameter_recoveryIU.m`) in Supplementary Figure 2 of the manuscript. Results from the analysis are plotted via `plot_param_recovery.m`. Output of simulations run via `Glaze_parameter_recovery.m` and `Glaze_parameter_recoveryIU.m` is included in the `sims` subdirectory.
### pupil:
Scripts for preprocessing and analysing pupillometric and eye-tracking data.
- First step is to convert native .idf files from SMI eye-tracker to .txt using the manufacturer's 'IDF Converter' software utility. 
-	`convertSMI.m`: reads eye-tracker data from .txt files into Matlab.
-	`loop_procASC.m`: cleans data read in via previous function and pulls trial information.
-	`loop_interpSS.m`: interpolates pupil and gaze position time-series.
-	`analyse_pupil.m`: analyses pupil data for a single subject.
-	`plot_output.m`: creates a variety of figures for visualising results of pupil analyses, including many that feature as figure panels in the manuscript.
