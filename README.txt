In this fMRI study, we investigated how beliefs about sensory uncertainty influence early visual information processing. We used unreliable feedback interventions (in which invalid feedback was given on 50% of trials) to induce such beliefs and used reliable feedback as a control condition. The experiment used a between-group design, and fMRI data from 29 healthy human participants (n= 15 in the unreliable feedback group, 14 in the control/ reliable feedback group) were obtained. Participants performed a challenging visual orientation discrimination task by identifying visual gratings presented in circular annuli (inner r = 1.32°, outer r = 6.69°) as clockwise or counter-clockwise with reference to one of two implicit diagonal references (45° or 135°). The stimuli deviated from the diagonals by an angle defined prior to the main part of the experiment using a staircase procedure. The main part of the experiment had three phases (1) a baseline pre-intervention phase (no fb), (2) an intervention phase with reliable/unreliable fb and (3) a post-intervention phase (no fb). The main analyses of interest is the comparison of pre/post intervention changes in behaviour and pattern distinctness between the reliable and unreliable feedback groups. 

Here the following data and codes are given:

BEHAVIOURAL
behavioural_performance.m : code showing the estimation of task performance from the raw data (data_across_blocks_and_subjects.mat). Run-wise estimates of task performance provided in tbl_DV1.mat.

Data organisation in data_across_blocks_and_subjects.mat :
Aalldata is a cell comprsing raw data. Rows: subject indices (6-35. Reliable group: 7,9,..., unreliable group: 6,8,...) and columns: run number (1-24).

Each cell consists of 2D matrices with rows as the number of trials (with responses) and columns as the following: (run number (1-8), trial number (1-32),stimulus type (1-4), response type after correcting for the random response mapping (1-2), response mapping (1-2), trial onset time, stimulus onset time, response onset time.

FMRI
fmri_cvMANOVA_code.m : shows the invocation of the cvMANOVA algorithm for the ROI and searchlight analyses

fmri_extract_anat_visual_voxels.m : creation of the V1 mask from the probabilistic map.

fmri_univariate_analysis_V1.m : code showing univariate analysis of V1 voxels.

cvMANOVA_V1results.mat: pattern distinctness estimates for each test phase and orientation pair (col1: pair1/pre, col2: pair2/pre, col3: pair1/post, col4: pair2/post, col5: fbtype = reliable, 2=unreliable feedback, col6 = threshold angular deviation of stimuli from the 45° or 135° references).
results_univ_act.csv: overall activity within the V1 voxels used for multivariate analyses

Eyetracking
data_eyetracking.mat - eyetracking data after preprocessing and baseline correction. Data is only from the stimulus presentation windows on each trial (approx 10 data points per trial, spaced 20ms apart. Stimulus display time 200ms)
data_across_blocks_and_subjects_incl_no_resp.mat - behavioural data across participants (n= 30) and runs (pre,intervention and post = 24 runs) including the missing trials used for the multivariate analyses of the eyetracking data
Decoding_accuracy_backup.csv - decoding accuracy at each test phase (pre/post) and for each participant. fbtype 1: reliable, 2: unreliable. Fourth column has subject numbers comparable to the behavioural and fMRI data.
eye_DataAnalysis.m - eyetracking data processing to obtain fixation percentage and to decode stimuli from eyetracking data similar to the method reported by Thielen et al., 2019, Scientific Reports

Debriefing data
Debriefing_ind.xlsx - has the percentages and ratings (converted to integers, scale mentioned in the file) used to estimate reliability on fb, awareness of feedback manipulation and task motivation.
