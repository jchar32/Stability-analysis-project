Project: The number of steps for stable real-world, unsupervised walking data using a shoe-worn inertial sensor

## Overview
The dataset and code relate to a project from my dissertation involving biomechanical gait data rising from a foot-worn inertial sensor.
Participants completed a week of unsupervised community walking. Biomechanical outcomes were extracted from the data and a stability
analysis was performed to identify how many steps were required before the outcome stabilized. 

## Components
Data
- files are structured as ROW: participant observation , COL: gait or demographic outcome

1) gait_metric_means.csv
mean values for each gait outcome over all the collected steps from each participant

2) percent_similar_data.csv
The percentage of similar data when stability was reached for each outcome and participant.

3) steps2stability_data.csv
- The number of steps when stability was reached for each outcome and participant.

Code
1) mainAnalysis.R
script ingests data files, cleans data, and performs the primary analysis of the project (generalized linear mixed effects model).

## Outcome Dictionary
fpa_classif - foot progression angle during midstance
peakacc_norm_HS	- Euclidean norm of linear acceleration at heel strike
peakacc_norm_TO	- Euclidean norm of alinear acceleration at toe-off
peakgyr_norm_HS	- Euclidean norm of angular velocity at heel strike
peakgyr_norm_TO	- Euclidean norm of angular velocity at toe-off
stance_time - time from heel strike to toe off
foot_strike_angle - angle of foot at heel strike in sagittal plane. 
group - participant grouping indicator: HA=healthy adult, KOA = knee ostoearthritis
ID - numeric id value. int: 0-30
Age - age: units years. float: positve
Sex - sex. categorical: M=male, F=female
BMI - body mass index: units kg/m^2. float: positive 
Klgrade - severity of structure signs of osteoarthritis. int: Empty, 2, 3, 4

## Usage
1) Open mainAnalysis.R
2) Run script section by section. NOTE: ensure the first section installs and loads all necessary libraries
3) Data tables are output into /data/ folder