/bash

-bash scripts to build and train predictor

/input

-parsed data stored for later use

/ouput

-predictions and test outputs

/scripts

-python scripts for data parsing, input generation, predictor building and training.

	-SS_parser: Functions to extract data from raw file, return numpy arrays to be used as input for sklearn svm. Sequence vecctors and labels are returned in separate arrays. Vectorsisation is based on simple dictionary with no MSA.
	-MS_parser: Same as SS_parser but vectorisation is based on PSSM derived from PSIBLAST.
	-cv_set_gen: Function to generate cross validation sets. Currently not used as sklearn.cross_val_score is prefered.
	-Other scripts: calling and executing functions. Check bash/runall.sh for complete description of steps.