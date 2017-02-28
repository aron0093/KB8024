/input

-parsed data stored for later use

/ouput

-predictions and test outputs

/scripts

-python scripts for data parsing, input generation, predictor building and training.

	-dense_data_parser: Functions to extract data from raw file, return numpy arrays to be used as input for sklearn svm. Sequence 		 vectors and labels are returned in separate arrays. Vectorsisation can be either pssm based or single sequence.
	-cv_set_gen: Function to generate cross validation sets. Currently not used as sklearn.cross_val_score is prefered.
	-cross_val_score: Scoring various SVM categories.
	-run_all: Final script to process data and build the model.
	-test: staging area for testing WIP code.
