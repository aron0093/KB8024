/input

-parsed data stored for later use

    -pssms: Contains all pssm from psiblast.
    -Window_3: Contains all proteins vectorised with a 3 residue frame. 

/ouput

-predictions and test outputs

    -Window_kernel: Contains data and graphs for comparison.
    -model: Contains trained models
    -Linear_grid_search: Contains output results for grid search optimisation of C parameter for LinearSVC

/scripts

-python scripts for data parsing, input generation, predictor building and training.

	-dense_data_parser: Functions to extract data from raw file, return numpy arrays to be used as input for sklearn svm. Sequence vectors and labels are returned in separate arrays. Vectorsisation is hardcoded and doesn't use MSA. Added a function to split data at protein level for cross validation.
	-pssm_data_parser: Generate vectors as dense_data_parser but based on pssms.
	-cv_set_gen: Function to generate cross validation sets. Currently not used as sklearn.cross_val_score is prefered. Shifted to Obsolete
	-pssm_func: Contains functions to generate pssm using psiblast on multiple computers over ssh.
	-pssm_ssh: Script to run psiblast by dividing data and executing on all computers.
	-window_kernel_CV: Scoring various SVM kernels and window sizes.
	-Linear_grid_search: Grid search to optimise C parameter using LinearSVC.
	-model_builder: Final script to process data and build the model.
	-test: staging area for testing WIP code.
