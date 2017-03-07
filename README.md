# Project in Molecular Life Science
Part of the Master's course in Molecular Techniques in Life Science at Stockholm University

##OBJECTIVE: 
 Develop a SVM predictor for that predicts, signal peptide, transmembrane and globular regions.
 
##METHODOLOGY:
*	Data parsing and feature extraction. Include multiple window size functionality.
*	Creating SVM architecture using scikit-learn and choosing initial parameters.
*	Train and score SVM using cross validation for single sequence data per sequence.
*	Optimize performance by modifying parameters.
*	Extract evolutionary relationships using PSI-BLAST and modify features.
*	Train and score SVM using cross validation for multiple sequence per sequence.
*	Optimize performance by modifying parameters.
*	Extract the data from 50 other proteins and test the predictor.
*	Analyse the results, compare previous results and also with the state of the art.
*	Analyse performance compared to other ML methods.

##FUNCTIONALITY:
*   Parse raw data from standard format into pandas dataframe supporting custom window_sizes.
*   Generate PSSM using psiblast locally on multiple CPUs over ssh.
*   Generate support vectors using PSSM or default encoding.
*   Compare different kernels, window_sizes and parametes to suit your data. Customised for large, unbalanced datasets for multiclass classification.
*   Final model is made by bagging 10 LinearSVC estimators though the script is easily modifiable to use any classifier.
*   Model is stored and can be accessed via bin. Scores and comparisons with other classifiers are available.

##OPERATIONALITY:
*   Run SignalP_predict.py file to choose between models created in the SignalP project.

## DEVELOPMENTAL STEPS:
*   Data parsing for defaulte encoding
*   Generate PSSMS.
*   Data parsering for MSA encoding.
*   Comparing performance of various kernels and window sizes using default encoding.
*   Optimising C via grid search using LinearSVC and optimum window size
*   Comparing default encoding with MSA encoding.
*   Further optimising window size.
*   Using BaggingClassifier to boost performance
*   Comparing with RandomForestClassifier
*   Testing on 50 randomn proteins.
