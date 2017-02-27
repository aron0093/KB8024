# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 13:11:47 2017

@author: Revant Gupta
"""

def cv_set_gen(X, Y, K, randomise=False):

    import numpy as np

    assert len(X) == len(Y)
    
    if randomise:
        randomize = np.arange(len(X))
        np.random.shuffle(randomize)
        x = X[randomize]
        y = Y[randomize]
        
    else:
        x = X
        y = Y

    X_sets = np.zeros(0)
    Y_sets = np.zeros(0)

    if len(x) % K ==0:
        X_sets = np.array(np.split(x, K, axis=0))
        Y_sets = np.array(np.split(y, K, axis=0))
    else:
        print("Cant divide sets into %d sets, choose a divisible number"%(K))
        return
                          
    for k in range(K):
        train_X = np.concatenate(([sets for i, sets in enumerate(X_sets) if i!=k]),axis=0)
        val_X = X_sets[k]
        train_Y = np.concatenate(([sets for i, sets in enumerate(Y_sets) if i!=k]),axis=0)
        val_Y = Y_sets[k]

	yield train_X, val_X, train_Y, val_Y


