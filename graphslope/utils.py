from __future__ import division

import numpy as np

def FDR(D, xhat, x0, tol=1e-6):
    dh = D.T.dot(xhat)
    d0 = D.T.dot(x0)
    supph = np.abs(dh) > tol
    supp0_c = np.abs(d0) <= tol
    false_and_pred = supph * supp0_c
    FDetect = np.count_nonzero(false_and_pred)
    Detect = np.count_nonzero(supph)
    if Detect < 1.:
        FDR = 0.
    else:
        FDR = FDetect / Detect
    return FDR


def TDR(D, xhat, x0, tol=1e-6):
    dh = D.T.dot(xhat)
    d0 = D.T.dot(x0)
    supph = np.abs(dh) > tol
    supp0 = np.abs(d0) > tol
    true_and_pred = supph * supp0
    TDetect = np.count_nonzero(true_and_pred)
    Detect = np.count_nonzero(supp0)
    if Detect < 1.:
        TDR = 0.
    else:
        TDR = TDetect / Detect
    return TDR

