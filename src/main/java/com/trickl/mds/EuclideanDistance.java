package com.trickl.mds;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

public class EuclideanDistance {
    
    public DoubleMatrix2D getSimilarities(DoubleMatrix2D X) {
        DoubleMatrix2D R = DoubleFactory2D.dense.make(X.rows(), X.rows());
        
        for (int i = 0; i < X.rows(); ++i) {
            for (int j = 0; j < X.rows(); ++j) {
                if (i == j) {
                    R.set(i, i, 0);
                }
                else {
                    double d2 = 0;
                    for (int k = 0; k < X.columns(); ++k) {
                        d2 += Math.pow(X.get(i, k) - X.get(j, k), 2);
                    }

                    R.set(i, j, Math.sqrt(d2));
                }
            }        
        }
        return R;
    }    
}
