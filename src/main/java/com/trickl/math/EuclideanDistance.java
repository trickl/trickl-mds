package com.trickl.math;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class EuclideanDistance implements DistanceMeasure {
    
    @Override
    public DoubleMatrix2D getDistances(DoubleMatrix2D X) {
        DoubleMatrix2D R = DoubleFactory2D.dense.make(X.rows(), X.rows());
        
        for (int i = 0; i < X.rows(); ++i) {
            for (int j = i; j < X.rows(); ++j) {
                if (i == j) {
                    R.set(i, i, 0);
                }
                else {
                    double distance = getDistance(X.viewRow(i), X.viewRow(j));
                    R.set(i, j, distance);
                    R.set(j, i, distance);                    
                }
            }        
        }
        return R;
    }    
    
    @Override
    public double getDistance(DoubleMatrix1D X, DoubleMatrix1D Y) {        
        double d2 = 0;
        for (int k = 0; k < X.size(); ++k) {
            d2 += Math.pow(X.get(k) - Y.get(k), 2);
        }
        return  Math.sqrt(d2);
    }  
}
