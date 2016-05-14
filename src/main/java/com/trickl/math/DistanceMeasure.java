package com.trickl.math;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public interface DistanceMeasure {

    double getDistance(DoubleMatrix1D X, DoubleMatrix1D Y);

    DoubleMatrix2D getDistances(DoubleMatrix2D X);    
}
