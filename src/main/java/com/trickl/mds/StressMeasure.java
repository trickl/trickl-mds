package com.trickl.mds;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

public class StressMeasure {
    
    public enum Type {

        Raw,
        Normalized,
        Kruskal
    }
    
    public StressMeasure(Type type) {
        this.type = type;
    }
    
    private final Type type;
    
    public double calculate(DoubleMatrix2D X, DoubleMatrix2D R) {
       DoubleMatrix2D W = DoubleFactory2D.dense.make(R.rows(), R.columns());
       W.assign(1.0);
       return calculate(X, R, W);
    }

    public double calculate(DoubleMatrix2D X, DoubleMatrix2D R, DoubleMatrix2D W) {

        if (X.rows() != R.rows()) {
            throw new IllegalArgumentException("X must have the same number of points as R.");
        }
        if (W.rows() != R.rows() || W.columns() != R.columns()) {
            throw new IllegalArgumentException("W must have the same dimensions as R.");
        }

        double stress = 0;
        double normalised_denominator = 0;
        double kruskal_denominator = 0;

        for (int i = 0; i < R.rows(); ++i) {
            for (int j = 0; j < R.rows(); ++j) {
                if (i > j) {
                    // Actual distance
                    double d2 = 0;

                    for (int k = 0; k < X.columns(); ++k) {
                        d2 += Math.pow(X.get(i, k) - X.get(j, k), 2.);
                    }

                    stress += Math.pow(Math.sqrt(d2) - R.get(i, j), 2.0) * W.get(i, j);
                    normalised_denominator += Math.pow(R.get(i, j), 2.0) * W.get(i, j);
                    kruskal_denominator += d2 * W.get(i, j);
                }
            }
        }

        if (type == Type.Normalized) {
            stress /= normalised_denominator;
        }
        if (type == Type.Kruskal) {
            stress /= kruskal_denominator;
        }
        
        stress = Math.sqrt(stress);        

        return stress;
    }
}
