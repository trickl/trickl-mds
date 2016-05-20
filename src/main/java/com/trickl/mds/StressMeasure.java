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
    private double stress = 0;
    private double normalised_denominator = 0;
    private double kruskal_denominator = 0;
    
    public double calculate(DoubleMatrix2D X, DoubleMatrix2D R) {
       DoubleMatrix2D W = R.like(R.rows(), R.columns());
       R.forEachNonZero((int i, int j, double value) -> {
            W.set(i, j, 1.0);
            return value;
       });
       return calculate(X, R, W);
    }

    public double calculate(DoubleMatrix2D X, DoubleMatrix2D R, DoubleMatrix2D W) {

        if (X.rows() != R.rows()) {
            throw new IllegalArgumentException("X must have the same number of points as R.");
        }
        if (W.rows() != R.rows() || W.columns() != R.columns()) {
            throw new IllegalArgumentException("W must have the same dimensions as R.");
        }

        stress = 0;
        normalised_denominator = 0;
        kruskal_denominator = 0;

        R.forEachNonZero((int i, int j, double value) -> {
            if (i > j) {
                // Actual distance
                double d2 = 0;

                for (int k = 0; k < X.columns(); ++k) {
                    d2 += Math.pow(X.get(i, k) - X.get(j, k), 2.);
                }

                stress += Math.pow(Math.sqrt(d2) - value, 2.0) * W.get(i, j);
                normalised_denominator += Math.pow(value, 2.0) * W.get(i, j);
                kruskal_denominator += d2 * W.get(i, j);
            }
            
            return value;
        });
        

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
