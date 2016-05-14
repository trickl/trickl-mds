package com.trickl.math;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import java.util.function.UnaryOperator;

public class CenteringMatrix {

    private static DoubleMatrix1D getSum(DoubleMatrix2D mat, UnaryOperator<Double> func, boolean rows) {
        final DoubleMatrix1D sum = mat.like1D(rows ? mat.columns() : mat.rows());
        mat.forEachNonZero((int first, int second, double value) -> {
            int index = rows ? second : first;
            sum.setQuick(index, sum.getQuick(index) + func.apply(value));
            return value;
        });
        return sum;
    }

    private static int[] getNonZeroCount(DoubleMatrix2D mat, boolean rows) {
        final int[] nonZeroCount = new int[rows ? mat.columns() : mat.rows()];
        mat.forEachNonZero((int first, int second, double value) -> {
            int index = rows ? second : first;
            nonZeroCount[index] += 1;
            return value;
        });
        return nonZeroCount;
    }
    
    // Double centring, so origin is centroid of all points
    // https://en.wikipedia.org/wiki/Centering_matrix
    // Also See http://forrest.psych.unc.edu/teaching/p230/Torgerson.pdf
    public static DoubleMatrix2D getDoubleCentered(DoubleMatrix2D input) {
        return getDoubleCentered(input, value -> value);
    }

    // Double centring, so origin is centroid of all points
    // https://en.wikipedia.org/wiki/Centering_matrix
    // Also See http://forrest.psych.unc.edu/teaching/p230/Torgerson.pdf
    public static DoubleMatrix2D getDoubleCentered(DoubleMatrix2D input, UnaryOperator<Double> transform) {
        
        final DoubleMatrix2D doubleCenteredInput = input.like();
        if (input instanceof SparseDoubleMatrix2D) {            
            // Sparse implementation
            final DoubleMatrix1D rowSum = getSum(input, transform, true);
            final int[] rowNonZeroCount = getNonZeroCount(input, true);
            final DoubleMatrix1D columnSum = getSum(input, transform, false);
            final int[] columnNonZeroCount = getNonZeroCount(input, false);
            final double matSum = rowSum.zSum();
            int matNonZeroCount = 0;
            for (int count : rowNonZeroCount) {
                matNonZeroCount += count;
            }
            final int matNzc = matNonZeroCount;        
            input.forEachNonZero((int first, int second, double value) -> {
                double d = transform.apply(value);
                // Subtract row average
                d -= rowSum.get(second) / (double) (rowNonZeroCount[second]);
                // Subtract column average
                d -= columnSum.get(first) / (double) (columnNonZeroCount[first]);
                // Add matrix average
                d += matSum / (double) matNzc;
                d /= -2;
                doubleCenteredInput.set(first, second, d);
                return value;
            });
        }
        else {
            // Dense implementation
            DoubleMatrix2D D = DoubleFactory2D.dense.make(input.rows(), input.rows());
            D.assign(input, (double x, double y) -> transform.apply(y));
            DoubleMatrix2D centeringMatrix = DoubleFactory2D.dense.make(input.rows(), input.rows()); 
            DoubleMatrix2D I = DoubleFactory2D.dense.identity(input.rows());
            for (int i = 0; i < input.rows(); ++i)
            {
               for (int j = 0; j < input.rows(); ++j)
               {
                  centeringMatrix.set(i, j, I.get(i, j) - (1.0 / (double) input.rows()));
               }
            }

            // Double center            
            DoubleMatrix2D singleCenteredInput = input.like();
            D.zMult(centeringMatrix, singleCenteredInput);                         
            centeringMatrix.zMult(singleCenteredInput, doubleCenteredInput, -0.5, 0, false, false);
        }
        
        return doubleCenteredInput;
    }    
}
