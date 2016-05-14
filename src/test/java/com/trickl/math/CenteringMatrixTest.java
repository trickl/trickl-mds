package com.trickl.math;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import java.util.function.UnaryOperator;
import org.junit.Assert;
import org.junit.Test;

public class CenteringMatrixTest {
    
    @Test
    public void testDenseDoubleCenteringNoTransform() {        
        DoubleMatrix2D input = DoubleFactory2D.dense.make(new double[][] {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
        });
        DoubleMatrix2D result = CenteringMatrix.getDoubleCentered(input);        
        assertIsCentered(result);
    }
    
    @Test
    public void testSparseDoubleCenteringNoTransform() {        
        DoubleMatrix2D input = DoubleFactory2D.sparse.make(new double[][] {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
        });
        DoubleMatrix2D result = CenteringMatrix.getDoubleCentered(input);        
        assertIsCentered(result);
    }
    
    @Test
    public void testDenseDoubleCenteringSquareTransform() {        
        DoubleMatrix2D input = DoubleFactory2D.dense.make(new double[][] {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
        });
        DoubleMatrix2D result = CenteringMatrix.getDoubleCentered(input, value -> value * value);        
        assertIsCentered(result);
    }
    
    @Test
    public void testSparseDoubleCenteringSquareTransform() {        
        DoubleMatrix2D input = DoubleFactory2D.sparse.make(new double[][] {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
        });
        DoubleMatrix2D result = CenteringMatrix.getDoubleCentered(input, value -> value * value);        
        assertIsCentered(result);
    }
    
    public void assertIsCentered(DoubleMatrix2D mat) {
        double tolerance = 1e-6;
        for (int row = 0; row < mat.rows(); row++) {
            Assert.assertEquals("Sum of row " + row + " is not zero.", 0, mat.viewRow(row).zSum(), tolerance);
        }
        for (int column = 0; column < mat.columns(); column++) {
            Assert.assertEquals("Sum of column " + column + " is not zero.", 0, mat.viewColumn(column).zSum(), tolerance);
        }
    }
}
