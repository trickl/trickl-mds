package com.trickl.math.lanczos;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import com.trickl.math.lanczos.LanczosSolver.ErrorInfo;
import com.trickl.math.lanczos.LanczosSolver.Info;
import com.trickl.math.lanczos.iteration.LanczosIteration;
import com.trickl.math.lanczos.iteration.LanczosIterationFixed;
import com.trickl.matrixunit.MatrixAssert;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.Pair;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

public class LanczosSolverTest {
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testDense2x2EigenvaluesCase1() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        double[][] data = new double[][] {
        {3, 2},
        {2, 3},        
        };
        
        double[] eigenvalues = getEigenvalues(DoubleFactory2D.dense.make(data));
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {1, 5}, eigenvalues, 1e-4);
    }    
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testDense2x2EigenvectorsCase1() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        double[][] data = new double[][] {
        {3, 2},
        {2, 3},        
        };
        
        DoubleMatrix2D eigenvectors = getEigenvectors(DoubleFactory2D.dense.make(data));
        
        double[][] expectedEigenvectors = new double[][] {
        {Math.sqrt(2)/2, Math.sqrt(2)/2},
        {-Math.sqrt(2)/2, Math.sqrt(2)/2},        
        };
        
        MatrixAssert.assertEquals(DoubleFactory2D.dense.make(expectedEigenvectors), eigenvectors, 1e-4);
    }  
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testDense2x2MEigenvaluesCase2() {
        
        // Chosen for simplicity, has eigenvalues 1 and 3
        double[][] data = new double[][] {
        {2, 1},
        {1, 2},        
        };
        
        double[] eigenvalues = getEigenvalues(DoubleFactory2D.dense.make(data));
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {1, 3}, eigenvalues, 1e-4);
    }   
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testDense2x2EigenvectorsCase2() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        double[][] data = new double[][] {
        {2, 1},
        {1, 2},        
        };
        
        DoubleMatrix2D eigenvectors = getEigenvectors(DoubleFactory2D.dense.make(data));
        
        double[][] expectedEigenvectors = new double[][] {
        {Math.sqrt(2)/2, Math.sqrt(2)/2},
        {-Math.sqrt(2)/2, Math.sqrt(2)/2},        
        };
        
        MatrixAssert.assertEquals(DoubleFactory2D.dense.make(expectedEigenvectors), eigenvectors, 1e-4);
    } 
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testDense3x3EigenvaluesCase1() {
        
        // Chosen for simplicity, has eigenvalues 11, 2, 1
        double[][] data = new double[][] {
        {2, 0, 0},
        {0, 3, 4},        
        {0, 4, 9},        
        };
        
        double[] eigenvalues = getEigenvalues(DoubleFactory2D.dense.make(data));
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {1, 2, 11}, eigenvalues, 1e-4);
    }   
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testDense3x3igenvectorsCase1() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        double[][] data = new double[][] {
        {2, 0, 0},
        {0, 3, 4},        
        {0, 4, 9},        
        };
        
        DoubleMatrix2D eigenvectors = getEigenvectors(DoubleFactory2D.dense.make(data));
        
        double[][] expectedEigenvectors = new double[][] {
        {0,                 -1,      0},
        {2 / Math.sqrt(5),  0,      1 / Math.sqrt(5)},        
        {-1 / Math.sqrt(5), 0,      2 / Math.sqrt(5)},       
        };
        
        MatrixAssert.assertEquals(DoubleFactory2D.dense.make(expectedEigenvectors), eigenvectors, 1e-4);
    } 
           
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testSparse3x3EigenvaluesCase1() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        DoubleMatrix2D matrix = DoubleFactory2D.sparse.make(3, 3);
        matrix.set(0, 0, 2);
        matrix.set(1, 1, 3);
        matrix.set(2, 2, 9);
        matrix.set(1, 2, 4);
        matrix.set(2, 1, 4);
                
        double[] eigenvalues = getEigenvalues(matrix);
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {1, 2, 11}, eigenvalues, 1e-4);
    }   
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testSparse3x3EigenvectorsCase1() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        DoubleMatrix2D matrix = DoubleFactory2D.sparse.make(3, 3);
        matrix.set(0, 0, 2);
        matrix.set(1, 1, 3);
        matrix.set(2, 2, 9);
        matrix.set(1, 2, 4);
        matrix.set(2, 1, 4);
        
        DoubleMatrix2D eigenvectors = getEigenvectors(matrix);
        
        double[][] expectedEigenvectors = new double[][] {
        {0,                 -1,      0},
        {2 / Math.sqrt(5),  0,      1 / Math.sqrt(5)},        
        {-1 / Math.sqrt(5), 0,      2 / Math.sqrt(5)},       
        };
        
        MatrixAssert.assertEquals(DoubleFactory2D.dense.make(expectedEigenvectors), eigenvectors, 1e-4);
    } 
    
    private static double[] getEigenvalues(DoubleMatrix2D mat) {
        
        int maxIterations = 20;
        LanczosIteration iter = new LanczosIterationFixed(maxIterations);
        RandomGenerator randomGenerator = new MersenneTwister(123456789);
        LanczosSolver solver = new LanczosSolver(mat);
        solver.calculateEigenvalues(iter, randomGenerator);
        double[] eigenvalues = solver.getEigenvalues();
        
        // Check convergence
        double errorTolerance = 1e-8;
        for (double error: solver.getErrors()) {
            Assert.assertTrue("Error " + error + " is larger than tolerance " + errorTolerance, error < errorTolerance);
        }
        return eigenvalues;
    }
    
    private static DoubleMatrix2D getEigenvectors(DoubleMatrix2D mat) {
        
        int maxIterations = 20;
        LanczosIteration iter = new LanczosIterationFixed(maxIterations);
        RandomGenerator randomGenerator = new MersenneTwister(123456789);
        LanczosSolver solver = new LanczosSolver(mat);    
        solver.calculateEigenvalues(iter, randomGenerator);
        double[] eigenvalues = solver.getEigenvalues();
        randomGenerator = new MersenneTwister(123456789);
        Pair<DoubleMatrix2D, Info> eigenvectorsWithInfo = solver.getEigenvectorsWithInfo(eigenvalues, randomGenerator, maxIterations);          
        Info eigenvectorInfo = eigenvectorsWithInfo.getSecond();
        
        for (int i = 0; i < eigenvectorInfo.size(); ++i) {
            Assert.assertEquals("Position " + i + " has unexpected result " + eigenvectorInfo.error_info(i), ErrorInfo.OK, eigenvectorInfo.error_info(i));
        }
        
        return eigenvectorsWithInfo.getFirst();
    }
}
