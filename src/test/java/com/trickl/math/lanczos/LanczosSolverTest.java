package com.trickl.math.lanczos;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import com.trickl.math.lanczos.iteration.LanczosIteration;
import com.trickl.math.lanczos.iteration.LanczosIterationFixed;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
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
        
        double[] eigenvalues = getEigenvalues(data);
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {5, 1}, eigenvalues, 1e-4);
    }    
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testDense2x2MEigenvaluesCase2() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        double[][] data = new double[][] {
        {2, 1},
        {1, 2},        
        };
        
        double[] eigenvalues = getEigenvalues(data);
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {3, 1}, eigenvalues, 1e-4);
    }   
    
    // See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Eigenvalues_and_the_characteristic_polynomial
    @Test
    public void testDense3x3EigenvaluesCase1() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        double[][] data = new double[][] {
        {2, 0, 0},
        {0, 3, 4},        
        {0, 4, 9},        
        };
        
        double[] eigenvalues = getEigenvalues(data);
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {11, 2, 1}, eigenvalues, 1e-4);
    }   
    
    private static double[] getEigenvalues(double[][] data) {
        DoubleMatrix2D mat= DoubleFactory2D.dense.make(data);
        
        LanczosIteration iter = new LanczosIterationFixed(20);
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
}
