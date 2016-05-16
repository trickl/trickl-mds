package com.trickl.math.lanczos;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;
import com.trickl.math.lanczos.iteration.LanczosIteration;
import com.trickl.math.lanczos.iteration.LanczosIterationFixed;
import java.util.List;
import org.junit.Assert;
import org.junit.Test;

public class LanczosSolverTest {
    
    @Test
    public void testDenseSmallMatrixEigenvalues() {
        // Copied from https://en.wikipedia.org/wiki/Singular_value_decomposition#Example
        double[][] data = new double[][] {
        {1, 0, 0, 0, 2},
        {0, 0, 3, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 2, 0, 0, 0},        
        };
        
        DoubleMatrix2D mat= DoubleFactory2D.dense.make(data);
        
        LanczosIteration iter = new LanczosIterationFixed(20);
        RandomEngine randomEngine = new MersenneTwister(123456789);
        LanczosSolver solver = new LanczosSolver(mat);
        solver.calculateEigenvalues(iter, randomEngine);
        double[] eigenvalues = solver.getEigenvalues();
        
        // Check convergence
        double errorTolerance = 1e-8;
        for (double error: solver.getErrors()) {
            Assert.assertTrue("Error " + error + " is larger than tolerance " + errorTolerance, error < errorTolerance);
        }
        
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {2, 3, 2.2361}, eigenvalues, 1e-4);
    }    
}
