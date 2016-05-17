package com.trickl.math.lanczos;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import com.trickl.math.lanczos.iteration.LanczosIteration;
import com.trickl.math.lanczos.iteration.LanczosIterationFixed;
import java.util.List;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

public class LanczosSolverTest {
    
    @Test
    public void testDenseSmallMatrixEigenvalues() {
        
        // Chosen for simplicity, has eigenvalues 1 and 5
        double[][] data = new double[][] {
        {3, 2},
        {2, 3},        
        };
        
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
        
        // Check converged to correct eigenvalues
        Assert.assertArrayEquals("Eigenvalues not as expected", new double[] {5, 1}, eigenvalues, 1e-4);
    }    
}
