package com.trickl.mds;

import com.trickl.math.CenteringMatrix;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import com.trickl.math.lanczos.LanczosSolver;
import com.trickl.math.lanczos.iteration.LanczosIteration;
import com.trickl.math.lanczos.iteration.LanczosIterationNHighest;
import com.trickl.pca.EigenspaceModel;
import java.util.Arrays;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;

public class SparseMds implements EigenspaceModel {
    
    private final int p;
    private final DoubleMatrix2D R;
    private final DoubleMatrix2D X;
    
    private final DoubleMatrix2D mean;
    private DoubleMatrix2D eigenvectors;
    private double[] eigenvalues;
    private double relTolerance;
    private double absTolerance;
    
    public SparseMds(DoubleMatrix2D R, int p) {
       this (R, p, 1e-8, 1e-10);
    }
    
    public SparseMds(DoubleMatrix2D R, int p, double relTolerance, double absTolerance) {
        if (R.rows() != R.columns()) throw new IllegalArgumentException("Relation matrix must be square.");
        this.p = p;
        this.R = R;   
        this.relTolerance = relTolerance;
        this.absTolerance = absTolerance;
        this.X = DoubleFactory2D.dense.make(R.rows(), p);
        this.mean = X.like(X.rows() > 0 ? 1 : 0, X.columns());
        
        solve();
    }    
    
    private void solve() {
        // First need to construct double centred distance matrix B of scalar products
        int n = R.rows();
        
        DoubleMatrix2D B = CenteringMatrix.getDoubleCentered(R, (value) -> value * value);
        
        int maxIterations = R.rows();
        LanczosIteration iter = new LanczosIterationNHighest(maxIterations, p, relTolerance, absTolerance);
        RandomGenerator randomGenerator = new MersenneTwister(123456789);
        LanczosSolver solver = new LanczosSolver(B, absTolerance);    
        solver.calculateEigenvalues(iter, randomGenerator);
        double[] allEigenvalues = solver.getEigenvalues();
        eigenvalues = Arrays.copyOfRange(allEigenvalues, allEigenvalues.length - p, allEigenvalues.length);
        randomGenerator = new MersenneTwister(123456789);
        eigenvectors = solver.getEigenvectors(eigenvalues, randomGenerator, maxIterations);   
        
        // Can now use B to solve for the co-ordinates X = EL^0.5 (XX' = B)        
        for (int i = 0; i < p; i++)
        {     
           double l = Math.sqrt(eigenvalues[i]);
           for (int j = 0; j < n; ++j)
           {
               X.set(j, i, l * eigenvectors.get(j, i));
           }
        }
        
        final double weight = 1. /  X.rows();        
        X.forEachNonZero((int first, int second, double value) -> {
            mean.setQuick(0, second, mean.getQuick(0, second) + value * weight);
            return value;
        });
        
    }
        
    public DoubleMatrix2D getReducedSpace() 
    {
        return X;
    }

    @Override
    public DoubleMatrix1D getEigenvalues() {
        // Not generally true, SVD is not an eigenvalue decomposition        
        if (eigenvalues == null) return null;
        DoubleMatrix1D eigenvaluesMatrix = DoubleFactory1D.dense.make(eigenvalues.length);
        eigenvaluesMatrix.assign(eigenvalues);
        return eigenvaluesMatrix;
    }

    @Override
    public DoubleMatrix2D getEigenvectors() {
        return eigenvectors; 
    }

    @Override
    public DoubleMatrix2D getMean() {
        if (eigenvectors == null) return null;
        return mean;
    }

    @Override
    public double getWeight() {
        return R.rows();
    }
}
