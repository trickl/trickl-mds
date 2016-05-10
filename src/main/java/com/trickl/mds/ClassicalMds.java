package com.trickl.mds;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import com.trickl.pca.EigenspaceModel;

// Torgerson Classical MDS (1952)
// R is a n x n relational matrix of similarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// X are the projected points
public class ClassicalMds implements EigenspaceModel {
    
    private int p;
    private DoubleMatrix2D R;
    private DoubleMatrix2D X;
    private DoubleMatrix2D mean;
    private SingularValueDecomposition svd;
    
    public ClassicalMds(DoubleMatrix2D R, int p) {
        if (R.rows() != R.columns()) throw new IllegalArgumentException("Relation matrix must be square.");
        this.p = p;
        this.R = R;   
        
        solve();
    }
        
    private void solve() {
        // First need to construct double centred distance matrix B of scalar products
        int n = R.rows();

        // Need d^2 matrix where d is a Euclidean distance
        DoubleMatrix2D D = DoubleFactory2D.dense.make(R.rows(), R.columns());
        for (int i = 0; i < n; ++i)
        {
           for (int j = 0; j < n; ++j)
           {
             D.set(i, j,  Math.pow(R.get(i, j), 2.0));
           }
        }

        // Double centring, so origin is centroid of all points
        DoubleMatrix2D J = DoubleFactory2D.dense.make(n, n); // Double centering matrix
        DoubleMatrix2D I = DoubleFactory2D.dense.identity(n);
        for (int i = 0; i < n; ++i)
        {
           for (int j = 0; j < n; ++j)
           {
              J.set(i, j, I.get(i, j) - (1.0 / (double) n));
           }
        }
        
        DoubleMatrix2D DJ = DoubleFactory2D.dense.make(D.rows(), J.columns());
        D.zMult(J, DJ);
        DoubleMatrix2D B = DoubleFactory2D.dense.make(J.rows(), DJ.columns());
        J.zMult(DJ, B, -0.5, 0, false, false);
        
        // Use svd to find B = ELE'
        svd = new SingularValueDecomposition(B);

        // Can now use B to solve for the co-ordinates X = EL^0.5 (XX' = B)
        X = DoubleFactory2D.dense.make(n, p);
        for (int i = 0; i < p; i++)
        {     
           double l = Math.sqrt(svd.getSingularValues()[i]);
           for (int j = 0; j < n; ++j)
           {
               X.set(j, i, l * svd.getU().get(j, i));
           }
        }
        
        final double weight = 1. /  X.rows();
        mean = X.like(X.rows() > 0 ? 1 : 0, X.columns());
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
        if (svd == null) return null;
        double[] singularValues = svd.getSingularValues();
        DoubleMatrix1D eigenvalues = DoubleFactory1D.dense.make(singularValues.length);
        eigenvalues.assign(singularValues);
        return eigenvalues;
    }

    @Override
    public DoubleMatrix2D getEigenvectors() {
        // Not generally true, SVD is not an eigenvalue decomposition
        if (svd == null) return null;
        return svd.getU(); 
    }

    @Override
    public DoubleMatrix2D getMean() {
        if (svd == null) return null;
        return mean;
    }

    @Override
    public double getWeight() {
        return R.rows();
    }
}