package com.trickl.mds;

import com.trickl.math.CenteringMatrix;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import com.trickl.pca.EigenspaceModel;
import java.util.function.UnaryOperator;

public class SparseMds implements EigenspaceModel {
    
    private final int p;
    private final DoubleMatrix2D R;
    private final DoubleMatrix2D X;
    private final DoubleMatrix2D mean;
    private SingularValueDecomposition svd;
    
    public SparseMds(DoubleMatrix2D R, int p) {
        if (R.rows() != R.columns()) throw new IllegalArgumentException("Relation matrix must be square.");
        this.p = p;
        this.R = R;   
        this.X = DoubleFactory2D.dense.make(R.rows(), p);
        this.mean = X.like(X.rows() > 0 ? 1 : 0, X.columns());
        
        solve();
    }
    
    
    private void solve() {
        
        DoubleMatrix2D doubleCenteredInput = CenteringMatrix.getDoubleCentered(R, (value) -> value * value);
        
/*
   // Use Iterative Eigensolver Template Library to solve for p largest eigenvectors
   typedef ietl::wrapper_vectorspace<boost::numeric::ublas::vector<double> > Vecspace;
   typedef boost::lagged_fibonacci607 Gen; 
   Vecspace vec(D.size1());
   Gen mygen;
   ietl::lanczos<symmetric_matrix<double, lower, column_major>, Vecspace> lanczos(D, vec);

   // Creation of an iteration object:    
   unsigned int max_iter = D.size1();
   double rel_tol = 500 * numeric_limits<double>::epsilon();
   double abs_tol = numeric_limits<double>::epsilon();  

   std::vector<double> eigenvalues;
   std::vector<boost::numeric::ublas::vector<double> > eigenvectors; // for storing the eigen vectors. 
   ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residual, status).
   ietl::lanczos_iteration_nhighest<double> iter(max_iter, max_p, rel_tol, abs_tol);

   try
   {
      lanczos.calculate_eigenvalues(iter,mygen);
      eigenvalues = lanczos.eigenvalues();

      unsigned int p = 0;
      double eigenvalue_total_sum = accumulate(eigenvalues.begin(), eigenvalues.end(), 0.);
      double eigenvalue_partial_sum = 0;
      for (std::vector<double>::const_reverse_iterator itr = eigenvalues.rbegin(), end = eigenvalues.rend();
           itr != end &&  (eigenvalue_partial_sum  / eigenvalue_total_sum) < (1 - variance_tolerance) && p < max_p; ++itr, ++p)
      {
         double eigenvalue = *itr;
         if (eigenvalue < 0) break;
         else eigenvalue_partial_sum += *itr;
      }
       
      lanczos.eigenvectors(eigenvalues.rbegin(), eigenvalues.rbegin() + p, back_inserter(eigenvectors), info, mygen); 
 
      // Can now use double centered matrix to solve for the co-ordinates X = EL^0.5 (XX' = B)
      std::vector<double> sqrt_eigenvalues;
      transform(eigenvalues.rbegin(), eigenvalues.rbegin() + p, back_inserter(sqrt_eigenvalues), ptr_fun<double, double>(sqrt));

      X = zero_matrix<double>(R.size1(), p);
      for (unsigned int j = 0; j < p; ++j)
      {
         matrix_column<matrix<double, column_major> > column(X, j);
         column += eigenvectors[j] * sqrt_eigenvalues[j];
      }
   }
   catch (runtime_error& e)
   {
      cerr << e.what() << endl;
   }
   */
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
