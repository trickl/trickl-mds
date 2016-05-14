package com.trickl.mds;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import com.trickl.pca.EigenspaceModel;

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
        /*
    }
   // Double centring, so origin is centroid of all points
   // Avoid centering matrix as it is slow and inefficient
   // Iterator only traverses one side of the diagonal, hence semi sums
   std::vector<pair<double, unsigned int> > semi_rowsums(R.size2()); // < sum, count >
   std::vector<pair<double, unsigned int> > semi_colsums(R.size1()); // < sum, count >
   fill(semi_rowsums.begin(), semi_rowsums.end(), pair<double, unsigned int>(0., 0));
   fill(semi_colsums.begin(), semi_colsums.end(), pair<double, unsigned int>(0., 0));

   for (symmetric_matrix<double, lower, column_major>::const_iterator1 i_itr = R.begin1(), i_end = R.end1(); i_itr != i_end; ++i_itr)
   {
      for (symmetric_matrix<double, lower, column_major>::const_iterator2 j_itr = i_itr.begin(), j_end = i_itr.end(); j_itr != j_end; ++j_itr)
      {
         unsigned int i = i_itr.index1(), j = j_itr.index2(); 
         if (i >= j) // Include the diagonal, lower triangle
         {  
            double value = R(i, j) * R(i, j);
            semi_rowsums[i].first += value;
            semi_rowsums[i].second++;
            semi_colsums[j].first += value;
            semi_colsums[j].second++;
         }
      } 
   }
   double semi_matsum = 0;
   unsigned int semi_matsize = 0; // Effective for sparse matrices
   for (std::vector<pair<double, unsigned int> >::const_iterator itr = semi_rowsums.begin(), end = semi_rowsums.end(); itr != end; itr++) 
   {
      semi_matsum += itr->first;
      semi_matsize += itr->second;
   }
   double matsum = semi_matsum * 2.;
   unsigned int matsize = semi_matsize * 2 - R.size1();

   symmetric_matrix<double, lower, column_major> D(R.size1());
   for (symmetric_matrix<double, lower, column_major>::iterator1 i_itr = D.begin1(), i_end = D.end1(); i_itr != i_end; ++i_itr)
   {
      for (symmetric_matrix<double, lower, column_major>::iterator2 j_itr = i_itr.begin(), j_end = i_itr.end(); j_itr != j_end; ++j_itr)
      {
         unsigned int i = i_itr.index1(), j = j_itr.index2(); 
         if (i >= j) // Include the diagonal, lower triangle
         { 
            // Distance squared
            double d = R(i, j);
            d *= d;

            // Subtract row average
            d -= ((semi_rowsums[i].first + semi_colsums[i].first) /
                       double(semi_rowsums[i].second + semi_colsums[i].second - 1));
            // Subtract column average
            d -= ((semi_colsums[j].first + semi_rowsums[j].first) /
                       double(semi_colsums[j].second + semi_rowsums[j].second - 1));
            // Add matrix average
            d += matsum / matsize;
            d /= -2.;
            D(i, j) = d;
         } 
      }
   }

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
