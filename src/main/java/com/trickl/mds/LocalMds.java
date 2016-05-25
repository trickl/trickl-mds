package com.trickl.mds;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.PlusMult;
import com.trickl.math.Sorting;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.MathArrays;

// Local MDS, designed for large data sets
// Coded from "A Fast Approximation to MDS"
// Paper by Tynia Yang1, Jinze Liu1, Leonard McMillan1, and Wei Wang1
// University of Chapel Hill at North Carolina, Chapel Hill NC 27599, USA
//
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// W is a n x n weight matrix for each of the dissimlarities
// p is the dimensionality of the target space
// X are the projected points
// s is the number of submatrices to divide R 
public class LocalMds {

    private final int p;
    private final int subdivisions;
    private final DoubleMatrix2D R;
    private final DoubleMatrix2D X;
    private final RandomGenerator randomGenerator;

    protected static final Logger log = Logger.getLogger(LocalMds.class.getCanonicalName());

    public LocalMds(DoubleMatrix2D R, int p, int subdivisions, RandomGenerator randomGenerator) {
        if (R.rows() != R.columns()) {
            throw new IllegalArgumentException("Relation matrix must be square.");
        }
        this.p = p;
        this.subdivisions = subdivisions;
        this.randomGenerator = randomGenerator;
        this.R = R;
        this.X = DoubleFactory2D.dense.make(R.rows(), p);
        
        solve();
    }

    private void solve() {

        int n = R.rows();

        // Divide the dissimilarity matrix into submatrices
        int alignmentSampleSize = p * 2;
        int[] alignmentIndices = new int[subdivisions * alignmentSampleSize];
        int localWidth = (int) (Math.ceil((double) n / (double) subdivisions)); // Must be > sample size

        for (int k = 0; k < subdivisions; ++k) {
            int localStart = (int) Math.floor((double) k * (double) n / (double) subdivisions);

            DoubleMatrix2D D = R.viewPart(localStart, localStart, localWidth, localWidth);

            // Sample  p * 2 points from each submatrix      
            int[] localPoints = new int[localWidth];
            for (int i = 0; i < localWidth; ++i) {
                localPoints[i] = localStart + i;
            }

            MathArrays.shuffle(localPoints, randomGenerator);

            // Just take the first sample size points from the range
            System.arraycopy(localPoints, 0, alignmentIndices, k * alignmentSampleSize, alignmentSampleSize);
        }
        Arrays.sort(alignmentIndices);

        // Construct alignment matrix from the sampled points
        DoubleMatrix2D M = DoubleFactory2D.dense.make(alignmentIndices.length, alignmentIndices.length);

        for (int i = 0; i < alignmentIndices.length; ++i) {
            for (int j = i; j < alignmentIndices.length; ++j) {
                if (i == j) {
                    M.set(i, j, 0);
                } else {

                    double distance = R.get(alignmentIndices[i], alignmentIndices[j]);
                    M.set(i, j, distance);
                    M.set(j, i, distance);
                }
            }
        }

        // Run MDS on the alignment matrix to produce an alignment mapping A
        ClassicalMds alignmentMds = new ClassicalMds(M, p);
        DoubleMatrix2D alignmentSpace = alignmentMds.getReducedSpace();

        // Run MDS on each sub matrix
        for (int k = 0; k < subdivisions; ++k) {
            log.log(Level.INFO, "Processing submatrix {0} of {1}", new Object[]{k + 1, subdivisions});

            int localStart = (int) Math.floor((double) k * (double) n / (double) subdivisions);
            DoubleMatrix2D D = R.viewPart(localStart, localStart, localWidth, localWidth);

            ClassicalMds localMds = new ClassicalMds(D, p);
            DoubleMatrix2D localSpace = localMds.getReducedSpace();

            
            int localAlignmentStart = Sorting.UpperBound(alignmentIndices, localStart - 1);
            int localAlignmentEnd = Sorting.LowerBound(alignmentIndices, localStart + localWidth);
            int[] localAlignmentIndices = new int[localAlignmentEnd - localAlignmentStart];
            System.arraycopy(alignmentIndices, localAlignmentStart, localAlignmentIndices, 0, localAlignmentEnd - localAlignmentStart);

            DoubleMatrix2D localAlignmentSpace = alignmentSpace.viewPart(localAlignmentStart, 0, localAlignmentIndices.length, alignmentSpace.columns());
            DoubleMatrix2D Dmdsi = DoubleFactory2D.dense.make(localAlignmentIndices.length, p + 1);
            {
                int i = 0;
                for (int index : localAlignmentIndices) {
                    double[] Dmdsrow = new double[p + 1];
                    System.arraycopy(localSpace.viewRow(index - localStart).toArray(), 0, Dmdsrow, 0, p);
                    Dmdsrow[p] = 1; // Allow for translation between D and M
                    Dmdsi.viewRow(i).assign(Dmdsrow);
                    i++;
                }
            }

             // Use SVD to find pseudoinverse of D for calculation of affine mapping A
            SingularValueDecomposition svd = new SingularValueDecomposition(Dmdsi);
            DoubleMatrix2D U = svd.getU();
            double[] S = svd.getSingularValues();

            // A = (V * S^-1 * U^T) * Mmds where (U * S * V^T) is the SVD solution of Dmds
            // First calc S^-1 * U^T
            DoubleMatrix2D SUt = DoubleFactory2D.dense.make(Dmdsi.columns(), Dmdsi.rows());
            for (int i = 0; i < Dmdsi.columns(); i++) {
                DoubleMatrix1D Ucol = U.viewColumn(i);
                DoubleMatrix1D SUtrow = SUt.viewRow(i);
                SUtrow.assign(0);
                SUtrow.assign(Ucol, PlusMult.plusDiv(S[i]));
            }

            // Get V from transpose Vt
            DoubleMatrix2D V = svd.getV();

            DoubleMatrix2D Dinv = DoubleFactory2D.dense.make(Dmdsi.columns(), Dmdsi.rows());
            V.zMult(SUt, Dinv);

            DoubleMatrix2D A = DoubleFactory2D.dense.make(Dmdsi.columns(), localAlignmentSpace.columns());
            Dinv.zMult(localAlignmentSpace, A);

            // Last row is translation vector,  store it separately
            DoubleMatrix1D translation = A.viewRow(p);
            A = A.viewPart(0, 0, p, p);

            // Use alignment mapping to transform all points D(i) into aligned M space.
            for (int i = 0; i < localSpace.rows(); i++) {
                DoubleMatrix1D Dmdsrow = localSpace.viewRow(i);
                DoubleMatrix1D Xrow = X.viewRow(localStart + i);

                // TODO: Handle edge overlaps in submatrices by averaging
                A.zMult(Dmdsrow, Xrow, 1.0, 1.0, true);
                Xrow.assign(translation, PlusMult.plusMult(1.0));
            }
        }
    }    
    
    public DoubleMatrix2D getReducedSpace() 
    {
        return X;
    }
}
