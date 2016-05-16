package com.trickl.math.lanczos;

import com.trickl.math.lanczos.TridiagonalMatrix;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.PlusMult;
import cern.jet.random.engine.RandomEngine;
import com.trickl.math.Bounds;
import com.trickl.math.lanczos.iteration.LanczosIteration;
import com.trickl.math.lanczos.iteration.LanczosIterationFixed;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.util.Pair;

/*
 * Converted to Java from the IETL library (lanczos.h) 
 * 
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */
public class LanczosSolver extends TridiagonalMatrix {

    public enum ErrorInfo {

        OK, NO_EIGENVALUE, NOT_CALCULATED
    }

    public static class Info {

        private final int[] m1;
        private final int[] m2;
        private final List<Integer> ma;
        private final List<Double> eigenvalue;
        private final List<Double> residuum;
        private final List<ErrorInfo> status;

        public Info(int[] m1, int[] m2, List<Integer> ma,
                List<Double> eigenvalue, List<Double> residuum,
                List<ErrorInfo>  status) {
            this.m1 = m1;
            this.m2 = m2;
            this.ma = ma;
            this.eigenvalue = eigenvalue;
            this.residuum = residuum;
            this.status = status;
        }

        public int m1(int i) {
            return m1[i];
        }

        public int m2(int i) {
            return m2[i];
        }

        public int ma(int i) {
            return ma.get(i);
        }

        public int size() {
            return m1.length;
        }

        public double eigenvalue(int i) {
            return eigenvalue.get(i);
        }

        public double residual(int i) {
            return residuum.get(i);
        }

        public ErrorInfo error_info(int i) {
            return status.get(i);
        }
    }

    private final Algebra alg = new Algebra();
    private final DoubleMatrix2D matrix;
    private DoubleMatrix1D startvector;
    private DoubleMatrix1D vec2;
    int vec2index; 
    int[] M1, M2;
    List<Integer> Ma = new LinkedList<>();

    public LanczosSolver(DoubleMatrix2D matrix) {
        this.matrix = matrix;
        this.startvector = matrix.like1D(Math.max(matrix.rows(), matrix.columns()));
        this.vec2 = matrix.like1D(startvector.size());
        this.vec2index = 0;
    }

    public void calculateEigenvalues(LanczosIteration iter, RandomEngine randomEngine) {
        generateTMatrix(iter, randomEngine);
    }

    public void getMoreEigenvalues(LanczosIteration iter) {
        generateTMatrix(iter);
    }

    public DoubleMatrix2D getEigenvectors(int start, int end, double[] evals, Info info, RandomEngine randomEngine, int maxIterations) {

        DoubleMatrix1D vec3 = startvector.like(startvector.size()); // a temporary vector.
        List<DoubleMatrix1D> eigvectors = new LinkedList<>(); // contains ritz vectors.
        List<List<Double>> Tvectors = new LinkedList<>(); // contains

                                             // eigenvectors of T matrix.
        // calculation of eigen vectors of T matrix(consists of alphas & betas):    
        int n1 = alpha.length;
        double mamax, error, lambda;
        Pair<Double, Double> a_and_b;
        int ma = 0, deltam;
        int nth, maMax = 0, count;
        List<Double> eigenval_a = new LinkedList<>();
        List<Double> residuum = new LinkedList<>();
        List<ErrorInfo> status = new LinkedList<>();

        findM1M2(start, end, evals);
        int M1_itr = 0;
        int M2_itr = 0;

        while (start != end) {

            int maxcount = 10;
            lambda = 0;
            count = 0;
            ErrorInfo errInf = ErrorInfo.OK;

            // calculation of ma starts:
            if (M1[M1_itr] != 0 && M2[M2_itr] != 0) {
                ma = (3 * (M1[M1_itr]) + M2[M2_itr]) / 4 + 1;
                deltam = ((3 * (M1[M1_itr]) + 5 * (M2[M2_itr])) / 8 + 1 - ma) / 10 + 1;
            } else if (M1[M1_itr] != 0 && M2[M2_itr] == 0) {
                ma = (5 * (M1[M1_itr])) / 4 + 1;
                mamax = Math.min((11 * n1) / 8 + 12, (13 * (M1[M1_itr])) / 8 + 1);
                deltam = (int) ((mamax - ma) / 10) + 1;
                if (maxIterations > 0) {
                    maxcount = maxIterations / deltam;
                }
            } else {
                errInf = ErrorInfo.NO_EIGENVALUE;
                deltam = 0;
                ma = 0;
            } // calculation of ma ends.

            eigvectors.add(DoubleFactory1D.dense.make(startvector.size()));
            // new ritz vector is being added in eigvectors.

            List<Double> Tvector = new LinkedList<>();
            Tvectors.add(Tvector);
            // new T matrix vector is being added in Tvectors.

            if (ma == 0) {
               eigvectors.get(eigvectors.size() - 1).assign(0);
            }

            if (ma != 0) {

                double[] eval, z;
                do {
                    if (ma > alpha.length) { // size of T matrix is to be increased.
                        LanczosIteration iter = new LanczosIterationFixed(ma);

                        generateTMatrix(iter);
                    }

                    count++;

                    // on return, z contains all orthonormal eigen vectors of T matrix.
                    EigenDecomposition eigenDecomposition = new EigenDecomposition(alpha, beta);
                    eval = eigenDecomposition.getRealEigenvalues();
                    
                    // search for the value of nth starts, where nth is the nth eigen vector in z.            
                    for (nth = ma - 1; nth >= 0; nth--) {
                        if (Math.abs(eval[nth] - evals[start]) <= thold) {
                            break;
                        }
                    }
                    
                    z = eigenDecomposition.getEigenvector(nth).toArray();

                    // search for the value of ith ends, where ith is the ith eigen vector in z.            
                    if (nth == -1) {
                        error = 0;
                        ma = 0;
                        eigvectors.get(eigvectors.size() - 1).assign(0);
                        errInf = ErrorInfo.NO_EIGENVALUE;
                    } else {
                        error = Math.abs(beta[ma - 1] * z[ma - 1]); // beta[ma - 1] = betaMplusOne.
                        if (error > error_tol) {
                            ma += deltam;
                        }
                    } // end of else
                } while (error > error_tol && count < maxcount);

                if (error > error_tol) {
                    eigvectors.get(eigvectors.size() - 1).assign(0);
                    errInf = ErrorInfo.NOT_CALCULATED;
                } else { // if error is small enough.
                    if (ma != 0) {
                        for (int i = 0; i < ma; i++) {
                            (Tvectors.get(Tvectors.size() - 1)).add(z[i]);
                        }
                        if (ma > maMax) {
                            maMax = ma;
                        }
                        lambda = eval[nth];
                    } // end of if(ma != 0), inner.
                } // end of else{//if error is small enough.
            } // end of if(ma != 0).

            eigenval_a.add(lambda); // for Info object.
            Ma.add(ma); // for Info object.
            status.add(errInf);
            start++;
            M1_itr++;
            M2_itr++;
        } // end of while(in_eigvals_start !=  in_eigvals_end)

      // basis transformation of eigen vectors of T. These vectors are good 
        // approximation of eigen vectors of actual matrix.  
        int eigenvectors_itr;
        int Tvectors_itr;

        a_and_b = makeFirstStep(randomEngine);
        if (a_and_b.getFirst() != alpha[0] || a_and_b.getSecond() != beta[0]) {
            throw new RuntimeException("T-matrix problem at first step");
        }

        eigenvectors_itr = 0;
        Tvectors_itr = 0;
        while (eigenvectors_itr != end) {
            if (!Tvectors.get(Tvectors_itr).isEmpty()) {
                DoubleMatrix1D eigvector = eigvectors.get(eigenvectors_itr);
                eigvector.assign(0);
                eigvector.assign(startvector, PlusMult.plusMult(Tvectors.get(Tvectors_itr).get(0)));
                eigvector.assign(vec2, PlusMult.plusMult(Tvectors.get(Tvectors_itr).get(1)));
            }
            eigenvectors_itr++;
            Tvectors_itr++;
        }
        vec2index = 2;
        for (int j = 2; j < maMax; j++) {
            a_and_b = makeStep(j - 1, vec3);
            if (a_and_b.getFirst() != alpha[j - 1] || a_and_b.getSecond() != beta[j - 1]) {
                throw new RuntimeException("T-matrix problem");
            }

            ++vec2index;
            eigenvectors_itr = 0;
            Tvectors_itr = 0;
            while (eigenvectors_itr != end) {
                if (Tvectors.get(Tvectors_itr).size() > j) {
                    DoubleMatrix1D eigvector = eigvectors.get(eigenvectors_itr);
                    eigvector.assign(vec2, PlusMult.plusMult(Tvectors.get(Tvectors_itr).get(j)));
                }
                // vec2 is being added in one vector of eigvectors.
                eigenvectors_itr++;
                Tvectors_itr++;
            } // end of while loop.      
        } // end of for(int j = 2; j < maMax; j++).
        // end of basis transformation.  // end of basis transformation.  

        // copying to the output iterator & residuum calculation starts:    
        int i = 0;
        DoubleMatrix2D eigenvectors = DoubleFactory2D.dense.make(startvector.size(), end - start);
        
        for (eigenvectors_itr = start;
                eigenvectors_itr != end; eigenvectors_itr++) {
            DoubleMatrix1D eigvector = eigvectors.get(eigenvectors_itr);
            eigenvectors.viewColumn(i).assign(eigvector);
            matrix.zMult(eigvector, vec3);
            vec3.assign(eigvector, PlusMult.minusMult(eigenval_a.get(i)));

            // now vec3 is (A*v - eigenval_a*v); *eigenvectors_itr) is being added in vec3.
            residuum.add(alg.norm2(vec3));
            i++;
        } // copying to the output iterator ends.    
        
        info = new Info(M1, M2, Ma, eigenval_a, residuum, status);
        return eigenvectors;
    }

    private void findM1M2(int start, int end, double[] eigenvalues) {

        int m2counter = 0;
        int n = 1;
        int pos = start;
        M1 = new int[end - start];
        Arrays.fill(M1, 0);
        M2 = new int[end - start];
        Arrays.fill(M2, 0);

        while (m2counter < (end - start) && (n < alpha.length)) {
            n++; // n++ == 2, at first time in this loop.

            EigenDecomposition eigenDecomposition = new EigenDecomposition(alpha, beta);
            double[] eval = eigenDecomposition.getRealEigenvalues();

            int M1_itr = 0;
            int M2_itr = 0;

            while (pos != end) {
                if (M1[M1_itr] == 0 || M2[M2_itr] == 0) {
                    double eigenvalue = eigenvalues[pos];
                    int ub = Bounds.Lower(0, eval.length, index -> eval[index] < eigenvalue + thold);
                    int lb = Bounds.Upper(0, eval.length, index -> eval[index] > eigenvalue - thold);

                    if (M1[M1_itr] == 0 && ub - lb >= 1) {
                        M1[M1_itr] = n;
                    }
                    if (M2[M2_itr] == 0 && ub - lb >= 2) {
                        M2[M2_itr] = n;
                        m2counter++;
                    }
                }
                pos++;
                M1_itr++;
                M2_itr++;
            }
        }
    }

    private void generateTMatrix(LanczosIteration iter, RandomEngine randomEngine) {

        if (alpha.length == 0) {
            Pair<Double, Double> ab = makeFirstStep(randomEngine);
            add(ab);
            vec2index = 1;
        }
        generateTMatrix(iter);
    }

    private void generateTMatrix(LanczosIteration iter) {
        DoubleMatrix1D vec3 = startvector.like(startvector.size());
        for (int j = 0; j < vec2index; j++) {
            iter.next();
        }
        if (alpha.length == 0) {
            throw new RuntimeException("TridiagonalMatrix error, size of TridiagonalMatrix is zero, more_eigenvalues() cannot be called.");
        }

        while (!iter.isFinished(this)) {
            Pair<Double, Double> ab = makeStep(vec2index, vec3);

            if (vec2index == alpha.length) {
                add(ab);
            }
            iter.next();
            ++vec2index;
        }
    }

    private Pair<Double, Double> makeFirstStep(RandomEngine randomEngine) {
        double a, b;
        
        startvector.assign(value -> randomEngine.nextDouble());
        double startvectorMag = alg.norm2(startvector);
        startvector.assign(value -> value / startvectorMag);
        matrix.zMult(startvector, vec2);
        a = startvector.zDotProduct(vec2);
        vec2.assign(startvector, PlusMult.minusDiv(a));
        b = alg.norm2(vec2);
        vec2.assign(value -> value / b);
        return new Pair(a, b);
    }

    private Pair<Double, Double> makeStep(int j, DoubleMatrix1D vec3) {

        matrix.zMult(vec2, vec3);
        double a = vec2.zDotProduct(vec3);
        vec3.assign(vec2, PlusMult.minusDiv(a));
        vec3.assign(startvector, PlusMult.minusDiv(beta[j - 1]));
        double b = alg.norm2(vec3);
        vec3.assign(value -> value / b);
        startvector = vec2;
        vec2 = vec3;

        return new Pair(a, b);
    }
}
