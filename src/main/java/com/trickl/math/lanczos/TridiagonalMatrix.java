package com.trickl.math.lanczos;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import org.apache.commons.lang.ArrayUtils;
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
public class TridiagonalMatrix {

    protected static final double epsilon = 1e-16;
    protected double[] alpha = new double[0];
    protected double[] beta = new double[0];
    protected final double errorTolerance;
    protected double threshold;

    private boolean computed;
    private double multipleTolerance;
    protected final List<Double> err = new LinkedList<>();
    protected final List<Double> err_noghost = new LinkedList<>();
    protected final List<Double> eigval_distinct = new LinkedList<>(); // distinct eigen values.
    protected final List<Double> eigval_distinct_noghost = new LinkedList<>(); // distinct eigen values 
    private final List<Integer> multiplicty = new LinkedList<>();
    private final List<Integer> multiplicty_noghost = new LinkedList<>();
    private double alpha_max;
    private double beta_max;
    private double beta_min;
    
    public TridiagonalMatrix() {
        this(Math.pow(epsilon, 2. / 3.));
    }

    public TridiagonalMatrix(double errorTolerance) {
        this.errorTolerance = errorTolerance;
    }

    public double[] getEigenvalues() {
        return getEigenvalues(true);
    }

    public double[] getEigenvalues(boolean discard_ghosts) {
        if (!computed) {
            compute();
        }
        
        List<Double> eigenvalueList;
        if (discard_ghosts) {
            eigenvalueList = eigval_distinct_noghost;
        } else {
            eigenvalueList =  eigval_distinct;
        }
        
        return ArrayUtils.toPrimitive(eigenvalueList.toArray(new Double[eigenvalueList.size()]));
    }

    public double[] getErrors() {
        return getErrors(true);
    }

    public double[] getErrors(boolean discard_ghosts) {
        if (!computed) {
            compute();
        }
        
        List<Double> errorList;
        if (discard_ghosts) {
            errorList = err_noghost;
        } else {
            errorList = err;
        }
        
        return ArrayUtils.toPrimitive(errorList.toArray(new Double[errorList.size()]));
    }

    protected List<Integer> getMultiplicities() {
        return getMultiplicities(true);
    }

    protected List<Integer> getMultiplicities(boolean discard_ghosts) {
        if (!computed) {
            compute();
        }
        if (discard_ghosts) {
            return multiplicty_noghost;
        } else {
            return multiplicty;
        }
    }
    
    
    protected void add(Pair<Double, Double> a_and_b) {
        add(a_and_b.getFirst(), a_and_b.getSecond());
    }

    protected void add(double a, double b) {
        computed = false;

        double[] alphaTmp = alpha;
        alpha = new double[alphaTmp.length + 1];
        System.arraycopy(alphaTmp, 0, alpha, 0, alphaTmp.length);
        alpha[alphaTmp.length] = a;

        double[] betaTmp = beta;
        beta = new double[betaTmp.length + 1];
        System.arraycopy(betaTmp, 0, beta, 0, betaTmp.length);
        beta[betaTmp.length] = b;

        if (alpha.length == 1) {
            alpha_max = a;
            beta_min = beta_max = b;
        } else {
            if (a > alpha_max) {
                alpha_max = a;
            }
            if (b > beta_max) {
                beta_max = b;
            }
            if (b < beta_min) {
                beta_min = b;
            }
        }
    }

    protected void compute() {
        err.clear();
        eigval_distinct.clear();
        multiplicty.clear();

        err_noghost.clear();
        eigval_distinct_noghost.clear();
        multiplicty_noghost.clear();

        computed = true;
        int n = alpha.length;
        EigenDecomposition eigenDecomposition = new EigenDecomposition(alpha, beta);
        double[] eval = eigenDecomposition.getRealEigenvalues();
        Arrays.sort(eval); // Consistent with IETL
    
        // tolerance values:
        multipleTolerance = Math.max(alpha_max, beta_max) * 2 * epsilon * (1000 + n);
        threshold = Math.max(eval[0], eval[n - 1]);
        threshold = Math.max(errorTolerance * threshold, 5 * multipleTolerance);

    // error estimates of eigen values starts:    
        // the unique eigen values selection, their multiplicities and corresponding errors calculation follows:
        double temp = eval[0];
        eigval_distinct.add(eval[0]);
        int multiple = 1;

        for (int i = 1; i < n; i++) {
            double[] eigenvector = eigenDecomposition.getEigenvector(eval.length - i).toArray();
            if (Math.abs(eval[i] - temp) > threshold) {
                eigval_distinct.add(eval[i]);
                temp = eval[i];
                multiplicty.add(multiple);
                if (multiple > 1) {
                    err.add(0.);
                } else {
                    err.add(Math.abs(beta[beta.length - 1] * eigenvector[n - 1])); // *beta.rbegin() = betaMplusOne.
                }
                multiple = 1;
            } else {
                multiple++;
            }
        }

        // for last eigen value.
        multiplicty.add(multiple);
        if (multiple > 1) {
            err.add(.0);
        } else {
            double[] eigenvector = eigenDecomposition.getEigenvector(eval.length - n).toArray();
            err.add(Math.abs(beta[beta.length - 1] * eigenvector[n - 1])); // *beta.rbegin() = betaMplusOne.
        }
        // the unique eigen values selection, their multiplicities and corresponding errors calculation ends.
        // ghosts calculations starts:
        double[] beta_g = Arrays.copyOfRange(beta, 1, beta.length);
        double[] alpha_g = Arrays.copyOfRange(alpha, 1, alpha.length);

        eigenDecomposition = new EigenDecomposition(alpha_g, beta_g);
        double[] eval_g = eigenDecomposition.getRealEigenvalues();
        Arrays.sort(eval_g); // Consistent with IETL
    
        int i = 0, t2 = 0;

        for (double eigval : eigval_distinct) {            
            if (multiplicty.get(i) == 1) { // test of spuriousness for the eigenvalues whose multiplicity is one.
                for (int j = t2; j < n - 1; j++, t2++) { // since size of reduced matrix is n-1
                    if (eval_g[j] - eigval >= multipleTolerance) {
                        break;
                    }

                    if (Math.abs(eigval - eval_g[j]) < multipleTolerance) {
                        multiplicty.set(i, 0);
                        err.set(i, .0); // if eigen value is a ghost => error calculation not required, 0=> ignore error.
                        t2++;
                        break;
                    }
                }
            }
            i++;
        }

        i = 0;
        for (double eigval : eigval_distinct) {            
            if (multiplicty.get(i) != 0) {
                eigval_distinct_noghost.add(eigval);
                multiplicty_noghost.add(multiplicty.get(i));
                err_noghost.add(err.get(i));
            }
            i++;
        }
    }
}
