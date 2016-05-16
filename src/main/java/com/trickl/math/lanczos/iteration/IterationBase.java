package com.trickl.math.lanczos.iteration;

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
public class IterationBase implements Iteration {
    
    protected final int maxIterations;
    protected final double relativeTolerance;
    protected final double absoluteTolerance;
    protected int iteration;
    protected int errorCode;        
    protected String errorMessage;

    public IterationBase(int maxIterations, double relativeTolerance, double absoluteTolerance)
    {
        this.iteration = 0;
        this.errorCode = 0;
        this.maxIterations = maxIterations;
        this.relativeTolerance = relativeTolerance;
        this.absoluteTolerance = absoluteTolerance;
    }

    @Override
    public boolean isFinished(double r, double lambda) {
        if (isConverged(r, lambda)) {
            return true;
        } else if (iteration < maxIterations) {
            return false;
        } else {
            fail(1, "maximum number of iterations exceeded");
            return true;
        }
    }

    boolean isConverged(double r, double lambda) {
        return r <= relativeTolerance * Math.abs(lambda) || r < absoluteTolerance; // relative or absolute tolerance.
    }

    @Override
    public void next() {
        ++iteration;
    }

    @Override
    public boolean isFirst() {
        return iteration == 0;
    }

    @Override
    public int getErrorCode() {
        return errorCode;
    }

    @Override
    public int getIteration() {
        return iteration;
    }

    @Override
    public double getRelativeTolerance() {
        return relativeTolerance;
    }

    @Override
    public double getAbsoluteTolerance() {
        return absoluteTolerance;
    }

    @Override
    public int getMaxIterations() {
        return maxIterations;
    }

    protected void fail(int errorCode) {
        this.errorCode = errorCode;
    }

    protected void fail(int errorCode, String errorMessage) {
        this.errorCode = errorCode;
        this.errorMessage = errorMessage;
    }
    
}
