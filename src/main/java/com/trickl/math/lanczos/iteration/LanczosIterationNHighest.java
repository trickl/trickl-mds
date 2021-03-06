package com.trickl.math.lanczos.iteration;

import com.trickl.math.lanczos.TridiagonalMatrix;
import java.util.List;

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
public class LanczosIterationNHighest extends LanczosIterationBase {
    private final int n;
    
    public LanczosIterationNHighest(int maxIterations, int n, double relativeTolerance, double absoluteTolerance) {
        super(maxIterations, relativeTolerance, absoluteTolerance);
        this.n = n;
    }
    
    @Override
    public boolean isConverged(TridiagonalMatrix matrix) {
        if(getIteration() > 1) {
            double[] errors = matrix.getErrors();
            double[] eigenvalues = matrix.getEigenvalues();  
        if(eigenvalues.length < n)
          return false;
        else { 
          for(int i = 0; i < n; i++)
            if (errors[errors.length - i - 1] > Math.max(getAbsoluteTolerance(), getRelativeTolerance() * Math.abs(eigenvalues[eigenvalues.length - i - 1])))
	          return false;
	      return true;
	    }
      }
        
      return false;
    }
};