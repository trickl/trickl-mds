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
public class LanczosIterationNLowest extends LanczosIterationBase {
    private final int n;
    
    public LanczosIterationNLowest(int maxIterations, int n, double relativeTolerance, double absoluteTolerance) {
        super(maxIterations, relativeTolerance, absoluteTolerance);
        this.n = n;
    }
    
    @Override
    public boolean isConverged(TridiagonalMatrix matrix) {
        if(getIteration() > 1) {
            List<Double> errors = matrix.errors();
            List<Double> eigenvalues = matrix.eigenvalues();  
        if(eigenvalues.size() < n)
          return false;
        else { 
          for(int i = 0; i < n; i++)
            if (errors.get(i) > Math.max(getAbsoluteTolerance(), getRelativeTolerance() * Math.abs(eigenvalues.get(i))))
	          return false;
	      return true;
	    }
      }
      return false;
    }
};    