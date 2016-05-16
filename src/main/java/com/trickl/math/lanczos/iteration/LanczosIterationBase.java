package com.trickl.math.lanczos.iteration;

import com.trickl.math.lanczos.TridiagonalMatrix;

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
abstract public class LanczosIterationBase extends IterationBase implements LanczosIteration {
    
        
    public LanczosIterationBase(int maxIterations, double relativeTolerance, double absoluteTolerance)
    {
        super(maxIterations, relativeTolerance, absoluteTolerance);
    }
    
    @Override
    public boolean isFinished(TridiagonalMatrix tmatrix) {
      if (isConverged(tmatrix))
	return true;
      else if (iteration < maxIterations)
	return false;
      else {
	fail (1, "maximum number of iterations exceeded");
	return true;
      }
    }
        
    @Override
    abstract public boolean isConverged(TridiagonalMatrix tmatrix);    
  };
    