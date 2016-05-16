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
public class LanczosIterationFixed extends LanczosIterationBase {
       
    public LanczosIterationFixed(int maxIterations) {
        super(maxIterations, 0, 0);
    }
    
    @Override
    public boolean isConverged(TridiagonalMatrix matrix) {
      return false;
    }
};    