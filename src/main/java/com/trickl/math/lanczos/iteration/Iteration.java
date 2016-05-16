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
public interface Iteration {

    boolean isFirst();

    double getAbsoluteTolerance();

    int getErrorCode();

    int getMaxIterations();

    double getRelativeTolerance();

    boolean isFinished(double r, double lambda);

    int getIteration();

    void next();   
}
