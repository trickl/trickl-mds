package com.trickl.mds;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.jet.random.engine.MersenneTwister;
import com.trickl.dataset.InclinedPlane3D;
import java.awt.Rectangle;
import java.io.IOException;
import org.junit.Assert;
import org.junit.Test;

public class ClassicalMdsTest {
    
    public ClassicalMdsTest() {
    }

    @Test
   public void testInclinedPlaneReductionNoNoise() throws IOException {
      DoubleMatrix1D normal = new DenseDoubleMatrix1D(3);
      normal.assign(new double[]{.0, .0, 1.0});

      InclinedPlane3D inclinedPlane = new InclinedPlane3D();
      inclinedPlane.setRandomEngine(new MersenneTwister(123456789));
      inclinedPlane.setNormal(normal);
      inclinedPlane.setBounds(new Rectangle(-5, -5, 10, 10));
      inclinedPlane.setNoiseStd(0.0);
      DoubleMatrix2D data = inclinedPlane.generate(100);

      // Generate similarity data
      EuclideanDistance distance = new EuclideanDistance();
      DoubleMatrix2D R = distance.getDistances(data);
      
      ClassicalMds mds = new ClassicalMds(R, 2);
      DoubleMatrix2D X = mds.getReducedSpace();
      
      // Check the stress is below a threshold
      StressMeasure stressMeasure = new StressMeasure(StressMeasure.Type.Kruskal);
      double stress = stressMeasure.calculate(X, R);
      Assert.assertTrue("The stress is too large given no noise in the data.", stress < 1e-6);
   }   
   
   @Test
   public void testInclinedPlaneReductionSmallNoise() throws IOException {
      DoubleMatrix1D normal = new DenseDoubleMatrix1D(3);
      normal.assign(new double[]{.0, .0, 1.0});

      InclinedPlane3D inclinedPlane = new InclinedPlane3D();
      inclinedPlane.setRandomEngine(new MersenneTwister(123456789));
      inclinedPlane.setNormal(normal);
      inclinedPlane.setBounds(new Rectangle(-5, -5, 10, 10));
      inclinedPlane.setNoiseStd(0.5);
      DoubleMatrix2D data = inclinedPlane.generate(100);

      // Generate similarity data
      EuclideanDistance distance = new EuclideanDistance();
      DoubleMatrix2D R = distance.getDistances(data);
      
      ClassicalMds mds = new ClassicalMds(R, 2);
      DoubleMatrix2D X = mds.getReducedSpace();
      
      // Check the stress is below a threshold
      StressMeasure stressMeasure = new StressMeasure(StressMeasure.Type.Kruskal);
      double stress = stressMeasure.calculate(X, R);
      Assert.assertTrue("The stress is too small given the amount of noise in the data.", stress > 1e-6);
      Assert.assertTrue("The stress is too large given the amount of noise in the data.", stress < 0.5);
   }   
}
