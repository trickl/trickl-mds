package com.trickl.mds;

import com.trickl.math.EuclideanDistance;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import org.apache.commons.math3.random.MersenneTwister;
import com.trickl.dataset.InclinedPlane3D;
import java.awt.Rectangle;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import org.junit.Assert;
import org.junit.Test;
import org.springframework.util.StringUtils;

public class LocalMdsTest {
    
    public LocalMdsTest() {
    }

    @Test
    public void testInclinedPlaneReductionNoNoise() throws IOException {
      DoubleMatrix1D normal = new DenseDoubleMatrix1D(3);
      normal.assign(new double[]{.0, .0, 1.0});

      InclinedPlane3D inclinedPlane = new InclinedPlane3D();
      inclinedPlane.setRandomGenerator(new MersenneTwister(123456789));
      inclinedPlane.setNormal(normal);
      inclinedPlane.setBounds(new Rectangle(-5, -5, 10, 10));
      inclinedPlane.setNoiseStd(0.0);
      DoubleMatrix2D data = inclinedPlane.generate(100);
      
      try (FileWriter writer = new FileWriter("inclined_plane.dat")) {
            writeToCsv(data, writer);
        }

      // Generate similarity data
      EuclideanDistance distance = new EuclideanDistance();
      DoubleMatrix2D R = distance.getDistances(data);
      
      LocalMds mds = new LocalMds(R, 2, 4, new MersenneTwister(123456789));
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
      inclinedPlane.setRandomGenerator(new MersenneTwister(123456789));
      inclinedPlane.setNormal(normal);
      inclinedPlane.setBounds(new Rectangle(-5, -5, 10, 10));
      inclinedPlane.setNoiseStd(0.5);
      DoubleMatrix2D data = inclinedPlane.generate(100);
      
      
        try (FileWriter writer = new FileWriter("inclined_plane_with_noise.dat")) {
            writeToCsv(data, writer);
        }
        

      // Generate similarity data
      EuclideanDistance distance = new EuclideanDistance();
      DoubleMatrix2D R = distance.getDistances(data);
      
      LocalMds mds = new LocalMds(R, 2, 4, new MersenneTwister(123456789));
      DoubleMatrix2D X = mds.getReducedSpace();
      
      // Check the stress is below a threshold
      StressMeasure stressMeasure = new StressMeasure(StressMeasure.Type.Kruskal);
      double stress = stressMeasure.calculate(X, R);
      Assert.assertTrue("The stress is too small given the amount of noise in the data.", stress > 1e-6);
      Assert.assertTrue("The stress is too large given the amount of noise in the data.", stress < 0.5);
   }   
   
   public void writeToCsv(DoubleMatrix2D mat, Writer writer) {
        try (PrintWriter printWriter = new PrintWriter(writer)) {
            for (int i = 0; i < mat.rows(); ++i) {
                StringBuffer row = new StringBuffer();
                List<String> stringValues = new ArrayList<>();
                for (double value : mat.viewRow(i).toArray()) {
                    stringValues.add(String.format("%.5f", value));
                }
                row.append(StringUtils.collectionToCommaDelimitedString(stringValues));
                printWriter.println(row);
            }
        }
    }
}
