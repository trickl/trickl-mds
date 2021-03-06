package com.trickl.mds;

import com.trickl.math.EuclideanDistance;
import com.trickl.math.DistanceMeasure;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import com.trickl.dataset.SwissRoll3D;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.random.MersenneTwister;
import org.junit.Test;
import org.springframework.util.StringUtils;

public class IsomapTest {

    @Test
    public void testSwissRoll() throws IOException {
        DoubleMatrix1D normal = new DenseDoubleMatrix1D(3);
        normal.assign(new double[]{.0, .0, 1.0});

        // Use the swiss roll generator to create a swiss roll
        SwissRoll3D swissRoll = new SwissRoll3D();
        swissRoll.setRandomGenerator(new MersenneTwister(123456789));
        swissRoll.setNormal(normal);
        swissRoll.setRevolutions(1.5);
        swissRoll.setNoiseStd(0);
        swissRoll.setRadius(2);
        swissRoll.setDepth(1);
        DoubleMatrix2D data = swissRoll.generate(250);
       /* 
        try (FileWriter writer = new FileWriter("swiss_roll.dat")) {
            writeToCsv(data, writer);
        }
         */       

        // Convert into similarity data
        DistanceMeasure distanceMeasure = new EuclideanDistance();
        DoubleMatrix2D R = distanceMeasure.getDistances(data);

        Isomap isomap = new Isomap(R, 2);
        DoubleMatrix2D S = isomap.getMappedRelations();

        // Run mapped data through classical mds and project into two dimensions
        ClassicalMds mds = new ClassicalMds(S, 2);
        DoubleMatrix2D X = mds.getReducedSpace();
        /*
        try (FileWriter writer = new FileWriter("swiss_roll_isomap.dat")) {
            writeToCsv(X, writer);
        }
          */       
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
