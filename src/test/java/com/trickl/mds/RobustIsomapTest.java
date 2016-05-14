package com.trickl.mds;

import com.trickl.math.EuclideanDistance;
import com.trickl.math.DistanceMeasure;
import com.trickl.util.function.IsAbove;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.jet.random.engine.MersenneTwister;
import com.trickl.dataset.SwissRoll3D;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import org.junit.Ignore;
import org.junit.Test;
import org.springframework.util.StringUtils;

public class RobustIsomapTest {
    
    public static DoubleMatrix2D createSimpleSquareGrid() {
        
        DoubleMatrix2D data = DoubleFactory2D.dense.make(9, 2); 
        
        //         1.0     1.0
        //      0 ----- 1-------2
        //      |       |       |
        // 0.9  |       | <0.8  | 0.7
        //      3-------4-------5
        //      | ^1.0  | ^1.0  |
        // 0.9  |       | <0.8  | 0.7
        //      6-------7------ 8
        //         1.0     1.0
        
        data.set(0, 0, 0.0);
        data.set(0, 1, 0.0);
        data.set(1, 0, 1.0);
        data.set(1, 1, 0.0);
        data.set(2, 0, 2.0);
        data.set(2, 1, 0.0);
        data.set(3, 0, 0.0);
        data.set(3, 1, 0.9);
        data.set(4, 0, 1.0);
        data.set(4, 1, 0.8);
        data.set(5, 0, 2.0);
        data.set(5, 1, 0.7);
        data.set(6, 0, 0.0);
        data.set(6, 1, 1.8);
        data.set(7, 0, 1.0);
        data.set(7, 1, 1.6);
        data.set(8, 0, 2.0);
        data.set(8, 1, 1.4);
        
        return data;
    }
    
    @Test
    @Ignore
    public void testSimpleSquareGraph() throws IOException {
        DoubleMatrix2D data = createSimpleSquareGrid();
        
        /*
        try (FileWriter writer = new FileWriter("square_graph.dat")) {
            writeToCsv(data, writer);
        }
          */      
        // Convert into similarity data
        DistanceMeasure distanceMeasure = new EuclideanDistance();
        DoubleMatrix2D R = distanceMeasure.getDistances(data);

        RobustIsomap isomap = new RobustIsomap(R, 4, 1, (flows) -> new IsAbove(10));
        DoubleMatrix2D S = isomap.getMappedRelations();

        // Run mapped data through classical mds and project into two dimensions
        ClassicalMds mds = new ClassicalMds(S, 2);
        DoubleMatrix2D X = mds.getReducedSpace();
        /*
        try (FileWriter writer = new FileWriter("square_graph_rkisomap.dat")) {
            writeToCsv(X, writer);
        }
          */      
    }

    @Test
    public void testSwissRoll() throws IOException {
        DoubleMatrix1D normal = new DenseDoubleMatrix1D(3);
        normal.assign(new double[]{.0, .0, 1.0});

        // Use the swiss roll generator to create a swiss roll
        SwissRoll3D swissRoll = new SwissRoll3D();
        swissRoll.setRandomEngine(new MersenneTwister(123456789));
        swissRoll.setNormal(normal);
        swissRoll.setRevolutions(1.5);
        swissRoll.setNoiseStd(0.4);
        swissRoll.setRadius(2);
        swissRoll.setDepth(1);
        DoubleMatrix2D data = swissRoll.generate(250);
        
        // Add an erroneous point connecting the manifold
        data.viewRow(0).assign(new double[] {-0.8, 0.25, 0.5});
        
        /*
        try (FileWriter writer = new FileWriter("swiss_roll_with_noise.dat")) {
            writeToCsv(data, writer);
        }
        */
                
        // Convert into similarity data
        DistanceMeasure distanceMeasure = new EuclideanDistance();
        DoubleMatrix2D R = distanceMeasure.getDistances(data);

        RobustIsomap isomap = new RobustIsomap(R, 5, 2, (flows) -> new IsAbove(250));
        DoubleMatrix2D S = isomap.getMappedRelations();

        // Run mapped data through classical mds and project into two dimensions
        ClassicalMds mds = new ClassicalMds(S, 2);
        DoubleMatrix2D X = mds.getReducedSpace();
        
        /*
        try (FileWriter writer = new FileWriter("swiss_roll_rkisomap.dat")) {
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
