package com.trickl.mds;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import com.vividsolutions.jts.algorithm.ConvexHull;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import org.apache.commons.lang.NotImplementedException;

/**
 * Defined as sqrt(area) / (raw stress + 2 * mean(distance))
 * A stress measure than penalises folding
 * @author tgee
 */
public class PlanarityMeasure {
    
    public enum Type {
        Convex,
        Concave
    }
    
    public PlanarityMeasure(Type type) {
        this.type = type;
    }
    
    private final Type type;
    
    public double calculate(DoubleMatrix2D X, DoubleMatrix2D R) {
       DoubleMatrix2D W = R.like(R.rows(), R.columns());
       R.forEachNonZero((int i, int j, double value) -> {
            W.set(i, j, 1.0);
            return value;
       });
       return calculate(X, R, W);
    }

    public double calculate(DoubleMatrix2D X, DoubleMatrix2D R, DoubleMatrix2D W) {

        if (X.columns() != 2) {
            throw new IllegalArgumentException("Planarity can only be calculated for 2D embeddings");
        }
        
        StressMeasure rawStressMeasure = new StressMeasure(StressMeasure.Type.Raw);
        GeometryFactory geometryFactory = new GeometryFactory();
        Geometry hull;
        if (type == Type.Concave) {
            // TODO: Note trickl graphs has a good Voronoi implementation already. This is half of the work of
            // of concave hull algorithm.
            throw new NotImplementedException("Still to implement concave stress, requires a concave hull algoritm.");
        }
        else {
            final Coordinate[] coords =new Coordinate[X.rows()];
            for (int i = 0; i < X.rows(); i++) {
                coords[i] = new Coordinate(X.get(i, 0), X.get(i, 1));
            }
            ConvexHull convexHull = new ConvexHull(coords, geometryFactory);                       
            hull = convexHull.getConvexHull();
        }
        
        double rawStress = rawStressMeasure.calculate(X, R, W);
        
        double area = hull.getArea();
        
        double meanDistance = R.zSum() / (double) R.cardinality();
        
        return Math.sqrt(area) / (rawStress + (2 * meanDistance));
    }
}
