package com.trickl.mds;

import cern.colt.matrix.DoubleFactory2D;
import com.trickl.graph.ShortestPaths;
import cern.colt.matrix.DoubleMatrix2D;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.IntPredicate;
import java.util.function.UnaryOperator;
import java.util.logging.Logger;
import org.jgrapht.DirectedGraph;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.WeightedGraph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.jgrapht.traverse.ClosestFirstIterator;

// J. B. Tenenbaum, V. de Silva, and J. C. Langford ISOMAP Algorithm
// C++ Version built using Josh Tenenbaum's MatLab version as reference
//
// Uses graph shortest-path distances to discover manifolds defined by closest points
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// k is the size of the nearest neighbourhood considered
// m is the p-norm distance measure for the graph, 1 = Manhattan, 2 = Euclidean
// X are the projected points
public class Isomap {

    protected static final Logger log = Logger.getLogger(Isomap.class.getCanonicalName());
    
    protected final int m;    
    protected final DoubleMatrix2D R;
    protected DoubleMatrix2D S;
    
    // Only calculate full graph distances for the indexes where this evaluates to true
    // This can reduce computational complexity drastically.
    protected final IntPredicate samplePredicate;
    
    public Isomap(DoubleMatrix2D R, int m) {
        this(R, m, index -> true); // Sample every vertex
    }

    public Isomap(DoubleMatrix2D R, int m, IntPredicate samplePredicate) {
        if (R.rows() != R.columns()) {
            throw new IllegalArgumentException("Relation matrix must be square.");
        }
        this.m = m;
        this.R = R;
        this.samplePredicate = samplePredicate;
    }
    
    public DoubleMatrix2D getMappedRelations() {
        if (S == null) {
            S = R.like();
            solve();
        }
        return S;
    }
    
    public static DoubleMatrix2D kNearest(DoubleMatrix2D R, int k) {
        // Preprocess
        Map<Integer, TreeMap<Double, Integer>> vertexDistances = new HashMap<>();
        R.forEachNonZero((i, j, value) -> {
            if (!vertexDistances.containsKey(i)) {
                vertexDistances.put(i, new TreeMap<>());
            }
            
            TreeMap<Double, Integer> distanceMap = vertexDistances.get(i);
            distanceMap.put(value, j);

            return value;
        });
        
        DoubleMatrix2D R2 = DoubleFactory2D.sparse.make(R.rows(), R.columns());
        for (Map.Entry<Integer, TreeMap<Double, Integer>> vertexDistanceEntry : vertexDistances.entrySet()) {
            int i = vertexDistanceEntry.getKey();
            int pos = 0;
            for (Map.Entry<Double, Integer> distances : vertexDistanceEntry.getValue().entrySet()) {
                if (pos++ >= k) {
                    break;
                }
                int j = distances.getValue();
                R2.set(i, j, distances.getKey());
            }
        }
        
        return R2;
    };
    
    protected <E> void connectGraph(WeightedGraph<Integer, E> weightedGraph) {
        
        int n = R.rows(); // Number of points
        
        // Create a neighbourhood graph using all the known dissimilarities
        // Undirected as we assume symmetrical dissimlarities and we need a weight for each edge
        for (int i = 0; i < n; ++i) {
            weightedGraph.addVertex(i);
        }
        
        R.forEachNonZero((i, j, value) -> {
            // Note undirected and no parallel edges, so this allows flow in and out         
            E edge = weightedGraph.addEdge(i, j);
            if (edge != null) {
                weightedGraph.setEdgeWeight(edge, Math.pow(value, m));
            }
            return value;
        });              
    }
    
    protected void checkIsSinglyConnected(ConnectivityInspector ci) {
        List connectedSets = ci.connectedSets();
        if (connectedSets.size() > 1) {
            throw new IllegalArgumentException("Isomap requires a single connected set.");
        }
    }
    
    protected static <E> void populateDistanceMatrixFromGraph(WeightedGraph<Integer, E> weightedGraph, DoubleMatrix2D S, 
            UnaryOperator<Double> distanceTransform,
            IntPredicate samplePredicate) {
        weightedGraph.vertexSet().stream().forEach((i) -> {
            if (samplePredicate.test(i)) {
                ShortestPaths<Integer, E> shortestPaths = new ShortestPaths<>(weightedGraph, i);
                shortestPaths.getDistances().entrySet().stream().forEach((vertexDistance) -> {
                    S.set(i, vertexDistance.getKey(), distanceTransform.apply(vertexDistance.getValue()));
                });
            }
        });
    }

    private void solve() {

        int n = R.rows(); // Number of points
        
        // Create a graph
        SimpleWeightedGraph<Integer, DefaultWeightedEdge> weightedGraph = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        connectGraph(weightedGraph);
        checkIsSinglyConnected(new ConnectivityInspector(weightedGraph));

        // Create a distance matrix of the shortest path between two nodes
        populateDistanceMatrixFromGraph(weightedGraph, S, value -> Math.pow(value, 1.0 / m), samplePredicate);
    }
}
