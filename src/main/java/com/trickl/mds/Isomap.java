package com.trickl.mds;

import com.trickl.graph.ShortestPaths;
import cern.colt.matrix.DoubleMatrix2D;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;
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

    private final int k;
    private final int m;
    private final DoubleMatrix2D S;
    private final DoubleMatrix2D R;

    public Isomap(DoubleMatrix2D R, int k, int m) {
        if (R.rows() != R.columns()) {
            throw new IllegalArgumentException("Relation matrix must be square.");
        }
        this.k = k;
        this.m = m;
        this.R = R;
        this.S = R.like();

        solve();
    }
    
    public DoubleMatrix2D getMappedRelations() {
        return S;
    }

    private void solve() {

        int n = R.rows(); // Number of points

        // Create a neighbourhood graph using all the known dissimilarities
        // Undirected as we assume symmetrical dissimlarities and we need a weight for each edge
        SimpleWeightedGraph<Integer, DefaultWeightedEdge> weightedGraph = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        for (int i = 0; i < n; ++i) {
            weightedGraph.addVertex(i);
        }

        for (int i = 0; i < n; ++i) {

            // Sort the nearest neighbours to this point
            TreeMap<Double, Integer> distanceVertexMap = new TreeMap<>();
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    distanceVertexMap.put(Math.pow(R.get(i, j), m), j);
                }
            }

            int pos = 0; // Add the k nearest neighbours to the graph
            for (Map.Entry<Double, Integer> itr : distanceVertexMap.entrySet()) {
                if (pos++ >= k) {
                    break;
                }

                // Note undirected and no parallel edges, so this allows flow in and out         
                DefaultWeightedEdge edge = weightedGraph.addEdge(i, itr.getValue());
                if (edge != null) {
                    weightedGraph.setEdgeWeight(edge, itr.getKey());
                }
            }
        }

        // computes all the connected components of the graph
        ConnectivityInspector ci = new ConnectivityInspector(weightedGraph);
        List connectedSets = ci.connectedSets();
        if (connectedSets.size() > 1) {
            throw new IllegalArgumentException("Isomap requires a single connected set.");
        }

        // Create a distance matrix of the shortest path between two nodes
        weightedGraph.vertexSet().stream().forEach((i) -> {
            ShortestPaths<Integer, DefaultWeightedEdge> shortestPaths = new ShortestPaths<>(weightedGraph, i);
            shortestPaths.getDistances().entrySet().stream().forEach((vertexDistance) -> {
                S.set(i, vertexDistance.getKey(), Math.pow(vertexDistance.getValue(), 1.0 / m));
            });
        });
    }
}
