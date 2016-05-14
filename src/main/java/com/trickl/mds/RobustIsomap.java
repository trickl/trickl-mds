package com.trickl.mds;

import com.trickl.graph.FlowEdge;
import com.trickl.graph.ShortestPaths;
import cern.colt.matrix.DoubleMatrix2D;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Function;
import java.util.function.IntPredicate;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;

// Based on the edge flow idea that allows removing critical outliers for robustness
// Note this implementation does not implement the Kernel feature.
//
// Robust Kernel ISOMAP Algorithm
// Implementation based on:
// Heeyoul Choi & Seungjin Choi 2007
// Department of Computer Science, Pohang University of Science and Technology,
// San 31 Hyoja-dong, Nam-gu, Pohang 790-784, Korea  
// This adds a method of removing critical outliers for robustness
// and a generalisation property by treating the Isomap like a Mercer kernel machine
//
// Uses graph shortest-path distances to discover manifolds defined by closest points
// R is a n x n relational matrix of dissimilarities (assumed to be Euclidean distances)
// p is the dimensionality of the target space 
// k is the size of the nearest neighbourhood considered
// X are the projected points
// m is the p-norm distance measure for the graph, 1 = Manhattan, 2 = Euclidean
// flow_interolance prevents too much flow going through a single node, the lower the intolerance, the more likely for short-circuits
public class RobustIsomap {

    protected static final Logger log = Logger.getLogger(Isomap.class.getCanonicalName());

    private final int k;
    private final int m;
    private final Function<int[], IntPredicate> outlierPredicateFunction;
    private final DoubleMatrix2D S;
    private final DoubleMatrix2D R;

    // TODO: Generalise the methodology for detecting outliers
    public RobustIsomap(DoubleMatrix2D R, int k, int m, Function<int[], IntPredicate> outlierPredicateFunction) {
        if (R.rows() != R.columns()) {
            throw new IllegalArgumentException("Relation matrix must be square.");
        }
        this.k = k;
        this.m = m;
        this.outlierPredicateFunction = outlierPredicateFunction;
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
        // Directed, as we create both an in flow and out flow for each edge (for robustness)
        SimpleDirectedWeightedGraph<Integer, FlowEdge> weightedGraph = new SimpleDirectedWeightedGraph<>(FlowEdge.class);
        for (int i = 0; i < n; ++i) {
            weightedGraph.addVertex(i);
        }

        for (int i = 0; i < n; ++i) {

            // Sort the nearest neighbours to this point
            List<Pair<Integer, Double>> vertexDistanceMap = new ArrayList<>(n - 1);
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    vertexDistanceMap.add(new Pair<>(j, Math.pow(R.get(i, j), m)));
                }
            }

            vertexDistanceMap.sort((Pair<Integer, Double> lhs, Pair<Integer, Double> rhs)
                    -> lhs.getSecond().compareTo(rhs.getSecond())
            );

            int pos = 0; // Add the k nearest neighbours to the graph
            for (Pair<Integer, Double> itr : vertexDistanceMap) {
                if (pos++ >= k) {
                    break;
                }

                // Note directed and no parallel edges, so this allows flow in and out
                FlowEdge edge = weightedGraph.addEdge(i, itr.getKey());
                if (edge != null) {
                    weightedGraph.setEdgeWeight(edge, itr.getValue());
                    edge.zeroFlow();
                }
            }
        }

        // computes all the connected components of the graph
        ConnectivityInspector ci = new ConnectivityInspector(weightedGraph);
        List connectedSets = ci.connectedSets();
        if (connectedSets.size() > 1) {
            throw new IllegalArgumentException("Isomap requires a single connected set.");
        }

        boolean isRobust = false;
        int maxRobustnessIterations = 3;
        for (int robustnessIteration = 0; isRobust == false && robustnessIteration <= maxRobustnessIterations; ++robustnessIteration) {

            // Create a distance and flow map of the shortest path between two nodes
            weightedGraph.vertexSet().stream().forEach((i) -> {
                ShortestPaths<Integer, FlowEdge> shortestPaths = new ShortestPaths<>(weightedGraph, i, (com.trickl.graph.FlowEdge edge) -> {
                    edge.incrementFlow();
                });
                shortestPaths.getDistances().entrySet().stream().forEach((vertexDistance) -> {
                    S.set(i, vertexDistance.getKey(), Math.pow(vertexDistance.getValue(), 1.0 / m));
                });
            });

            // Summarise the edge flows
            final int[] edgeFlows = new int[weightedGraph.edgeSet().size()];
            int edgeIndex = 0;
            for (FlowEdge edge : weightedGraph.edgeSet()) {
                edgeFlows[edgeIndex++] = edge.getFlow();
            }

            // Detect outliers
            IntPredicate isOutlier = outlierPredicateFunction.apply(edgeFlows);

            // Remove outliers, as long as it does not disconnect the graph
            List<FlowEdge> edgesPendingRemoval = new LinkedList<>();
            weightedGraph.edgeSet().stream().forEach((com.trickl.graph.FlowEdge edge) -> {
                if (isOutlier.test(edge.getFlow())) {
                    if (weightedGraph.outDegreeOf(weightedGraph.getEdgeSource(edge)) > 1
                            && weightedGraph.outDegreeOf(weightedGraph.getEdgeTarget(edge)) > 1) {
                        edgesPendingRemoval.add(edge);
                    }
                }
            });

            if (!edgesPendingRemoval.isEmpty()) {
                log.log(Level.INFO, "Iteration:{0} failed robust test.", new Object[]{robustnessIteration});
                // Safe removal without invalidating the iterator
                edgesPendingRemoval.stream().forEach((com.trickl.graph.FlowEdge edge) -> {
                    weightedGraph.removeEdge(edge);
                });

                // Reset the flow for all remaining edges
                edgesPendingRemoval.stream().forEach((com.trickl.graph.FlowEdge edge) -> {
                    edge.zeroFlow();
                });
            } else {
                log.log(Level.INFO, "Iteration:{0} passed robust test.", new Object[]{robustnessIteration});
                isRobust = true;
            }
        }
    }
}
