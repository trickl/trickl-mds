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
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;
import org.jgrapht.graph.SimpleWeightedGraph;

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
public class RobustIsomap extends Isomap {

    private final Function<int[], IntPredicate> outlierPredicateFunction;

    public RobustIsomap(DoubleMatrix2D R, int m, Function<int[], IntPredicate> outlierPredicateFunction) {
        this(R, m, outlierPredicateFunction, index -> true);
    }
    
    public RobustIsomap(DoubleMatrix2D R, int m, Function<int[], IntPredicate> outlierPredicateFunction, IntPredicate samplePredicate) {
        super(R, m, samplePredicate);
        this.outlierPredicateFunction = outlierPredicateFunction;
    }

    private void solve() {

        int n = R.rows(); // Number of points

        // Create a neighbourhood graph using all the known dissimilarities        
        // Directed, as we create both an in flow and out flow for each edge (for robustness)
        SimpleDirectedWeightedGraph<Integer, FlowEdge> weightedGraph = new SimpleDirectedWeightedGraph<>(FlowEdge.class);
        connectGraph(weightedGraph);
        checkIsSinglyConnected(new ConnectivityInspector(weightedGraph));

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
