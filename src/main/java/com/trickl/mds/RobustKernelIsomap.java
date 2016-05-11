package com.trickl.mds;

import cern.colt.matrix.DoubleMatrix2D;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.jgrapht.traverse.ClosestFirstIterator;

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
public class RobustKernelIsomap {
    
    public static class FlowEdge extends DefaultWeightedEdge {
        private int flow = 0;
        
        public int getFlow() {
            return flow;
        }
                       
        public void zeroFlow() {
            this.flow = 0;
        }
        
        public void incrementFlow() {
            this.flow += 1;
        }
    }

    protected static final Logger log = Logger.getLogger(Isomap.class.getCanonicalName());

    private final int k;
    private final int m;
    private final double flowIntolerance;
    private final DoubleMatrix2D S;
    private final DoubleMatrix2D R;

    public RobustKernelIsomap(DoubleMatrix2D R, int k, int m, double flowIntolerance) {
        if (R.rows() != R.columns()) {
            throw new IllegalArgumentException("Relation matrix must be square.");
        }
        this.k = k;
        this.m = m;
        this.flowIntolerance = flowIntolerance;
        this.R = R;
        this.S = R.like();

        solve();
    }

    public DoubleMatrix2D getMappedRelations() {
        return S;
    }

    private void solve()  {

        int n = R.rows(); // Number of points

        // Create a neighbourhood graph using all the known dissimilarities
        // Undirected as we assume symmetrical dissimlarities and we need a weight for each edge
        SimpleWeightedGraph<Integer, FlowEdge> weightedGraph = new SimpleWeightedGraph<>(FlowEdge.class);
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
                FlowEdge edge = weightedGraph.addEdge(i, itr.getValue());
                if (edge != null) {
                    weightedGraph.setEdgeWeight(edge, itr.getKey());
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

            // Create a distance matrix of the shortest path between two nodes
            // using Dijkstra's shortest paths            
            for (int i = 0; i < n; ++i) {

                ClosestFirstIterator<Integer, FlowEdge> itr = new ClosestFirstIterator<>(weightedGraph, i, Double.MAX_VALUE);

                while (itr.hasNext()) {
                    Integer j = itr.next();
                    S.set(i, j, Math.pow(itr.getShortestPathLength(j), 1.0 / m));
                    FlowEdge edge = itr.getSpanningTreeEdge(j);
                    edge.incrementFlow();
                }
            }

            // Calculate the distribution of flow across edges
            SummaryStatistics flowSummaryStatistics = new SummaryStatistics();                        
            weightedGraph.edgeSet().stream().forEach((edge) -> {
                flowSummaryStatistics.addValue(edge.getFlow());
            });

            double flow_sigma = flowSummaryStatistics.getStandardDeviation();
            double flow_mean = flowSummaryStatistics.getMean();
            NormalDistribution normalDistribution = new NormalDistribution(flow_mean, flow_sigma);

            List<FlowEdge> edgesPendingRemoval = new LinkedList<>();
            weightedGraph.edgeSet().stream().forEach((edge) -> {
                // This edge has too much flow going through it, check that we don't orphan either vertex by removal
                double cumulative_density = normalDistribution.cumulativeProbability(edge.getFlow());
                if (flowIntolerance > 1 - cumulative_density) {
                    if (weightedGraph.degreeOf(weightedGraph.getEdgeSource(edge)) > 1 && 
                            weightedGraph.degreeOf(weightedGraph.getEdgeTarget(edge)) > 1) {
                        edgesPendingRemoval.add(edge);
                    }
                }
            });

            if (!edgesPendingRemoval.isEmpty()) {
                log.log(Level.INFO, "Iteration:{0} failed robust test.", new Object[] {robustnessIteration});
                // Safe removal without invalidating the iterator
                edgesPendingRemoval.stream().forEach((edge) -> {
                    weightedGraph.removeEdge(edge);
                });

                // Reset the flow for all remaining edges
                edgesPendingRemoval.stream().forEach((edge) -> {
                    edge.zeroFlow();
                });
            } else {
                log.log(Level.INFO, "Iteration:{0} passed robust test.", new Object[] {robustnessIteration});
                isRobust = true;
            }
        }
    }
}
