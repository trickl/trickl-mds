package com.trickl.mds;

import com.trickl.graph.ShortestPaths;
import java.util.HashMap;
import java.util.Map;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.junit.Assert;
import org.junit.Test;

public class ShortestPathsTest {
    
    public static SimpleWeightedGraph<Integer, DefaultWeightedEdge> createSimpleSquareGraph() {
        SimpleWeightedGraph<Integer, DefaultWeightedEdge> graph = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        for (int i = 0; i < 9; i++) {
            graph.addVertex(i);
        }
        
        //         1.0     1.0
        //      0 ----- 1-------2
        //      |       |       |
        // 0.9  |       | <0.8  | 0.7
        //      3-------4-------5
        //      | ^1.0  | ^1.0  |
        // 0.9  |       | <0.8  | 0.7
        //      6-------7------ 8
        //         1.0     1.0
        
        graph.addEdge(0, 1);
        graph.setEdgeWeight(graph.getEdge(0, 1), 1.0);
        graph.addEdge(0, 3);
        graph.setEdgeWeight(graph.getEdge(0, 3), 0.9);
        graph.addEdge(1, 2);
        graph.setEdgeWeight(graph.getEdge(1, 2), 1.0);
        graph.addEdge(1, 4);
        graph.setEdgeWeight(graph.getEdge(1, 4), 0.8);
        graph.addEdge(2, 5);
        graph.setEdgeWeight(graph.getEdge(2, 5), 0.7);
        graph.addEdge(3, 4);
        graph.setEdgeWeight(graph.getEdge(3, 4), 1.0);
        graph.addEdge(3, 6);
        graph.setEdgeWeight(graph.getEdge(3, 6), 0.9);
        graph.addEdge(4, 5);
        graph.setEdgeWeight(graph.getEdge(4, 5), 1.0);
        graph.addEdge(4, 7);
        graph.setEdgeWeight(graph.getEdge(4, 7), 0.8);
        graph.addEdge(5, 8);
        graph.setEdgeWeight(graph.getEdge(5, 8), 0.7);
        graph.addEdge(6, 7);
        graph.setEdgeWeight(graph.getEdge(6, 7), 1.0);
        graph.addEdge(7, 8);
        graph.setEdgeWeight(graph.getEdge(7, 8), 1.0);
        
        return graph;
    }
        
    @Test
    public void testGraphDistancesFromCorner() {
                
        SimpleWeightedGraph<Integer, DefaultWeightedEdge> graph = createSimpleSquareGraph();
        ShortestPaths<Integer, DefaultWeightedEdge> shortestPaths = new ShortestPaths<>(graph, 0);
        Map<Integer, Double> shortestDistances = shortestPaths.getDistances();
        Assert.assertEquals("Shortest distances not calculated for all nodes", 9, shortestDistances.size());
        HashMap<Integer, Double> expectedDistances = new HashMap<Integer, Double>() {{                       
           put(0, 0.0);  
           put(1, 1.0);  
           put(2, 2.0);  
           put(3, 0.9);  
           put(4, 1.8);  
           put(5, 2.7);  
           put(6, 1.8);  
           put(7, 2.6);  
           put(8, 3.4); 
        }};
        
        double tolerance = 1e-6;
        graph.vertexSet().stream().forEach((i) -> {
            Assert.assertEquals("Distance from vertex (" + i + ") not as expected",expectedDistances.get(i), shortestDistances.get(i), tolerance);
        });
    }
    
    @Test
    public void testEdgeFlowFromCorner() {
        
        SimpleWeightedGraph<Integer, DefaultWeightedEdge> graph = createSimpleSquareGraph();
        Map<DefaultWeightedEdge, Integer> edgeFlow = new HashMap<>();
        graph.edgeSet().stream().forEach((edge) -> {
            edgeFlow.put(edge, 0);
        });
        ShortestPaths<Integer, DefaultWeightedEdge> shortestPaths = new ShortestPaths<>(graph, 0, Double.MAX_VALUE, Integer.MAX_VALUE, edge -> {
            edgeFlow.put(edge, edgeFlow.get(edge) + 1);
        });
        shortestPaths.getDistances();
        
        HashMap<DefaultWeightedEdge, Integer> expectedFlow = new HashMap<DefaultWeightedEdge, Integer>() {{                       
           put(graph.getEdge(0, 1), 1);  
           put(graph.getEdge(0, 3), 1);  
           put(graph.getEdge(1, 2), 1);  
           put(graph.getEdge(1, 4), 1);  
           put(graph.getEdge(2, 5), 1);  
           put(graph.getEdge(3, 4), 0);  
           put(graph.getEdge(3, 6), 1);  
           put(graph.getEdge(4, 5), 0);  
           put(graph.getEdge(4, 7), 1);  
           put(graph.getEdge(5, 8), 1);  
           put(graph.getEdge(6, 7), 0);  
           put(graph.getEdge(7, 8), 0);  
        }};
        
        graph.edgeSet().stream().forEach((e) -> {
            Assert.assertEquals("Flow on edge (" + e + ") not as expected",expectedFlow.get(e), edgeFlow.get(e));
        });
    }
    
    @Test
    public void testGraphDistancesFromMiddle() {
                
        SimpleWeightedGraph<Integer, DefaultWeightedEdge> graph = createSimpleSquareGraph();
        ShortestPaths<Integer, DefaultWeightedEdge> shortestPaths = new ShortestPaths<>(graph, 4);
        Map<Integer, Double> shortestDistances = shortestPaths.getDistances();
        Assert.assertEquals("Shortest distances not calculated for all nodes", 9, shortestDistances.size());
        HashMap<Integer, Double> expectedDistances = new HashMap<Integer, Double>() {{
           put(0, 1.8);  
           put(1, 0.8);  
           put(2, 1.7);  
           put(3, 1.0);  
           put(4, 0.0);  
           put(5, 1.0);  
           put(6, 1.8);  
           put(7, 0.8);  
           put(8, 1.7);              
        }};
        
        double tolerance = 1e-6;
        graph.vertexSet().stream().forEach((i) -> {
            Assert.assertEquals("Distance from vertex (" + i + ") not as expected",expectedDistances.get(i), shortestDistances.get(i), tolerance);
        });
    }
    
}
