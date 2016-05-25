package com.trickl.graph;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Consumer;
import org.jgrapht.WeightedGraph;
import org.jgrapht.traverse.ClosestFirstIterator;

public class ShortestPaths<V, E> {
        
    private final WeightedGraph<V, E> graph;
    private final V startVertex;
    private final double maxDistance;
    private final int maxVertices;
    private final Consumer<E> onEdge;
    
    public ShortestPaths(WeightedGraph<V, E> graph, V startVertex) {
        this(graph, startVertex, Double.MAX_VALUE, Integer.MAX_VALUE);
    }
    
    public ShortestPaths(WeightedGraph<V, E> graph, V startVertex, double maxDistance, int maxVertices) {
        this(graph, startVertex, maxDistance, maxVertices, null);
    }
    
    public ShortestPaths(WeightedGraph<V, E> graph, V startVertex, double maxDistance, int maxVertices, Consumer<E> onEdge) {
        this.graph = graph;
        this.startVertex = startVertex;
        this.maxDistance = maxDistance;
        this.maxVertices = maxVertices;
        this.onEdge = onEdge;
    }

    public Map<V, Double> getDistances() {
        // Create a distance matrix of the shortest path between two nodes
        // using Dijkstra's shortest paths
        Map<V, Double> distances = new HashMap<>();
        ClosestFirstIterator<V, E> itr = new ClosestFirstIterator<>(graph, startVertex, maxDistance);
        int vertices = 0;
        while (itr.hasNext() && vertices++ < maxVertices) {
            V vj = itr.next();
            distances.put(vj, itr.getShortestPathLength(vj));

            if (onEdge != null) {
                E edge = itr.getSpanningTreeEdge(vj);
                if (edge != null) {
                    onEdge.accept(edge);
                }
            }
        }

        return distances;
    }
}
