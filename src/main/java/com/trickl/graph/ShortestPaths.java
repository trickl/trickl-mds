package com.trickl.graph;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Consumer;
import org.jgrapht.WeightedGraph;
import org.jgrapht.traverse.ClosestFirstIterator;

public class ShortestPaths<V, E> {
        
    private final WeightedGraph<V, E> graph;
    private final V startVertex;
    private final Consumer<E> onEdge;
    
    public ShortestPaths(WeightedGraph<V, E> graph, V startVertex) {
        this(graph, startVertex, null);
    }
    
    public ShortestPaths(WeightedGraph<V, E> graph, V startVertex, Consumer<E> onEdge) {
        this.graph = graph;
        this.startVertex = startVertex;
        this.onEdge = onEdge;
    }

    public Map<V, Double> getDistances() {
        // Create a distance matrix of the shortest path between two nodes
        // using Dijkstra's shortest paths
        Map<V, Double> distances = new HashMap<>();
        ClosestFirstIterator<V, E> itr = new ClosestFirstIterator<>(graph, startVertex, Double.MAX_VALUE);
        while (itr.hasNext()) {
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
