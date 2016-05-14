/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.trickl.graph;

import org.jgrapht.graph.DefaultWeightedEdge;

/**
 *
 * @author tgee
 */
public class FlowEdge extends DefaultWeightedEdge {
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

    @Override
    public String toString() {
        return super.toString() + " [ " + flow + "]";
    }
    
}
