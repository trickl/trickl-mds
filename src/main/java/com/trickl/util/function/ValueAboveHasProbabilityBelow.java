/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.trickl.util.function;

import java.util.function.IntPredicate;
import org.apache.commons.math3.distribution.RealDistribution;

/**
 *
 * @author tgee
 */
public class ValueAboveHasProbabilityBelow implements IntPredicate {
    private final RealDistribution distribution;
    private final double probability;

    public ValueAboveHasProbabilityBelow(RealDistribution distribution, double probability) {
        this.distribution = distribution;
        this.probability = probability;
    }

    @Override
    public boolean test(int value) {
        double cumulativeProbability = distribution.cumulativeProbability(value);
        return (1 - cumulativeProbability) < probability;
    }
    
}
