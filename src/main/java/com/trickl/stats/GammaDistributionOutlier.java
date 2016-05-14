/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.trickl.stats;

import com.trickl.util.function.ValueAboveHasProbabilityBelow;
import java.util.function.Function;
import java.util.function.IntPredicate;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

/**
 *
 * @author tgee
 */
public class GammaDistributionOutlier implements Function<int[], IntPredicate> {
    private final double probability;

    public GammaDistributionOutlier(double probability) {
        this.probability = probability;
    }

    @Override
    public IntPredicate apply(int[] edgeFlows) {
        // Calculate the distribution of flow across edges
        SummaryStatistics flowSummaryStatistics = new SummaryStatistics();
        for (int flow : edgeFlows) {
            flowSummaryStatistics.addValue(flow);
        }
        double flowVar = flowSummaryStatistics.getVariance();
        double flowMean = flowSummaryStatistics.getMean();
        double gammaShape = (flowMean * flowMean) / flowVar;
        double gammaScale = flowVar / flowMean;
        GammaDistribution gammaDistribution = new GammaDistribution(gammaShape, gammaScale);
        return new ValueAboveHasProbabilityBelow(gammaDistribution, probability);
    }
    
}
