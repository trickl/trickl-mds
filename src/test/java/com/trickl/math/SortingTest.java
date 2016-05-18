package com.trickl.math;

import org.junit.Test;
import static org.junit.Assert.*;

public class SortingTest {
    
    public SortingTest() {
    }

    @Test
    public void testLowerBounds() {
        int[] sortedInts = new int[] {10, 10, 10, 20, 20, 20, 30, 30};
        int lowerBound = Sorting.LowerBound(sortedInts, 20);
        assertEquals("Lower bound not as expected", 3, lowerBound);
    }

    @Test
    public void testUpperBounds() {
        int[] sortedInts = new int[] {10, 10, 10, 20, 20, 20, 30, 30};
        int upperBound = Sorting.UpperBound(sortedInts, 20);
        assertEquals("Upper bound not as expected", 6, upperBound);
    }
    
}
