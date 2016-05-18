package com.trickl.math;

import java.util.function.IntPredicate;

public class Sorting {
    
    public static int LowerBound(double[] arr, double value) {
       return LowerBound(0, arr.length, (index) -> arr[index] < value);
    }
    
    public static int LowerBound(int[] arr, int value) {
       return LowerBound(0, arr.length, (index) -> arr[index] < value);
    }
    
    public static <T extends Comparable, Object> int LowerBound(T[] arr, T value) {
       return LowerBound(0, arr.length, (index) -> arr[index].compareTo(value) < 0);
    }

    public static int LowerBound(int first, int last, IntPredicate isLower) {
        if (isLower == null) {
            throw new NullPointerException();
        }
        int len = last - first;
        while (len > 0) {
            int half = len / 2;
            int middle = first + half;
            if (isLower.test(middle)) {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        return first;
    }
    
    public static int UpperBound(double[] arr, double value) {
       return UpperBound(0, arr.length, (index) -> arr[index] > value);
    }
    
    public static int UpperBound(int[] arr, int value) {
       return UpperBound(0, arr.length, (index) -> arr[index] > value);
    }
    
    public static <T extends Comparable> int UpperBound(T[] arr, T value) {
       return UpperBound(0, arr.length, (index) -> arr[index].compareTo(value) > 0);
    }

    public static int UpperBound(int first, int last, IntPredicate isHigher) {
        if (isHigher == null) {
            throw new NullPointerException();
        }
        int len = last - first;
        while (len > 0) {
            int half = len / 2;
            int middle = first + half;
            if (isHigher.test(middle)) {
                len = half;
            } else {
                first = middle + 1;
                len -= half + 1;
            }
        }
        return first;
    }
}
