package com.trickl.math;

import java.util.function.IntPredicate;

public class Bounds {

    public static int Lower(int first, int last, IntPredicate isLower) {
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

    public static int Upper(int first, int last, IntPredicate isHigher) {
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
