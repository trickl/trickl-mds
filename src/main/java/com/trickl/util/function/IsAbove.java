package com.trickl.util.function;

import java.util.function.IntPredicate;

public class IsAbove implements IntPredicate {
    private final int cap;

    public IsAbove(int cap) {
        this.cap = cap;
    }

    @Override
    public boolean test(int value) {
        return value > cap;
    }    
}
