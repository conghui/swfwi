## SWC++ compiler issues
- we can only pass the `-O0` optimized flags to the `mpicxx -ver 5.421-sw-437 -host`. Compile it with `-O2` will creash the program.

- the follow codes in `src/common/parabola-vertex.cpp`

```
  if (std::abs(k2 - k1) < 0.001 * (std::max(std::abs(k2), std::abs(k1))) ||
      (xv == -std::numeric_limits<double>::quiet_NaN())) {
    WARNING() << "THE SET OF POINTS DON'T FIT PARABOLIC WELL, SET y TO y3 ON PURPOSE
    xv = x3;
    yv = y3;
  }
```


