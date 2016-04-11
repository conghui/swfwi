#!/bin/bash

# arguments explanation:
# 1st: NXMAX, total nx, including absorbing boundary and stencil boundary.
#      533 = 461 (actual nx) + 30*2 (absorbing boundary) + 6*2 (stencil boundary)
#
# 2nd: NZMAX, total nz, similar meaning as NXMAX
#      193 = 151 (actual nz) + 30 (absorbing boundary) + 6*2 (stencil boundary)
#
# 3nd: total number of samples (velocities)
#
# 4th: perturbation
#
./genVelPerturb 533 193 100 350
