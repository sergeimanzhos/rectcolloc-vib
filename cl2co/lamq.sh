#!/bin/bash
while kill -0 438 2>/dev/null; do sleep 30; done
for a in "0.25" "0.1" "1.0 3000" "1.0 2000"; do python3 run.py 5000 5000 lam $a; done
