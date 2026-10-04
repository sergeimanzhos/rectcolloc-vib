#!/bin/bash
while kill -0 217 2>/dev/null; do sleep 20; done
python3 run.py 5000 5000 cfg 2 0.95 1150 "2500,5000" t_f095
python3 run.py 5000 5000 cfg 2 1.05 1150 "2500,5000" t_f105
python3 run.py 5000 5000 cfg 2 1.0 900 "2500,5000" t_ec900
python3 run.py 5000 5000 cfg 2 1.0 1400 "2500,5000" t_ec1400
python3 run.py 5000 5000 cfg 3 1.0 1150 "2500,5000" t_p3
