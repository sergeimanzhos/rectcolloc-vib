#!/bin/bash
cd ~/lico2_work/lico2_project/rc; . ./env.sh
python tune100.py --workers 5
echo "tune100 done $(date)"
