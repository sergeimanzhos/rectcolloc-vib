#!/bin/bash
# DFT points for the CO2 collocation: 8000 fit points and 800 validation points with V <= 8000 cm-1. Log: pts.log
cd ~/lico2_work/lico2_project/co2_atz; . ../rc/env.sh
python co2_points.py fit 8000 8000 60 && python co2_points.py val 800 8000 60
echo "pts done $? $(date)"
