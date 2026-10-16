#!/bin/bash
cd ~/lico2_work/lico2_project/rc
while pgrep -f "dft2.py fit --tag atz" > /dev/null; do sleep 60; done
./harm100.sh > harm100.log 2>&1
