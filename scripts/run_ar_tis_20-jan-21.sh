#!/bin/sh





python ./scripts/run_tis_20614002001.py >  tis_out_20614002001.txt &
python ./scripts/run_tis_20614003001.py >  tis_out_20614003001.txt &
wait

echo "all orbit scripts finished"