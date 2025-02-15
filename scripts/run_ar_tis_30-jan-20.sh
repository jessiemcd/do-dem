#!/bin/sh





python ./scripts/run_tis_20515018001.py >  tis_out_20515018001.txt &
python ./scripts/run_tis_20515021001.py >  tis_out_20515021001.txt &
wait

echo "all orbit scripts finished"