#!/bin/sh





python ./scripts/run_tis_90710201001.py >  tis_out_90710201001.txt &
python ./scripts/run_tis_90710203001.py >  tis_out_90710203001.txt &
wait

echo "all orbit scripts finished"