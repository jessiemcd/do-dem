#!/bin/sh





python ./scripts/run_tis_90311211001.py >  tis_out_90311211001.txt &
python ./scripts/run_tis_90311212001.py >  tis_out_90311212001.txt &
python ./scripts/run_tis_90311213001.py >  tis_out_90311213001.txt &
wait

echo "all orbit scripts finished"