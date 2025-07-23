#!/bin/sh





python ./scripts/run_tis_80710201001.py >  tis_out_80710201001.txt &
python ./scripts/run_tis_80710202001.py >  tis_out_80710202001.txt &
python ./scripts/run_tis_80710203001.py >  tis_out_80710203001.txt &
python ./scripts/run_tis_80710204001.py >  tis_out_80710204001.txt &
python ./scripts/run_tis_80710205001.py >  tis_out_80710205001.txt &
python ./scripts/run_tis_80710206001.py >  tis_out_80710206001.txt &
python ./scripts/run_tis_80710207001.py >  tis_out_80710207001.txt &
python ./scripts/run_tis_80710208001.py >  tis_out_80710208001.txt &
wait

echo "all orbit scripts finished"