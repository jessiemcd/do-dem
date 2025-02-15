#!/bin/sh





python ./scripts/run_tis_80414201001.py >  tis_out_80414201001.txt &
python ./scripts/run_tis_80414202001.py >  tis_out_80414202001.txt &
python ./scripts/run_tis_80414203001.py >  tis_out_80414203001.txt &
wait

echo "all orbit scripts finished"