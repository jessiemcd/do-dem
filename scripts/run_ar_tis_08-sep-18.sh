#!/bin/sh





python ./run_tis_80414201001.py >  tis_out_80414201001.txt &
python ./run_tis_80414202001.py >  tis_out_80414202001.txt &
python ./run_tis_80414203001.py >  tis_out_80414203001.txt &
python ./run_tis_80415201001.py >  tis_out_80415201001.txt &
python ./run_tis_80415202001.py >  tis_out_80415202001.txt &
python ./run_tis_80415203001.py >  tis_out_80415203001.txt &
wait

echo "all orbit scripts finished"
