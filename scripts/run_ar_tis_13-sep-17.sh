#!/bin/sh





python ./scripts/run_tis_80310231001.py >  tis_out_80310231001.txt &
python ./scripts/run_tis_80310232001.py >  tis_out_80310232001.txt &
wait

echo "all orbit scripts finished"