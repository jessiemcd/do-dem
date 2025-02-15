#!/bin/sh





python ./run_tis_20616003001.py >  tis_out_20616003001.txt &
python ./run_tis_20616004001.py >  tis_out_20616004001.txt &
python ./run_tis_20616005001.py >  tis_out_20616005001.txt &
wait

echo "all orbit scripts finished"