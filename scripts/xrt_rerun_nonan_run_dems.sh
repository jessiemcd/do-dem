#!/bin/sh



python ./scripts/dems_01-sep-15.py >  ./scripts/dems_out_01-sep-15.txt &
python ./scripts/dems_11-sep-17.py >  ./scripts/dems_out_11-sep-17.txt &
python ./scripts/dems_13-sep-17.py >  ./scripts/dems_out_13-sep-17.txt &
python ./scripts/dems_29-may-18_1.py >  ./scripts/dems_out_29-may-18_1.txt &
python ./scripts/dems_29-apr-21.py >  ./scripts/dems_out_29-apr-21.txt &
python ./scripts/dems_30-jul-21_1.py >  ./scripts/dems_out_30-jul-21_1.txt &
wait

echo "all orbit scripts finished"