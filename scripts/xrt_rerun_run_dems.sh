#!/bin/sh



python ./scripts/dems_01-sep-15.py >  dems_out_01-sep-15.txt &
python ./scripts/dems_19-feb-16.py >  dems_out_19-feb-16.txt &
python ./scripts/dems_22-apr-16_2.py >  dems_out_22-apr-16_2.txt &
python ./scripts/dems_11-sep-17.py >  dems_out_11-sep-17.txt &
python ./scripts/dems_12-sep-17.py >  dems_out_12-sep-17.txt &
python ./scripts/dems_13-sep-17.py >  dems_out_13-sep-17.txt &
python ./scripts/dems_29-may-18_1.py >  dems_out_29-may-18_1.txt &
python ./scripts/dems_12-apr-19.py >  dems_out_12-apr-19.txt &
python ./scripts/dems_13-apr-19.py >  dems_out_13-apr-19.txt &
python ./scripts/dems_29-jan-20.py >  dems_out_29-jan-20.txt &
python ./scripts/dems_29-apr-21.py >  dems_out_29-apr-21.txt &
python ./scripts/dems_03-may-21_2.py >  dems_out_03-may-21_2.txt &
python ./scripts/dems_30-jul-21_1.py >  dems_out_30-jul-21_1.txt &
wait

echo "all orbit scripts finished"