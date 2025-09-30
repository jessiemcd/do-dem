#!/bin/sh



python ./scripts/dems_01-sep-15.py >  ./scripts/dems_out_01-sep-15.txt &
python ./scripts/dems_02-sep-15.py >  ./scripts/dems_out_02-sep-15.txt &
python ./scripts/dems_19-feb-16.py >  ./scripts/dems_out_19-feb-16.txt &
python ./scripts/dems_22-apr-16_1.py >  ./scripts/dems_out_22-apr-16_1.txt &
python ./scripts/dems_22-apr-16_2.py >  ./scripts/dems_out_22-apr-16_2.txt &
python ./scripts/dems_26-jul-16_1.py >  ./scripts/dems_out_26-jul-16_1.txt &
python ./scripts/dems_27-jul-16_1.py >  ./scripts/dems_out_27-jul-16_1.txt &
python ./scripts/dems_26-jul-16_2.py >  ./scripts/dems_out_26-jul-16_2.txt &
python ./scripts/dems_11-sep-17.py >  ./scripts/dems_out_11-sep-17.txt &
python ./scripts/dems_12-sep-17.py >  ./scripts/dems_out_12-sep-17.txt &
python ./scripts/dems_13-sep-17.py >  ./scripts/dems_out_13-sep-17.txt &
python ./scripts/dems_10-oct-17.py >  ./scripts/dems_out_10-oct-17.txt &
python ./scripts/dems_29-may-18_1.py >  ./scripts/dems_out_29-may-18_1.txt &
python ./scripts/dems_29-may-18_2.py >  ./scripts/dems_out_29-may-18_2.txt &
python ./scripts/dems_09-sep-18.py >  ./scripts/dems_out_09-sep-18.txt &
python ./scripts/dems_10-sep-18.py >  ./scripts/dems_out_10-sep-18.txt &
python ./scripts/dems_12-apr-19.py >  ./scripts/dems_out_12-apr-19.txt &
python ./scripts/dems_13-apr-19.py >  ./scripts/dems_out_13-apr-19.txt &
python ./scripts/dems_29-jan-20.py >  ./scripts/dems_out_29-jan-20.txt &
python ./scripts/dems_06-jun-20.py >  ./scripts/dems_out_06-jun-20.txt &
python ./scripts/dems_07-jun-20.py >  ./scripts/dems_out_07-jun-20.txt &
python ./scripts/dems_08-jun-20.py >  ./scripts/dems_out_08-jun-20.txt &
python ./scripts/dems_09-jun-20.py >  ./scripts/dems_out_09-jun-20.txt &
python ./scripts/dems_08-jan-21.py >  ./scripts/dems_out_08-jan-21.txt &
python ./scripts/dems_20-jan-21.py >  ./scripts/dems_out_20-jan-21.txt &
python ./scripts/dems_29-apr-21.py >  ./scripts/dems_out_29-apr-21.txt &
python ./scripts/dems_03-may-21_1.py >  ./scripts/dems_out_03-may-21_1.txt &
python ./scripts/dems_03-may-21_2.py >  ./scripts/dems_out_03-may-21_2.txt &
python ./scripts/dems_07-may-21.py >  ./scripts/dems_out_07-may-21.txt &
python ./scripts/dems_20-jul-21.py >  ./scripts/dems_out_20-jul-21.txt &
python ./scripts/dems_30-jul-21_1.py >  ./scripts/dems_out_30-jul-21_1.txt &
python ./scripts/dems_30-jul-21_2.py >  ./scripts/dems_out_30-jul-21_2.txt &
wait

echo "all orbit scripts finished"