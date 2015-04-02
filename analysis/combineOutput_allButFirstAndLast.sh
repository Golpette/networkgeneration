#!/bin/bash

## Script to combine multiple data files but ommitting first and last lines of file.
## (First line in our comparison output labels the columns, while the last line can me corrupted/incomplete if we cut simulations short)


ls z_*.dat | while read file; 

  do

      cat $file | head -n -1 | tail -n +2 >> an_ParameterSampling.dat
      
  done;
