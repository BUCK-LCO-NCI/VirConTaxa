#!/bin/tcsh

touch $2
grep -v "#" $1 >>$2

touch $2

