#!/bin/tcsh

touch $2

grep -v "#" $1 | cut -f1 | sort | uniq >>$2

touch $2
