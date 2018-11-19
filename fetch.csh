#!/bin/tcsh

touch $3

foreach i (`cat $1`)
samtools faidx $2 $i >>! $3
end


