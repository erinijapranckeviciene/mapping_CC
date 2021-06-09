./mapping_CC -s $1 -m $2 -w 1 -e 0 | gawk '{sum=0; for(i=2;i<=NF;i++) sum+=$i; print $1, sum/(NF-1);}' | gawk -f chr_max_pos.awk window=73 buffer=10000 > $3.pos
