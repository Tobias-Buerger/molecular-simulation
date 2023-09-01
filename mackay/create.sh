#! /bin/bash

for d in {2..20} ; do
   ./mackay $d 3 > cluster_$d.xyz
   x=$(head -n 1 cluster_$d.xyz)
   mv cluster_$d.xyz cluster_$x.xyz
done
