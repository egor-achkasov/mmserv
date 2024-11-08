#!/bin/fish

date
module load questa/10v7a

set root $PWD
set defconf config/default.mk
set result "$root/results.csv"

rm -f $defconf

echo 'NR_LANES;VLEN;datalayout;order;M;K;cycles' > $result

for nrl in 2 4 8 16
    for vl in 512 1024 2048 4096
        # Write config
        echo "nr_lanes ?= $nrl" > $defconf
        echo "vlen ?= $vl" >> $defconf
        # Compile binary
        cd apps
        rm -f gfdmrv/src/gfdrmv.c.o
        rm -f gfdmrv/data.S.o
        make bin/gfdmrv
        make bin/gfdmrv
        cd $root
        #Run simulation
        cd hardware
        app=gfdmrv timeout -s KILL 1h make simc
        grep "~~~~~" build/transcript
        if set resline (grep "~~~~~" build/transcript)
            echo $resline | sed "s/# ~~~~~/$nrl;$vl;/" >> $result
        else
            echo "Something went wrong for NR_LANES $nrl and VLEN $vl"
        end
        cd $root
    end
end

date
