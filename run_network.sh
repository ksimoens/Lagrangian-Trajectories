export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

skip_headers=1
i=0
while IFS=, read -r col1 col2 col3 col4 col5
do
    if ((skip_headers))
    then
        ((skip_headers--))
    else
        echo "I got:$col1|$col2|$col3"
        sed -i "1s/.*/x0 $col2/" config
        sed -i "2s/.*/y0 $col3/" config
        sed -i "4s/.*/write \/media\/kobe\/Windows\/spectrum\/network\/sim_network_test_$(printf "%04d" $col1)/" config
        ((i++))
        exit
        OMP_NUM_THREADS=4 ./main
    fi
done < ../healpix/test_grid.csv

#OMP_NUM_THREADS=4 /usr/bin/time ./simulation 