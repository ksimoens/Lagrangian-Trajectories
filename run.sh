export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

skip_headers=1
i=0
while IFS=, read -r col1 col2
do
    if ((skip_headers))
    then
        ((skip_headers--))
    else
        echo "I got:$col1|$col2"
        sed -i "1s/.*/x0 $col1/" config
        sed -i "2s/.*/y0 $col2/" config
        sed -i "4s/.*/write sim_sargasso_vel_$(printf "%02d" $i)/" config
        ((i++))
        OMP_NUM_THREADS=4 ./main
    fi
done < initial/sargasso_points_vel.csv

#OMP_NUM_THREADS=4 /usr/bin/time ./simulation 