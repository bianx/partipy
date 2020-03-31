# Runs

    ./main -m 1000 -s rk4 -t 40 -d 0.02 -e 10 -o punto -c gauss  < initial/circle > p
    ./main -m 1000 -s rk4 -t 40 -d 0.02 -e 10 -o punto -c chorin < initial/circle > p
    
    punto -D 2 -c -G p
