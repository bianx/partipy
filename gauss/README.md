# Install

    (cd ../lib/barnes-hut && make)
    make

# Runs

    ./main -m 1000 -s rk4 -t 40 -d 0.02 -e 10 -o punto -c gauss  -a n2 < initial/circle > p
    ./main -m 1000 -s rk4 -t 40 -d 0.02 -e 10 -o punto -c chorin -a n2 < initial/circle > p
    
    punto -D 2 -c -G p


Axisymetrization: Koumoutsakos, J. Comput. Phys. 1997

   ./main -m 100 -s rk4 -t 12.0 -d 0.04 -e 1 -o punto -c chorin -a n2 < initial/vorI  > a

   ./main -m 100 -s rk4 -t 1.0 -d 0.04 -e 1 -o punto -c chorin -a n2 < initial/vorI  > a

   ./main -r 20 -m 200 -s rk4 -t 10.0 -d 0.08 -e 1 -o punto -c chorin -a n2 < initial/vorI  > a

punto -D 2 -c -G a


# Store

