#!/bin/bash
mpicc fsecuential.c -o fsecuential
mpicc fcannon.c -o fcannon 
mpicc isecuential.c -o isecuential
mpicc icannon.c -o icannon
for i in $(seq 1 1 100)
do
   mpirun -np 4 ./fsecuential>>ftiemposs.txt
done
echo "Float Secuencial listo"
for i in $(seq 1 1 100)
do
   mpirun -np 4 ./fcannon >>ftiemposc.txt
done
echo "Float Cannon listo"
for i in $(seq 1 1 100)
do
   mpirun -np 4 ./isecuential >>itiemposs.txt
done
echo "Int Secuencial listo"
for i in $(seq 1 1 100)
do
   mpirun -np 4 ./icannon >>itiemposc.txt
done
echo "Int Cannon listo"
echo "Todo listo"
