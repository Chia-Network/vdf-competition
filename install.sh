#!/bin/sh
sudo apt-get install libgmp3-dev -y
sudo apt-get install libflint-dev -y
g++ -O3 vdf.cpp -lgmpxx -lgmp -lflint -lmpfr -lpthread
