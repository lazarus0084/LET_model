#!/bin/bash

INPUTFILE=input1.txt

rm *.o iitpave2.out
make iitpave2
./iitpave2.out $INPUTFILE