#!/bin/bash

INPUTFILE=input2.txt

rm *.o iitpave2.out
make iitpave2
./iitpave2.out $INPUTFILE