#!/usr/bin/bash

read -p "voulez vous rebuild le projet ? [Y/n]:" reponse 

cd ../build

if [ "$reponse" = "Y" ];
then
    cmake ..
    make
fi 
 
./myfem

