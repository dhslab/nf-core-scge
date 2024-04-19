#!/usr/bin/env bash

PWD=$(realpath .)

for i in *.orig;
do
    sed -e 's!REPLACEPATH!'$PWD'!g' $i > $(basename $i .orig)
done


