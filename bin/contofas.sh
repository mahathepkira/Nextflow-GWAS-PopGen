#!/bin/bash

awk '{ printf ">%s\n%s\n",$1,$2 }' $1
