#! /bin/bash

awk     '{if (($1 == "final") && ($2 == "lattice") ) {getline;getline;getline;print $0                  }}' $1 >> a.dat
awk     '{if (($1 == "final") && ($2 == "lattice") ) {getline;getline;getline;getline; print $0         }}' $1 >> b.dat
awk     '{if (($1 == "final") && ($2 == "lattice") ) {getline;getline;getline;getline;getline; print $0 }}' $1 >> c.dat

