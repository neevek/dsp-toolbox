#!/bin/bash

g++ -lpng -std=c++11 FourierTransform.cxx ConstantQTransform.cxx test_constantq.cc -o cqtgen
