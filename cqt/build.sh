#!/bin/bash

g++ -I cq-v1.1/cq -lpng -L. -lcq -std=c++11 FourierTransform.cxx ConstantQTransform.cxx test.cpp
