#!/bin/bash

g++ -lpng -std=c++11 FourierTransform.cxx ConstantQTransform.cxx ../util/png_image.cc cqtgen.cc -o cqtgen
