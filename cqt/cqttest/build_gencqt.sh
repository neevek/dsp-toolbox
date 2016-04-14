#!/bin/bash

g++ -lpng -std=c++11 FourierTransform.cxx ConstantQTransform.cxx png_image.cc gencqt.cc -o cqtgen
