#!/bin/bash
cd "$(dirname "$0")"
Compilateur -O3 -fp-model fast=2 -o Module_Numerical_Integration.o -c Module_Numerical_Integration.f90
Compilateur -O3 -fp-model fast=2 -o Module_Complementaire.o -c Module_Complementaire.f90
Compilateur -O3 -fp-model fast=2 -o Module_Conservation.o -c Module_Conservation.f90
Compilateur -O3 -fp-model fast=2 -o Module_Init_tmp.o -c Module_Init_tmp.f90
Compilateur -O3 -fp-model fast=2 -o Module_Output.o -c Module_Output.f90
Compilateur -O3 -fp-model fast=2 -o Module_Surface.o -c Module_Surface.f90
Compilateur -O3 -fp-model fast=2 -o Module_Thermal_GFD.o -c Module_Thermal_GFD.f90
Compilateur -O3 -fp-model fast=2 -o Module_Thermal_Newton.o -c Module_Thermal_Newton.f90
Compilateur -O3 -fp-model fast=2 -o Module_Thermal_Newton_OLD.o -c Module_Thermal_Newton_OLD.f90
Compilateur -O3 -fp-model fast=2 -o Module_Thermal.o -c Module_Thermal.f90
Compilateur -O3 -fp-model fast=2 -o Module_Thickness_GFD.o -c Module_Thickness_GFD.f90
Compilateur -O3 -fp-model fast=2 -o Module_Thickness_Newton.o -c Module_Thickness_Newton.f90
Compilateur -O3 -fp-model fast=2 -o Module_Thickness.o -c Module_Thickness.f90
Compilateur -O3 -fp-model fast=2 -o main.o -c main.f90
Compilateur -O3 -fp-model fast=2 -o run  main.o Module_Complementaire.o Module_Conservation.o Module_Init_tmp.o Module_Numerical_Integration.o Module_Output.o Module_Surface.o Module_Thermal.o Module_Thermal_GFD.o Module_Thermal_Newton.o Module_Thermal_Newton_OLD.o Module_Thickness.o Module_Thickness_GFD.o Module_Thickness_Newton.o

rm *.o
rm *.mod

