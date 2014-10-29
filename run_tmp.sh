#!/bin/bash

ifort -o Module_Numerical_Integration.o -c Module_Numerical_Integration.f90
ifort -o Module_Complementaire.o -c Module_Complementaire.f90
ifort -o Module_Conservation.o -c Module_Conservation.f90
ifort -o Module_Init_tmp.o -c Module_Init_tmp.f90
ifort -o Module_Output.o -c Module_Output.f90
ifort -o Module_Surface.o -c Module_Surface.f90
ifort -o Module_Thermal_GFD.o -c Module_Thermal_GFD.f90
ifort -o Module_Thermal_Newton.o -c Module_Thermal_Newton.f90
ifort -o Module_Thermal_Newton_OLD.o -c Module_Thermal_Newton_OLD.f90
ifort -o Module_Thermal.o -c Module_Thermal.f90
ifort -o Module_Thickness_GFD.o -c Module_Thickness_GFD.f90
ifort -o Module_Thickness_Newton.o -c Module_Thickness_Newton.f90
ifort -o Module_Thickness.o -c Module_Thickness.f90
ifort -o main.o -c main.f90
ifort -o E1D0_G0D0_N1D-3_P1D-3_D1D-2_C-3D-1_R1D0_S2D-2_Dr1D-2_Ep1D-4_Dt1D-6  main.o Module_Complementaire.o Module_Conservation.o Module_Init_tmp.o Module_Numerical_Integration.o Module_Output.o Module_Surface.o Module_Thermal.o Module_Thermal_GFD.o Module_Thermal_Newton.o Module_Thermal_Newton_OLD.o Module_Thickness.o Module_Thickness_GFD.o Module_Thickness_Newton.o

rm *.o
rm *.mod

