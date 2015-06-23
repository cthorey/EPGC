#!/bin/bash
cd "$(dirname "$0")"
Compilateur -o Mobility_Thermal_Skin_Rheology.o -c Mobility_Thermal_Skin_Rheology.f90
Compilateur -o Mobility_Thickness_Skin_Rheology.o -c Mobility_Thickness_Skin_Rheology.f90
Compilateur -o Module_Numerical_Integration.o -c Module_Numerical_Integration.f90
Compilateur -o Module_Complementaire.o -c Module_Complementaire.f90
Compilateur -o Module_Conservation.o -c Module_Conservation.f90
Compilateur -o Module_Init_tmp.o -c Module_Init_tmp.f90
Compilateur -o Module_Mobility.o -c Module_Mobility.f90
Compilateur -o Module_Output.o -c Module_Output.f90
Compilateur -o Module_Surface.o -c Module_Surface.f90
Compilateur -o Module_Thermal_GFD_GRAVI_Bercovici.o -c Module_Thermal_GFD_GRAVI_Bercovici.f90
Compilateur -o Module_Thermal_IntE_Newton_Bercovici.o -c Module_Thermal_IntE_Newton_Bercovici.f90
Compilateur -o Module_Thermal_Skin_GFD_Bercovici.o -c Module_Thermal_Skin_GFD_Bercovici.f90
Compilateur -o Module_Thermal_Skin_GFD_GRAVI.o -c Module_Thermal_Skin_GFD_GRAVI.f90
Compilateur -o Module_Thermal_Skin_Newton.o -c Module_Thermal_Skin_Newton.f90
Compilateur -o Module_Thermal_Skin_Newton_Arrhenius.o -c Module_Thermal_Skin_Newton_Arrhenius.f90
Compilateur -o Module_Thermal_Skin_Newton_Bercovici.o -c Module_Thermal_Skin_Newton_Bercovici.f90
Compilateur -o Module_Thermal_Skin_Newton_Roscoe.o -c Module_Thermal_Skin_Newton_Roscoe.f90
Compilateur -o Module_Thermal.o -c Module_Thermal.f90
Compilateur -o Module_Thickness_Inte_GFD_Bercovici.o -c Module_Thickness_Inte_GFD_Bercovici.f90
Compilateur -o Module_Thickness_Skin_Newton.o -c Module_Thickness_Skin_Newton.f90
Compilateur -o Module_Thickness_Skin_Newton_Arrhenius.o -c Module_Thickness_Skin_Newton_Arrhenius.f90
Compilateur -o Module_Thickness_Skin_Newton_Bercovici.o -c Module_Thickness_Skin_Newton_Bercovici.f90
Compilateur -o Module_Thickness_Newton_Gravi.o -c Module_Thickness_Newton_Gravi.f90
Compilateur -o Module_Thickness.o -c Module_Thickness.f90
Compilateur -o Module_Thickness_Skin_Newton_Roscoe.o -c Module_Thickness_Skin_Newton_Roscoe.f90
Compilateur -o lib_array.o -c lib_array.f90
Compilateur -o main.o -c main.f90
Compilateur -o run lib_array.o main.o Mobility_Thermal_Skin_Rheology.o Mobility_Thickness_Skin_Rheology.o Module_Complementaire.o Module_Conservation.o Module_Init_tmp.o Module_Mobility.o Module_Numerical_Integration.o Module_Output.o Module_Surface.o Module_Thermal.o Module_Thermal_GFD_GRAVI_Bercovici.o Module_Thermal_IntE_Newton_Bercovici.o Module_Thermal_Skin_GFD_Bercovici.o Module_Thermal_Skin_GFD_GRAVI.o Module_Thermal_Skin_Newton.o Module_Thermal_Skin_Newton_Arrhenius.o Module_Thermal_Skin_Newton_Bercovici.o Module_Thermal_Skin_Newton_Roscoe.o Module_Thickness.o Module_Thickness_Inte_GFD_Bercovici.o Module_Thickness_Newton_Gravi.o Module_Thickness_Skin_Newton.o Module_Thickness_Skin_Newton_Arrhenius.o Module_Thickness_Skin_Newton_Bercovici.o Module_Thickness_Skin_Newton_Roscoe.o

rm *.o
rm *.mod

