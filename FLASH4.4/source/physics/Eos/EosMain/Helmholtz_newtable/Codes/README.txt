{\rtf1\ansi\ansicpg1252\cocoartf1561\cocoasubrtf600
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 README for EOS Transpiling Project:\
\
BUILDING THE ENTROPY EOS TABLE:\
Entropic_EOS_BuildTable.py\
	call to EOS_minimum_energy_calculator.py\
	call to Entropic_EOS.py\
		call to Entropic_EOS_SingleFile.py\
			call to FermiDiracIntegralCalculator.py\
		call to TwoDNewtonRaphson.py\
			\
\
READING TABLE:\
EOS_ReadTable_CalcTDState.py\
	call to EOS_MLL_Solver.py\
\
\
Additional Notes:\
\
- Included in the folder is the TimmesEOS_Fortran code that was the basis for EntropicEOS_SingleFile.py\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 - Included in the folder is the fermi_dirac_with_derivatives_Fortran code that was the basis for FermiDiracIntegralCalculator.py}