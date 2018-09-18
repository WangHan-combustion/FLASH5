The table contains the following information in the following format

Row 1: number of variables	number of derivatives
Row 2: number of cells computed in each direction
Row 3: Flags for whether log is used in each direction
Row 4: energy array size	energy array
Row 5: density array size	density array
Row 6: abar array size	abar array
Row 7: zbar array size	zbar array	!!! Note the zbar vector is created with descending values unlike the other arrays. I will fix this.

Row 8-16:
Columns 0-48: the weights for entropy and two first derivatives 
Columns 49-52: the four hyper parameters  L1, L2, L3, L4, for energy, density, abar, and zbar respectively
