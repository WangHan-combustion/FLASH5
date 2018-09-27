The table contains the following information in the following format

Row 2: number of variables	number of derivatives
Row 3: number of cells computed in each direction
Row 4: Flags for whether log is used in each direction
Row 5: energy array size for this table 	energy array giving cell boundaries
Row 6: density array size for this table	density array giving cell boundaries
Row 7: Ye array size for this table		Ye array giving cell boundaries
Row 8: abar array size for this table		abar array giving cell boundaries
Row 9: degree of derivative for energy, density, Ye, and abar. Group this row of number in sets of four, Example 1 0 0 0 0 1 0 0 means 1000 and 0100 meaning the first derivative w.r.t energy is given first in the gaussian weight sets below and the first derivative w.r.t density, 0100, is given second in the gaussian weight sets below.

Rows 10-25:
Columns 0-48: the weights for entropy and two first derivatives 
Columns 49-52: the four hyperparameters  L1, L2, L3, L4, for energy, density, Ye, and sumy respectively

Note: Rows 10-25 are created using four(4) for loops. The ordering of the for loops follows the ordering of Rows 5-8: energy, density, Ye, and then abase 
