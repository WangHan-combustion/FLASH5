The table contains the following information in the following format

Row 2: total number of possible variables	number of derivatives	number of variables actually used in the table 
Row 3: number of cells computed in each direction
Row 4: Flags for whether log is used in each direction
Row 5: energy array size for this table		total energy array size		energy array giving cell boundary values
Row 6: density array size for this table	total density array size	density array giving cell boundary values
Row 7: Ye array size for this table		total Ye array size		Ye array giving cell boundary values
Row 8: abar array size for this table		total abar array size 		abar array giving cell boundary values

Row 9-10: degree of derivative for energy, density, Ye, and abar.
Example: 
1 0 0 0
0 1 0 0
meaning the first derivative w.r.t energy is given first in the gaussian weight sets below and the first derivative w.r.t density, 0100, is given second in the gaussian weight sets below.

Rows 11-r:
Columns 0-c: the weights for entropy and two first derivatives 
Columns c+1-c+4: the four hyperparameters  L1, L2, L3, L4, for energy, density, Ye, and abar respectively

Note: Rows 11-r are created using four(4) for loops. The ordering of the for loops are: abar, Ye, density, and energy so abar varies the slowest while energy varies the fastest 
