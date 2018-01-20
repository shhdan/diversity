/*
* zorder.cpp
*
*  Created on: 12/07/2016
*      Author: uqjgan
*/


/*
*  Compare the Z-order directly by the original coordinates.
*/
/* added by SHD*/
//#include "stdafx.h"

#include "zorder.h"

bool compZorderByCoord(const int* v1, const int* v2, int _dim) {
	int pos = 0;
	unsigned int curMax = 0;
	unsigned int temp = 0;
	for (int i = 0; i < _dim; i++) {
		temp = v1[i] ^ v2[i];
		if (curMax < temp && curMax < (temp ^ curMax)) {
			// Record the dimension where the different most significant bit lies.
			pos = i;
			curMax = temp;
		}
	}
	return v1[pos] < v2[pos];
}

