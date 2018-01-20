/*
* zorder.h
*
*  Created on: 12/07/2016
*      Author: uqjgan
*/

#ifndef ZORDER_H_
#define ZORDER_H_

/*
*  Compare the Z-order directly by the original coordinates.
*  Return:
*  	If the Z-order of v1 is <= that of v2, return TRUE. Otherwise, return FALSE.
*/
bool compZorderByCoord(const int* v1, const int* v2, int _dim);


#endif /* ZORDER_H_ */

