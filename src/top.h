#ifndef _HIER_FUNC_H
#define _HIER_FUNC_H

#include "stdint.h"

extern void top(char* query, char* db, uint32_t subQueryLength, uint32_t cmd, uint32_t queryLength, uint32_t dbLength, uint32_t hppsScore,
		int* similarityMatrix, short* directionMatrix);

#endif
