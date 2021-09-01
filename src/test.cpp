#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "top.h"
#include <iostream>

#define N 12

void top(char* query, char* db, uint32_t subQueryLength, uint8_t cmd,
		uint32_t queryLength, uint32_t dbLength, uint32_t hppsScore, int* similarityMatrix,
		short* directionMatrix);

int main(int argc, char** argv) {

	uint32_t subqueryLength = 3;
    char* query = "AAAATGGTTCCC";
    char* db = "AAAAAAAAAATT";
    printf("db and query are defined!\n");
    fflush(stdout);

    int* similarityMatrix = (int *) malloc (sizeof(int) * N * N);
    short* directionMatrix = (short *) malloc(sizeof(short) * N * N);
    memset(similarityMatrix, 0, sizeof(int) * N * N);
    memset(directionMatrix, 0, sizeof(short) * N * N);

    uint8_t cmd = 1;

	int index, i;

	top(query, db, subqueryLength, cmd, 12, 12, 1, similarityMatrix, directionMatrix);

	//cmd = 2;
	//top(query, db, subqueryLength, cmd, 12, 12);

	for (index = 0; index < N * N; index ++) {

		//printf("%d\t", similarityMatrix[index]);
		std::cout << similarityMatrix[index] << "\t";
		if (index % N == 11)
			//printf("\n");
			std::cout << std::endl;
	}

	return 0;
}
