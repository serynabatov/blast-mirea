#include <stdint.h>
#include <cstring>
#include "stdio.h"
#include <stdlib.h>
#include "blastn.h"

// directions codes
static const int CENTER = 0;
static const int NORTH = 1;
static const int NORTH_WEST = 2;
static const int WEST = 3;

// scores used for Smith Waterman similarity computation
static const short GAP_i = -1;
static const short GAP_d = -1;
static const short MATCH = 2;
static const short MISS_MATCH = -1;

uint32_t lengthSubQuery = 0;
SequencePositionMap querySeq[N];
SequencePositionMap dbSeq[N];
MatchPositionMap finalMatches[N];
uint32_t sequenceQueryLen;
uint32_t sequenceDbLen;
uint32_t mapLen;

// TODO: indexing the database (BLAST 6)
// TODO: fuck esort read again about recursion implementation in FPGA

unsigned int RSHash(char* str, uint32_t len) {
    unsigned int b = 378551;
    unsigned int a = 63689;
    unsigned int hash = 0;
    uint32_t i = 0;

    for (i = 0; i < len; i++){
#pragma HLS pipeline II=3
        hash = hash * a + str[i];
        a = a * b;
    }

    return (hash & 0x7FFFFFFF);
}

// TODO: think about easing the loop
bool search(char* key, uint32_t len, bool mode) {

    unsigned int hashIndex = RSHash(key, len);
	for (uint32_t i = 0; i < len; i++) {
		if (!mode) {
			if (hashIndex == RSHash(querySeq[i].string, len)) {
				return true;
			}
		} else {
			if (hashIndex == RSHash(dbSeq[i].string, len)) {
				return true;
			}
		}
	}
	return false;
}

void pushArrayElement(char* key, uint32_t len, uint32_t posLocation) {

	unsigned int hashIndex = RSHash(key, len);

	for (uint32_t i = 0; i < len; i++) {
		if (RSHash(querySeq[i].string, len) == hashIndex) {
#pragma HLS pipeline II=3
			uint32_t k = querySeq[i].positions.count;
			querySeq[i].positions.positions[k] = posLocation;
			k++;
			querySeq[i].positions.count = k;
			break;
		}
	}

}

void pushArrayElementDb(char* key, uint32_t len, uint32_t posLocation) {
	unsigned int hashIndex = RSHash(key, len);

	for (uint32_t i = 0; i < len; i++) {
		if (RSHash(dbSeq[i].string, len) == hashIndex) {
#pragma HLS pipeline II=3
			uint32_t k = dbSeq[i].positions.count;
			dbSeq[i].positions.positions[k] = posLocation;
			k++;
			dbSeq[i].positions.count = k;
			break;
		}
	}
}

uint16_t checkLetter(char a) {
	uint16_t ans;
	switch(a) {
	case 'A':
		ans = 0;
		break;
	case 'T':
		ans = 1;
		break;
	case 'G':
		ans = 2;
	case 'C':
		ans = 3;
	}
	return ans;
}


void sort() {
	SequencePositionMap temp[N];

	for (int width = 1; width < sequenceDbLen; width = 2 * width) {
		int f1 = 0;
		int f2 = width;
		int i2 = width;
		int i3 = 2 * width;
		if (i2 >= sequenceDbLen) i2 = sequenceDbLen;
		if (i3 >= sequenceDbLen) i3 = sequenceDbLen;

		for (int i = 0; i < sequenceDbLen; i ++) {
#pragma HLS pipeline II=1
			SequencePositionMap t1 = dbSeq[f1];
			SequencePositionMap t2 = dbSeq[f2];
			if ((f1 < i2 && t1.key <= t2.key) || f2 == i3) {
				temp[i] = t1;
				f1++;
			} else {
				assert(f2 < i3);
				temp[i] = t2;
				f2++;
			}
			if (f1 == i2 && f2 == i3) {
				f1 = i3;
				i2 += 2 * width;
				i3 += 2 * width;
				if (i2 >= sequenceDbLen) i2 = sequenceDbLen;
				if (i3 >= sequenceDbLen) i3 = sequenceDbLen;
				f2 = i2;

			}
		}

		// copy temp to dbSeq
		for (int i = 0; i < sequenceDbLen; i++) {
#pragma HLS pipeline II=1
			dbSeq[i] = temp[i];
		}
	}
}

void parseSequence(char* query, uint32_t len, bool mode) {
	char querySequences[N][N + 1];

	for (uint32_t i = 0; i < len; i++) {
	    char p[N + 1];
	    int c;

	    for (c = 0; c < lengthSubQuery; c++) {
#pragma HLS pipeline II=3
	      p[c] = query[i + c];
	    }

	    p[c] = '\0';

	    memcpy(querySequences[i], p, sizeof(char) * (lengthSubQuery + 1));
	}

	uint64_t k = 0;

	for (uint32_t i = 0; i < len; i++) {
		if (!search(querySequences[i], len, mode)) {
			if (!mode) {
				uint32_t powQuery = 1;
				uint32_t key = 0;
				memcpy(querySeq[k].string, querySequences[i],
						sizeof(char) * (lengthSubQuery + 1));
				querySeq[k].positions.positions[0] = i;
				querySeq[k].positions.count = 1;
				for (uint32_t j = 0; j < lengthSubQuery; j++) {
					key += powQuery * (uint32_t)(checkLetter(querySeq[k].string[j]));
					powQuery *= lengthSubQuery;
				}
				querySeq[k].key = key;
				k++;
			} else {
				uint32_t powQuery = 1;
				uint32_t key = 0;
				memcpy(dbSeq[k].string, querySequences[i],
						sizeof(char) * (lengthSubQuery + 1));
				dbSeq[k].positions.positions[0] = 1;
				dbSeq[k].positions.count = 1;
				for (uint32_t j = 0; j < lengthSubQuery; j++) {
					key += powQuery * (uint32_t)(checkLetter(dbSeq[k].string[j]));
					powQuery *= lengthSubQuery;
				}
				dbSeq[k].key = key;
				k++;
			}
		} else {
			if (!mode) {
				pushArrayElement(querySequences[i], len, i);
			} else {
				pushArrayElementDb(querySequences[i], len, i);
			}
		}
	}

	if (!mode) {
		sequenceQueryLen = k;
	} else {
		sequenceDbLen = k;
	}
}

int binarySearch(int l, int r, uint32_t x) {
	while (l <= r) {
		int mid = l + (r - l) / 2;

		if (dbSeq[mid].key == x)
			return mid;

		if (dbSeq[mid].key > x)
			r = mid - 1;
		else
			l = mid + 1;
	}

	// if there is no element
	return -1;
}

void match(uint32_t hppsScore, uint32_t len) {

	//char set[] = {'A', 'T', 'G', 'C'};
	uint32_t k = 0;

	// the original
	for (uint32_t i = 0; i < sequenceQueryLen; i++) {
#pragma HLS unroll factor=2
		int middle = binarySearch(0, sequenceDbLen - 1, querySeq[i].key);
		if ( middle != -1) {
			memcpy(finalMatches[k].string, querySeq[i].string, sizeof(querySeq[i].string));
			finalMatches[k].key = querySeq[i].key;
			memcpy(finalMatches[k].dbPosition.positions, dbSeq[middle].positions.positions,
					sizeof(dbSeq[middle].positions.positions));
			memcpy(finalMatches[k].queryPosition.positions, querySeq[i].positions.positions,
					sizeof(querySeq[i].positions.positions));
			finalMatches[k].dbPosition.count = dbSeq[middle].positions.count;
			finalMatches[k].queryPosition.count = querySeq[i].positions.count;
			k ++;
		}
	}
	mapLen = k;
}

// Smith Waterman Score
int prepareForScore(char* db, char* query, int* similarityMatrix, short* directionMatrix,
		uint32_t dbLength, uint32_t queryLength) {

	int index = 0;
	int i = 0;
	int j = 0;
	short dir = CENTER;
	short match = 0;
	int val = 0;
	int north = 0;
	int west = 0;
	int northwest = 0;
	int testVal = 0;
	int matrixSize = dbLength * queryLength;
	int maxValue[N];

    memset(maxValue, 0, sizeof(int) * N);

	int k = 0;

	for (int d = 0; d < mapLen; d++) {
		for (int e = 0; e < finalMatches[d].dbPosition.count; e ++) {
			for (int f = 0; f < finalMatches[d].queryPosition.count; f ++) {
				int max = 0;

				i = finalMatches[d].dbPosition.positions[e] % dbLength; // column index
				j = finalMatches[d].queryPosition.positions[f] / queryLength; //row index

				val = 0;
				west = 0;
				northwest = 0;
				testVal = 0;

				for (index = dbLength * i + j; index < matrixSize; index++) {
					dir = CENTER;
					val = 0;

					if (i == 0) {
						// first column
						west = 0;
						northwest = 0;
					} else {
						//printf("\n %d \n", northwest);

						north = similarityMatrix[index - dbLength];
						//printf("%d\n", similarityMatrix[15]);
						//printf("%d\n", similarityMatrix[uint32_t(index - 12)]);
						//printf("%d %d \n", north, index);
						match = (db[i] == query[j]) ? MATCH : MISS_MATCH;

//						if (index == 15)
//								printf("%d %d\n", northwest, match);

						testVal = northwest + (int)match;
//						Prints("test val - %d\n", testVal);
//
//						if (testVal < 0)
//							continue;

						if (testVal > val) {
							val = testVal;
							dir = NORTH;
						}

						testVal = north + (int)GAP_d;
//						printf("test val - %d\n", testVal);

////
//						if (testVal < 0)
//							continue;

						if (testVal > val) {
							val = testVal;
							dir = NORTH;
						}

						testVal = west + (int)GAP_i;
//						printf("test val - %d\n", testVal);

//						if (testVal < 0)
//							continue;


						if (testVal > val) {
							val = testVal;
							dir = WEST;
						}

//						printf("index/val %d %d\n", index, val);
						similarityMatrix[index] = val;
//						printf("index/sim %d %d\n", index, similarityMatrix[index]);
						directionMatrix[index] = dir;
						west = val;
						northwest = north;
//						printf("west - %d\n", west);

						if (val > max) {
							//maxIndex[0] = index;
							max = val;
						//	k++;
						}
					}

					i = index % dbLength;
					j = index / queryLength;
				}
				maxValue[k] = max;
				k++;
			}
		}
	}

	return index;
}

void top(char* query, char* db, uint32_t subQueryLength, uint8_t cmd,
		 uint32_t queryLength, uint32_t dbLength, uint32_t hppsScore,
		 int* similarityMatrix, short* directionMatrix) {
#pragma HLS INTERFACE s_axilite port=return

#pragma HLS INTERFACE s_axilite offset=slave port=query
#pragma HLS INTERFACE s_axilite offset=slave port=db
#pragma HLS INTERFACE s_axilite offset=slave port=subQueryLength
#pragma HLS INTERFACE s_axilite offset=slave port=cmd
#pragma HLS INTERFACE s_axilite offset=slave port=queryLength
#pragma HLS INTERFACE s_axilite offset=slave port=dbLength
#pragma HLS INTERFACE s_axilite offset=slave port=hppsScore
#pragma HLS INTERFACE s_axilite offset=slave port=similarityMatrix
#pragma HLS INTERFACE s_axilite offset=slave port=directionMatrix

#pragma HLS INTERFACE m_axi port=query bundle=gmem
#pragma HLS INTERFACE m_axi port=db bundle=gmem
#pragma HLS INTERFACE m_axi port=similarityMatrix bundle=gmem
#pragma HLS INTERFACE m_axi port=directionMatrix bundle=gmem

	switch(cmd) {
		// initialize arrays within the memory
		case 1: {
			lengthSubQuery = subQueryLength;
			uint32_t maxLengthQuery = queryLength - lengthSubQuery + 1;
			parseSequence(query, maxLengthQuery, false);
			uint32_t maxLengthDb = dbLength - lengthSubQuery + 1;
			parseSequence(db, maxLengthDb, true);
			match(hppsScore, subQueryLength);

			int ans = prepareForScore(db, query, similarityMatrix, directionMatrix, dbLength, queryLength);
			//printf("%d\n", ans);
			break;
		}
		// sort by keys
		case 2: {
//			sort();
//			for (uint32_t i = 0; i < mapLen; i++) {
//				printf("%s %d ", dbSeq[i].string, dbSeq[i].key);
//				for (uint32_t j = 0; j < dbSeq[i].positions.count; j++) {
//					printf("%d; ", dbSeq[i].positions.positions[j]);
//				}
//				printf("\n");
//			}
			break;
		}
		// search for the matches
		case 3: {
			match(hppsScore, subQueryLength);
			break;
		}
		case 4: {
			prepareForScore(db, query, similarityMatrix, directionMatrix, dbLength, queryLength);
			break;
		}
	}
}
