// BLASTN
#include<stdint.h>
#include<assert.h>

#define N 1024
#define M 256

#define A 0
#define T 1
#define G 2
#define C 3

typedef struct Position {
	uint32_t positions[N]; // [1, 3] - start positions
	uint32_t count; // number of position
} Position;

typedef struct SequencePositionMap {
	char string[N];   // TODO: delete it, it's for debug
	Position positions;
	uint32_t key; // k-mers value here
} SequencePositionMap;

typedef struct MatchPositionMap {
	char string[N];
	uint32_t key;
	Position dbPosition;
	Position queryPosition;
} MatchPositionMap;

extern void parseSequence (char* query, uint32_t len, bool mode);

// the subquery that our query will be divided
extern uint32_t lengthSubQuery;

extern SequencePositionMap querySeq[N];
extern uint32_t sequenceQueryLen;
extern SequencePositionMap dbSeq[N];
extern uint32_t sequenceDbLen;
extern void match(uint32_t hppsScore, uint32_t len);

extern void sort();

//	static uint32_t* extend(MatchPositionMap* matches,
//						char const* query,
//						char const* database,
//						int subqueryLength,
//						uint32_t length);
//
//	static float similarity(uint32_t* maxValues, uint32_t len);
//
//
extern int prepareForScore(char* db, char* query, int* similarityMatrix, short* directionMatrix, uint32_t dbLength, uint32_t queryLength);

extern char* substring(char* string, int position);

extern unsigned int RSHash(char const* str, uint32_t len);

extern bool search(char* key, uint32_t len, bool mode);

extern void pushArrayElement(char* key, uint32_t len, uint32_t posLocation);
extern void pushArrayElementDb(char* key, uint32_t len, uint32_t posLocation);
//	static void print(SequencePositionMap* head);
//
//	static void printMatch(MatchPositionMap* head);
//
//	static bool searchMatch(MatchPositionMap* seq, char const* key, uint32_t len);
//
//	static void pushMatch(MatchPositionMap** seq, char* key, uint32_t* posQueryLocation, uint32_t* posDbLocation,
//						 uint32_t dbCount, uint32_t queryCount);
//};
