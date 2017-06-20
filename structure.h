#include <iostream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <ctype.h>
#include <zlib.h>
#include <pthread.h>
#include <sys/time.h>

#define ReadChunkSize 1000

#define KmerSize 8
#define KmerPower 0x3FFF

using namespace std;

typedef unsigned char ubyte_t;
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef unsigned short uint16_t;
typedef unsigned long long uint64_t;
typedef unsigned long long bwtint_t;

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

typedef struct
{
	bwtint_t x[3];
} bwtintv_t;

typedef struct
{
	int len;
	int freq;
	bwtint_t* LocArr;
} bwtSearchResult_t;

typedef struct
{
	uint32_t wid; // word id
	uint32_t pos; // occurrence position
} KmerItem_t;

typedef struct
{
	int PosDiff;
	uint32_t rPos;
	uint32_t gPos;
} KmerPair_t;

typedef struct
{
	char* name; // chromosome name
	int64_t FowardLocation;
	int64_t ReverseLocation;
	int64_t len; // chromosome length
} Chromosome_t;

typedef struct
{
	bool bSimple;
	int rPos; // read position
	int64_t gPos; // genome position
	int rLen; // read block size
	int gLen; // genome block size
	int64_t PosDiff; // gPos-rPos
} SeedPair_t;

typedef struct
{
	bool bDir; // true:forward, false:reverse
	string CIGAR;
	int64_t gPos;
	int ChromosomeIdx;
} Coordinate_t;

typedef struct
{
	int Score; // seed score
	int64_t PosDiff;
	int PairedAlnCanIdx;
	vector<SeedPair_t> SeedVec;
} AlignmentCandidate_t;

typedef struct
{
	int AlnScore;
	int SamFlag; // sam flag
	int PairedAlnCanIdx;
	Coordinate_t coor;
} AlignmentReport_t;

typedef struct
{
	int rlen;
	char* header;
	char* seq;
	uint8_t* EncodeSeq;
	// aln report
	int mapq;
	int score;
	int sub_score;
	int CanNum;
	int iBestAlnCanIdx;
	AlignmentReport_t* AlnReportArr;
} ReadItem_t;

// Global variables
extern bwt_t *Refbwt;
extern bwaidx_t *RefIdx;
extern unsigned char nst_nt4_table[256];
extern int64_t GenomeSize, TwoGenomeSize;
extern vector<Chromosome_t> ChromosomeVec;
extern vector<string> ReadVec, ReadHeaderVec;
extern vector<int64_t> AccumulationLengthVec, PositionShiftPosVec;

extern const char* VersionStr;
extern map<int64_t, int> ChrLocMap;
extern bool bDebugMode, bPairEnd, bPacBioData, gzCompressed, FastQFormat, bMultiHit;
extern char *RefSequence, *GenomeFileName, *IndexFileName, *ReadFileName, *ReadFileName2, *OutputFileName;
extern int MaxInsertSize, iThreadNum, iChromsomeNum, WholeChromosomeNum, ChromosomeNumMinusOne, MaxGaps, MinSeedLength, OutputFileFormat;

// bwt_index.cpp
extern void RestoreReferenceInfo();
extern void bwa_idx_destroy(bwaidx_t *idx);
extern bwaidx_t *bwa_idx_load(const char *hint);

// bwt_search.cpp
extern bwtSearchResult_t BWT_Search(uint8_t* seq, int start, int stop);

// Mapping.cpp
extern void Mapping();
extern bool CheckPairedAlignmentCandidates(int64_t EstiDistance, vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2);

// AlignmentCandidates.cpp
extern void RemoveShortSeeds(vector<SeedPair_t>& SeedVec, int thr);
extern bool CompByPosDiff(const SeedPair_t& p1, const SeedPair_t& p2);
extern bool CompByGenomePos(const SeedPair_t& p1, const SeedPair_t& p2);
extern void IdentifyNormalPairs(int rlen, int glen, vector<SeedPair_t>& SeedVec);
extern vector<SeedPair_t> IdentifySeedPairs_FastMode(int rlen, uint8_t* EncodeSeq);
extern vector<SeedPair_t> IdentifySeedPairs_RescueMode(int rlen, uint8_t* EncodeSeq);
extern vector<SeedPair_t> IdentifySeedPairs_SensitiveMode(int rlen, uint8_t* EncodeSeq);
extern void GenMappingReport(bool bFirstRead, ReadItem_t& read, vector<AlignmentCandidate_t>& AlignmentVec);
extern vector<AlignmentCandidate_t> GenerateAlignmentCandidateForPacBioSeq(int rlen, vector<SeedPair_t>& SeedPairVec);
extern vector<AlignmentCandidate_t> GenerateAlignmentCandidateForIlluminaSeq(int rlen, vector<SeedPair_t>& SeedPairVec);
extern Coordinate_t GenCoordinateInfo(bool bFirstRead, int64_t gPos, int64_t end_gPos, vector<pair<int, char> >& cigar_vec);

// AlignmentRescue.cpp
extern bool RescueUnpairedAlignment(int EstDistance, ReadItem_t& r1, ReadItem_t& r2, vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2);

// GetData.cpp
extern bool CheckBWAIndexFiles(string IndexPrefix);
extern int GetNextChunk(bool bSepLibrary, FILE *file, FILE *file2, ReadItem_t* ReadArr);
extern int gzGetNextChunk(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr);

// tools.cpp
extern void ShowSeedLocationInfo(int64_t MyPos);
extern void SelfComplementarySeq(int len, char* rseq);
extern void ShowSeedInfo(vector<SeedPair_t>& SeedPairVec);
extern void GetComplementarySeq(int len, char* seq, char* rseq);
extern int CalFragPairIdenticalBases(int len, char* frag1, char* frag2);
extern int AddNewCigarElements(string& str1, string& str2, vector<pair<int, char> >& cigar_vec);
extern int ProcessHeadSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec);
extern int ProcessTailSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec);
extern int ProcessNormalSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec);
extern void ShowAlignmentCandidateInfo(bool bFirst, char* header, vector<AlignmentCandidate_t>& AlignmentVec);
//extern int ProcessSimpleSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec);

// KmerAnalysis.cpp
extern vector<KmerItem_t> CreateKmerVecFromReadSeq(int len, char* seq);
extern vector<KmerPair_t> IdentifyCommonKmers(int MaxShift, vector<KmerItem_t>& vec1, vector<KmerItem_t>& vec2);
extern vector<SeedPair_t> GenerateSimplePairsFromCommonKmers(int MinSeedLength, vector<KmerPair_t>& KmerPairVec);
extern vector<SeedPair_t> GenerateSimplePairsFromFragmentPair(int MaxDist, int len1, char* frag1, int len2, char* frag2);

// PairwiseAlignment.cpp
extern void PairwiseSequenceAlignment(int m, string& s1, int n, string& s2);
