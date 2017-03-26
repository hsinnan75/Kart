#include <cmath>
#include "structure.h"

#define MAPQ_COEF 30

FILE *output = stdout;
int64_t iDistance = 0;
time_t StartProcessTime;
bool bSepLibrary = false;
fstream ReadFileHandler1, ReadFileHandler2;
static pthread_mutex_t LibraryLock, OutputLock;
int iTotalReadNum = 0, iUnMapped = 0, iPaired = 0;

void SetSingleAlignmentFlag(ReadItem_t& read)
{
	int i;

	if (read.score > read.sub_score || bMultiHit == false) // unique mapping or bMultiHit=false
	{
		i = read.iBestAlnCanIdx;
		if (read.AlnReportArr[i].coor.bDir == false) read.AlnReportArr[i].iFrag = 0x10;
		else read.AlnReportArr[i].iFrag = 0;
	}
	else if(read.score > 0)
	{
		for (i = 0; i < read.CanNum; i++)
		{
			if (read.AlnReportArr[i].AlnScore > 0)
			{
				if (read.AlnReportArr[i].coor.bDir == false) read.AlnReportArr[i].iFrag = 0x10;
				else read.AlnReportArr[i].iFrag = 0;
			}
		}
	}
	else
	{
		read.AlnReportArr[0].iFrag = 0x4;
	}
}

void SetPairedAlignmentFlag(ReadItem_t& read1, ReadItem_t& read2)
{
	int i, j;

	//printf("read1:[%d, %d]:#%d, read2:[%d, %d] #%d\n", read1.score, read1.sub_score, read1.iBestAlnCanIdx, read2.score, read2.sub_score, read2.iBestAlnCanIdx); fflush(stdout);
	if (read1.score > read1.sub_score && read2.score > read2.sub_score) // unique mapping
	{
		i = read1.iBestAlnCanIdx;
		read1.AlnReportArr[i].iFrag = 0x41; // read1 is the first read in a pair

		j = read2.iBestAlnCanIdx;
		read2.AlnReportArr[j].iFrag = 0x81; // read2 is the second read in a pair

		if (j == read1.AlnReportArr[i].PairedAlnCanIdx) // reads are mapped in a proper pair
		{
			read1.AlnReportArr[i].iFrag |= 0x2;
			read2.AlnReportArr[j].iFrag |= 0x2;
		}
		read1.AlnReportArr[i].iFrag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);
		read2.AlnReportArr[j].iFrag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);
	}
	else
	{
		if (read1.score > read1.sub_score || bMultiHit == false) // unique mapping or bMultiHit=false
		{
			i = read1.iBestAlnCanIdx;
			read1.AlnReportArr[i].iFrag = 0x41; // read1 is the first read in a pair
			read1.AlnReportArr[i].iFrag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);
			if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0) read1.AlnReportArr[i].iFrag |= 0x2;// reads are mapped in a proper pair
			else read1.AlnReportArr[i].iFrag |= 0x8; // next segment unmapped
		}
		else if(read1.score > 0)
		{
			for (i = 0; i < read1.CanNum; i++)
			{
				if (read1.AlnReportArr[i].AlnScore > 0)
				{
					read1.AlnReportArr[i].iFrag = 0x41; // read1 is the first read in a pair
					read1.AlnReportArr[i].iFrag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);

					if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0) read1.AlnReportArr[i].iFrag |= 0x2;// reads are mapped in a proper pair
					else read1.AlnReportArr[i].iFrag |= 0x8; // next segment unmapped
				}
			}
		}
		else
		{
			read1.AlnReportArr[0].iFrag = 0x41; // read1 is the first read in a pair
			read1.AlnReportArr[0].iFrag |= 0x4;

			if (read2.score == 0) read1.AlnReportArr[0].iFrag |= 0x8; // next segment unmapped
			else read1.AlnReportArr[0].iFrag |= (read2.AlnReportArr[read2.iBestAlnCanIdx].coor.bDir ? 0x10 : 0x20);
		}

		if (read2.score > read2.sub_score || bMultiHit == false) // unique mapping or bMultiHit=false
		{
			j = read2.iBestAlnCanIdx;
			read2.AlnReportArr[j].iFrag = 0x81; // read2 is the second read in a pair
			read2.AlnReportArr[j].iFrag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);

			if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0) read2.AlnReportArr[j].iFrag |= 0x2;// reads are mapped in a proper pair
			else read2.AlnReportArr[j].iFrag |= 0x8; // next segment unmapped
		}
		else if (read2.score > 0)
		{
			for (j = 0; j < read2.CanNum; j++)
			{
				if (read2.AlnReportArr[j].AlnScore > 0)
				{
					read2.AlnReportArr[j].iFrag = 0x81; // read2 is the second read in a pair
					read2.AlnReportArr[j].iFrag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);

					if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0) read2.AlnReportArr[j].iFrag |= 0x2;// reads are mapped in a proper pair
					else read2.AlnReportArr[j].iFrag |= 0x8; // next segment unmapped
				}
			}
		}
		else
		{
			read2.AlnReportArr[0].iFrag = 0x81; // read2 is the second read in a pair
			read2.AlnReportArr[0].iFrag |= 0x4; // segment unmapped
			if (read1.score == 0) read2.AlnReportArr[0].iFrag |= 0x8; // next segment unmapped
			else read2.AlnReportArr[0].iFrag |= (read1.AlnReportArr[read1.iBestAlnCanIdx].coor.bDir ? 0x10 : 0x20);
		}
	}
}

void EvaluateMAPQ(ReadItem_t& read)
{
	//float f;

	if (read.score == 0 || read.score == read.sub_score) read.mapq = 0;
	else
	{
		if (read.sub_score == 0 || read.score - read.sub_score > 10) read.mapq = 60;
		else
		{
			read.mapq = (int)(MAPQ_COEF * (1 - (float)(read.score - read.sub_score)/read.score)*log(read.score) + 0.4999);
			if (read.mapq > 60) read.mapq = 60;

			//f = 1.0*read.score / read.rlen;
			//read.mapq = (f < 0.95 ? (int)(read.mapq * f * f) : read.mapq);
		}
	}
}

void OutputPairedSamFile(ReadItem_t& read1, ReadItem_t& read2)
{
	char *seq, *rseq;
	int i, j, dist = 0;

	if (read1.score == 0)
	{
		iUnMapped++;
		fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read1.header, read1.AlnReportArr[0].iFrag, read1.seq);
	}
	else
	{
		seq = read1.seq; rseq = NULL;
		for (i = read1.iBestAlnCanIdx; i < read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore > 0)
			{
				if (read1.AlnReportArr[i].coor.bDir == false && rseq == NULL)
				{
					rseq = new char[read1.rlen + 1]; rseq[read1.rlen] = '\0';
					GetComplementarySeq(read1.rlen, seq, rseq);
				}
				if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0)
				{
					dist = (int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen));
					if (i == read1.iBestAlnCanIdx)
					{
						iPaired += 2;
						if (abs(dist) < 10000) iDistance += abs(dist);
					}
					fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t=\t%ld\t%d\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header, read1.AlnReportArr[i].iFrag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), read2.AlnReportArr[j].coor.gPos, dist, (read1.AlnReportArr[i].coor.bDir ? seq : rseq), read1.rlen - read1.score, read1.score, read1.sub_score);
				}
				else fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header, read1.AlnReportArr[i].iFrag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (read1.AlnReportArr[i].coor.bDir ? seq : rseq), read1.rlen - read1.score, read1.score, read1.sub_score);
			}
			if (!bMultiHit) break;
		}
		if (rseq != NULL)
		{
			delete[] rseq;
			rseq = NULL;
		}
	}

	if (read2.score == 0)
	{
		iUnMapped++;
		fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read2.header, read2.AlnReportArr[0].iFrag, read2.seq);
	}
	else
	{
		rseq = read2.seq; seq = NULL;
		for (j = read2.iBestAlnCanIdx; j < read2.CanNum; j++)
		{
			if (read2.AlnReportArr[j].AlnScore > 0)
			{
				if (read2.AlnReportArr[j].coor.bDir == true && seq == NULL)
				{
					seq = new char[read2.rlen + 1]; seq[read2.rlen] = '\0';
					GetComplementarySeq(read2.rlen, rseq, seq);
				}
				if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0)
				{
					dist = 0 - ((int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen)));
					fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t=\t%ld\t%d\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header, read2.AlnReportArr[j].iFrag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), read1.AlnReportArr[i].coor.gPos, dist, (read2.AlnReportArr[j].coor.bDir ? seq : rseq), read2.rlen - read2.score, read2.score, read2.sub_score);
				}
				else fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header, read2.AlnReportArr[j].iFrag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (read2.AlnReportArr[j].coor.bDir ? seq : rseq), read2.rlen - read2.score, read2.score, read2.sub_score);
			}
			if (!bMultiHit) break;
		}
		if (seq != NULL)
		{
			delete[] seq;
			seq = NULL;
		}
	}
}

void OutputSingledSamFile(ReadItem_t& read)
{
	if (read.score == 0)
	{
		iUnMapped++;
		fprintf(output, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t*\tAS:i:0\tXS:i:0\n", read.header, read.AlnReportArr[0].iFrag, read.seq);
	}
	else
	{
		int i;
		char *seq, *rseq;

		seq = read.seq; rseq = NULL;
		for (i = read.iBestAlnCanIdx; i < read.CanNum; i++)
		{
			if (read.AlnReportArr[i].AlnScore == read.score)
			{
				if (read.AlnReportArr[i].coor.bDir == false && rseq == NULL)
				{
					rseq = new char[read.rlen + 1]; rseq[read.rlen] = '\0';
					GetComplementarySeq(read.rlen, seq, rseq);
				}
				fprintf(output, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%s\t*\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read.header, read.AlnReportArr[i].iFrag, ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].name, read.AlnReportArr[i].coor.gPos, read.mapq, read.AlnReportArr[i].coor.CIGAR.c_str(), (read.AlnReportArr[i].coor.bDir? seq: rseq), read.rlen - read.score, read.score, read.sub_score);
				if (!bMultiHit) break;
			}
		}
		if (rseq != NULL)
		{
			delete[] rseq;
			rseq = NULL;
		}
	}
}

void RemoveRedundantCandidates(vector<AlignmentCandidate_t>& AlignmentVec)
{
	int thr;
	vector<AlignmentCandidate_t>::iterator iter;

	if (AlignmentVec.size() <= 1) return;
	else
	{
		int score1, score2;

		score1 = score2 = 0;
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
		{
			if (iter->Score > score2)
			{
				if (iter->Score >= score1)
				{
					score2 = score1;
					score1 = iter->Score;
				}
				else score2 = iter->Score;
			}
		}
		if (bPacBioData || score1 == score2 || score1 - score2 > 20) thr = score1;
		else thr = score2;

		//if (bDebugMode) printf("Candidate score threshold = %d\n", thr);
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++) if (iter->Score < thr) iter->Score = 0;
	}
}

bool CheckPairedAlignmentCandidates(int64_t EstiDistance, vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	int64_t dist;
	bool bPairing = false;
	int i, j, best_mate, s, num1, num2;

	num1 = (int)AlignmentVec1.size(); num2 = (int)AlignmentVec2.size();

	if (num1*num2 > 1000)
	{
		RemoveRedundantCandidates(AlignmentVec1);
		RemoveRedundantCandidates(AlignmentVec2);
	}

	for (i = 0; i != num1; i++)
	{
		if (AlignmentVec1[i].Score == 0) continue;

		for (best_mate = -1, s = 0, j = 0; j != num2; j++)
		{
			if (AlignmentVec2[j].Score == 0 || AlignmentVec2[j].PosDiff < AlignmentVec1[i].PosDiff) continue;

			dist = AlignmentVec2[j].PosDiff - AlignmentVec1[i].PosDiff;
			//printf("#%d (s=%d) and #%d (s=%d) (dist=%ld / %d)\n", i+1, AlignmentVec1[i].Score, j+1, AlignmentVec2[j].Score, dist, EstiDistance), fflush(stdout);
			if (dist < EstiDistance)
			{
				if (AlignmentVec2[j].Score > s)
				{
					best_mate = j;
					s = AlignmentVec2[j].Score;
				}
				else if (AlignmentVec2[j].Score == s) best_mate = -1;
			}
		}
		if (s > 0 && best_mate != -1)
		{
			j = best_mate;
			if (AlignmentVec2[j].PairedAlnCanIdx == -1)
			{
				bPairing = true;
				AlignmentVec1[i].PairedAlnCanIdx = j;
				AlignmentVec2[j].PairedAlnCanIdx = i;
			}
			else if (AlignmentVec1[i].Score > AlignmentVec1[AlignmentVec2[j].PairedAlnCanIdx].Score)
			{
				AlignmentVec1[AlignmentVec2[j].PairedAlnCanIdx].PairedAlnCanIdx = -1;
				AlignmentVec1[i].PairedAlnCanIdx = j;
				AlignmentVec2[j].PairedAlnCanIdx = i;
			}
		}
	}
	return bPairing;
}

void RemoveUnMatedAlignmentCandidates(vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	int i, j, num1, num2;

	num1 = (int)AlignmentVec1.size(); num2 = (int)AlignmentVec2.size();

	for (i = 0; i != num1; i++)
	{
		if (AlignmentVec1[i].PairedAlnCanIdx == -1) AlignmentVec1[i].Score = 0;
		else
		{
			j = AlignmentVec1[i].PairedAlnCanIdx;
			AlignmentVec1[i].Score = AlignmentVec2[j].Score = AlignmentVec1[i].Score + AlignmentVec2[j].Score;
		}
	}
	for (j = 0; j != num2; j++) if (AlignmentVec2[j].PairedAlnCanIdx == -1) AlignmentVec2[j].Score = 0;

	if (bDebugMode)
	{
		for (i = 0; i != num1; i++)
		{
			if ((j = AlignmentVec1[i].PairedAlnCanIdx) != -1)
				printf("#%d(s=%d) and #%d(s=%d) are pairing\n", i + 1, AlignmentVec1[i].Score, j + 1, AlignmentVec2[j].Score);
		}
	}
}

void CheckPairedFinalAlignments(ReadItem_t& read1, ReadItem_t& read2)
{
	bool bMated;
	int64_t dist;
	int i, j, best_mate, s;

	//printf("BestIdx1=%d, BestIdx2=%d\n", read1.iBestAlnCanIdx + 1, read2.iBestAlnCanIdx + 1);
	bMated = read1.AlnReportArr[read1.iBestAlnCanIdx].PairedAlnCanIdx == read2.iBestAlnCanIdx ? true : false;
	if (!bMultiHit && bMated) return;
	if (!bMated && read1.score > 0 && read2.score > 0) // identify mated pairs
	{
		for (s = 0, i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore > 0 && (j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0)
			{
				bMated = true;
				if (s < read1.AlnReportArr[i].AlnScore + read2.AlnReportArr[j].AlnScore)
				{
					s = read1.AlnReportArr[i].AlnScore + read2.AlnReportArr[j].AlnScore;
					read1.iBestAlnCanIdx = i; read1.score = read1.AlnReportArr[i].AlnScore;
					read2.iBestAlnCanIdx = j; read2.score = read2.AlnReportArr[j].AlnScore;
				}
			}
		}
	}
	if(bMated)
	{
		for (i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore != read1.score || ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore != read2.score))
			{
				read1.AlnReportArr[i].AlnScore = 0;
				read1.AlnReportArr[i].PairedAlnCanIdx = -1;
				continue;
			}
		}
	}
	else // remove all mated info
	{
		for (i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].PairedAlnCanIdx != -1) read1.AlnReportArr[i].PairedAlnCanIdx = -1;
			if (read1.AlnReportArr[i].AlnScore > 0 && read1.AlnReportArr[i].AlnScore != read1.score) read1.AlnReportArr[i].AlnScore = 0;
		}
		for (j = 0; j != read2.CanNum; j++)
		{
			if (read2.AlnReportArr[j].PairedAlnCanIdx != -1) read2.AlnReportArr[j].PairedAlnCanIdx = -1;
			if (read2.AlnReportArr[j].AlnScore > 0 && read2.AlnReportArr[j].AlnScore != read2.score) read2.AlnReportArr[j].AlnScore = 0;
		}
	}
}

void *ReadMapping(void *arg)
{
	bool bReadPairing;
	int *my_id = (int*)arg;
	ReadItem_t* ReadArr = NULL;
	int i, j, max_dist, ReadNum, EstDistance;
	vector<SeedPair_t> SeedPairVec1, SeedPairVec2;
	vector<AlignmentCandidate_t> AlignmentVec1, AlignmentVec2;

	ReadArr = new ReadItem_t[ReadChunkSize];

	while (true)
	{
		pthread_mutex_lock(&LibraryLock);
		ReadNum = GetNextChunk(bSepLibrary, ReadFileHandler1, ReadFileHandler2, ReadArr);
		fprintf(stderr, "\r%d %s reads have been processed in %ld seconds...", iTotalReadNum, (bPairEnd? "paired-end":"singled-end"), (time(NULL) - StartProcessTime));
		pthread_mutex_unlock(&LibraryLock);
		
		if (ReadNum == 0) break;
		if (bPacBioData)
		{
			for (i = 0; i != ReadNum; i++)
			{
				if (bDebugMode) printf("\n\n\nMapping pacbio read#%d %s (len=%d):\n", i + 1, ReadArr[i].header, ReadArr[i].rlen);
				
				SeedPairVec1 = IdentifySeedPairs_SensitiveMode(ReadArr[i].rlen, ReadArr[i].EncodeSeq);
				AlignmentVec1 = GenerateAlignmentCandidateForPacBioSeq(ReadArr[i].rlen, SeedPairVec1);
				RemoveRedundantCandidates(AlignmentVec1);
				if (bDebugMode) ShowAlignmentCandidateInfo(1, ReadArr[i].header, AlignmentVec1);
				GenMappingReport(true, ReadArr[i], AlignmentVec1);

				SetSingleAlignmentFlag(ReadArr[i]); EvaluateMAPQ(ReadArr[i]);
				if (bDebugMode) printf("\nEnd of mapping for read#%s (len=%d)\n%s\n", ReadArr[i].header, ReadArr[i].rlen, string().assign(100, '=').c_str());
			}
		}
		else if (bPairEnd && ReadNum % 2 == 0)
		{
			if (iPaired >= 1000)
			{
				EstDistance = (int)(iDistance / (iPaired >> 2));
				EstDistance = EstDistance + (EstDistance >> 1);
			}
			else  EstDistance = MaxInsertSize;

			//if (bDebugMode) printf("iDistance = %ld, iPaired=%d, EstiDistance=%d\n", iDistance, iPaired, EstDistance);
			for (i = 0, j = 1; i != ReadNum; i += 2, j += 2)
			{
				//if (bDebugMode) printf("Mapping paired reads#%d %s (len=%d) and %s (len=%d):\n", i + 1, ReadArr[i].header, ReadArr[i].rlen, ReadArr[j].header, ReadArr[j].rlen);

				SeedPairVec1 = IdentifySeedPairs_FastMode(ReadArr[i].rlen, ReadArr[i].EncodeSeq);
				AlignmentVec1 = GenerateAlignmentCandidateForIlluminaSeq(ReadArr[i].rlen, SeedPairVec1);

				SeedPairVec2 = IdentifySeedPairs_FastMode(ReadArr[j].rlen, ReadArr[j].EncodeSeq);
				AlignmentVec2 = GenerateAlignmentCandidateForIlluminaSeq(ReadArr[j].rlen, SeedPairVec2);

				//if (bDebugMode)
				//{
				//	ShowAlignmentCandidateInfo(1, ReadArr[i].header, AlignmentVec1);
				//	ShowAlignmentCandidateInfo(0, ReadArr[j].header, AlignmentVec2);
				//}
				bReadPairing = CheckPairedAlignmentCandidates(EstDistance, AlignmentVec1, AlignmentVec2);
				if (!bReadPairing) bReadPairing = RescueUnpairedAlignment(EstDistance, ReadArr[i], ReadArr[j], AlignmentVec1, AlignmentVec2);
				if (bReadPairing) RemoveUnMatedAlignmentCandidates(AlignmentVec1, AlignmentVec2);

				RemoveRedundantCandidates(AlignmentVec1); RemoveRedundantCandidates(AlignmentVec2);
				//if (bDebugMode)
				//{
				//	ShowAlignmentCandidateInfo(1, ReadArr[i].header, AlignmentVec1);
				//	ShowAlignmentCandidateInfo(0, ReadArr[j].header, AlignmentVec2);
				//}
				GenMappingReport(true,  ReadArr[i], AlignmentVec1);
				GenMappingReport(false, ReadArr[j], AlignmentVec2);

				CheckPairedFinalAlignments(ReadArr[i], ReadArr[j]);

				SetPairedAlignmentFlag(ReadArr[i], ReadArr[j]);
				EvaluateMAPQ(ReadArr[i]); EvaluateMAPQ(ReadArr[j]);

				//if (bDebugMode) printf("\nEnd of mapping for read#%s\n%s\n", ReadArr[i].header, string().assign(100, '=').c_str());
			}
		}
		else
		{
			for (i = 0; i != ReadNum; i++)
			{
				if (bDebugMode) printf("Mapping single read#%d %s (len=%d):\n", i + 1, ReadArr[i].header, ReadArr[i].rlen);

				SeedPairVec1 = IdentifySeedPairs_FastMode(ReadArr[i].rlen, ReadArr[i].EncodeSeq);
				AlignmentVec1 = GenerateAlignmentCandidateForIlluminaSeq(ReadArr[i].rlen, SeedPairVec1);
				RemoveRedundantCandidates(AlignmentVec1); if (bDebugMode) ShowAlignmentCandidateInfo(1, ReadArr[i].header, AlignmentVec1);
				GenMappingReport(true, ReadArr[i], AlignmentVec1);

				SetSingleAlignmentFlag(ReadArr[i]); EvaluateMAPQ(ReadArr[i]);
				
				if (bDebugMode) printf("\nEnd of mapping for read#%s\n%s\n", ReadArr[i].header, string().assign(100, '=').c_str());
			}
		}
		pthread_mutex_lock(&OutputLock);
		iTotalReadNum += ReadNum;
		if (bPairEnd && ReadNum % 2 == 0) for (i = 0, j = 1; i != ReadNum; i += 2, j += 2) OutputPairedSamFile(ReadArr[i], ReadArr[j]);
		else for (i = 0; i != ReadNum; i++) OutputSingledSamFile(ReadArr[i]);
		fflush(output);
		pthread_mutex_unlock(&OutputLock);

		for (i = 0; i != ReadNum; i++)
		{
			delete[] ReadArr[i].header;
			delete[] ReadArr[i].seq;
			delete[] ReadArr[i].EncodeSeq;
			delete[] ReadArr[i].AlnReportArr;
		}
	}
	delete[] ReadArr;
	//fprintf(stderr, "\niSimpleLen=%.2f, iNormalLen=%.2f\n", 1.0*iSimpleLength / (iSimpleLength + iNormalLength), 1.0*iNormalLength / (iSimpleLength + iNormalLength));

	return (void*)(1);
}

void Mapping()
{
	int i;
	vector<int> vec(iThreadNum); for (i = 0; i < iThreadNum; i++) vec[i] = i;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];

	for (MinSeedLength = 13; MinSeedLength < 16; MinSeedLength++) if (TwoGenomeSize < pow(4, MinSeedLength)) break;

	ReadFileHandler1.open(ReadFileName, ios_base::in);
	if (ReadFileName2 != NULL)
	{
		bSepLibrary = bPairEnd = true;
		ReadFileHandler2.open(ReadFileName2, ios_base::in);
	}
	if (!bDebugMode && SamFileName != NULL) output = fopen(SamFileName, "w");

	if (bDebugMode) iThreadNum = 1;
	else
	{
		fprintf(output, "@PG\tPN:Kart\tVN:%s\n", VersionStr);
		for (i = 0; i < iChromsomeNum; i++) fprintf(output, "@SQ\tSN:%s\tLN:%ld\n", ChromosomeVec[i].name, ChromosomeVec[i].len);
	}
	StartProcessTime = time(NULL);
	for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, ReadMapping, &vec[i]);
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
	delete[] ThreadArr;

	ReadFileHandler1.close(); if (ReadFileName2 != NULL) ReadFileHandler2.close();
	fprintf(stderr, "\rAll the %d %s reads have been processed in %lld seconds.\n", iTotalReadNum, (bPairEnd? "paired-end":"single-end"), (long long)(time(NULL) - StartProcessTime));

	if (SamFileName != NULL) fclose(output);

	if(iTotalReadNum > 0)
	{
		if (bPairEnd) fprintf(stderr, "\t# of total mapped sequences = %d (sensitivity = %.2f%%)\n\t# of paired sequences = %d (%.2f%%), average insert size = %d\n", iTotalReadNum - iUnMapped, (int)(10000 * (1.0*(iTotalReadNum - iUnMapped) / iTotalReadNum) + 0.5) / 100.0, iPaired, (int)(10000 * (1.0*iPaired / iTotalReadNum) + 0.5) / 100.0, (iPaired > 1 ? (int)(iDistance / (iPaired >> 1)) : 0));
		else fprintf(stderr, "\t# of total mapped sequences = %d (sensitivity = %.2f%%)\n", iTotalReadNum - iUnMapped, (int)(10000 * (1.0*(iTotalReadNum - iUnMapped) / iTotalReadNum) + 0.5) / 100.0);
	}
}
