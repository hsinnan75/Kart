#include <cmath>
#include "structure.h"
#include "htslib/htslib/sam.h"
#include "sam_opts.h"
#include "htslib/htslib/kseq.h"
#include "htslib/htslib/kstring.h"

#define MAPQ_COEF 30
#define Max_MAPQ  60

FILE *sam_out = 0;
samFile *bam_out = 0;
int64_t iDistance = 0;
time_t StartProcessTime;
bool bSepLibrary = false;
bam_hdr_t *header = NULL;
FILE *ReadFileHandler1, *ReadFileHandler2;
gzFile gzReadFileHandler1, gzReadFileHandler2;
static pthread_mutex_t DataLock, LibraryLock, OutputLock;
int64_t iTotalReadNum = 0, iUniqueMapping = 0, iUnMapping = 0, iPaired = 0;

static bam_hdr_t *sam_hdr_sanitise(bam_hdr_t *h) 
{
	if (!h) return NULL;
	if (h->l_text == 0) return h;

	uint32_t i;
	char *cp = h->text;
	for (i = 0; i < h->l_text; i++) {
		// NB: l_text excludes terminating nul.  This finds early ones.
		if (cp[i] == 0) break;
	}
	if (i < h->l_text) { // Early nul found.  Complain if not just padding.
		uint32_t j = i;
		while (j < h->l_text && cp[j] == '\0') j++;
	}
	return h;
}

bam_hdr_t *SamHdr2BamHdr(kstring_t *str)
{
	bam_hdr_t *h = NULL;
	h = sam_hdr_parse(str->l, str->s);
	h->l_text = str->l; h->text = str->s;

	return sam_hdr_sanitise(h);
}

void SetSingleAlignmentFlag(ReadItem_t& read)
{
	int i;

	if (read.score > read.sub_score) // unique mapping
	{
		i = read.iBestAlnCanIdx;
		if (read.AlnReportArr[i].coor.bDir == false) read.AlnReportArr[i].SamFlag = 0x10;
		else read.AlnReportArr[i].SamFlag = 0;
	}
	else if(read.score > 0)
	{
		for (i = 0; i < read.CanNum; i++)
		{
			if (read.AlnReportArr[i].AlnScore > 0)
			{
				if (read.AlnReportArr[i].coor.bDir == false) read.AlnReportArr[i].SamFlag = 0x10;
				else read.AlnReportArr[i].SamFlag = 0;
			}
		}
	}
	else read.AlnReportArr[0].SamFlag = 0x4;
}

void SetPairedAlignmentFlag(ReadItem_t& read1, ReadItem_t& read2)
{
	int i, j;

	//printf("read1:[%d, %d]:#%d, read2:[%d, %d] #%d\n", read1.score, read1.sub_score, read1.iBestAlnCanIdx, read2.score, read2.sub_score, read2.iBestAlnCanIdx); fflush(stdout);
	if (read1.score > read1.sub_score && read2.score > read2.sub_score) // unique mapping
	{
		i = read1.iBestAlnCanIdx;
		read1.AlnReportArr[i].SamFlag = 0x41; // read1 is the first read in a pair

		j = read2.iBestAlnCanIdx;
		read2.AlnReportArr[j].SamFlag = 0x81; // read2 is the second read in a pair

		if (j == read1.AlnReportArr[i].PairedAlnCanIdx) // reads are mapped in a proper pair
		{
			read1.AlnReportArr[i].SamFlag |= 0x2;
			read2.AlnReportArr[j].SamFlag |= 0x2;
		}
		read1.AlnReportArr[i].SamFlag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);
		read2.AlnReportArr[j].SamFlag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);
	}
	else
	{
		if (read1.score > read1.sub_score) // unique mapping or bMultiHit=false
		{
			i = read1.iBestAlnCanIdx;
			read1.AlnReportArr[i].SamFlag = 0x41; // read1 is the first read in a pair
			read1.AlnReportArr[i].SamFlag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);
			if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0) read1.AlnReportArr[i].SamFlag |= 0x2;// reads are mapped in a proper pair
			else read1.AlnReportArr[i].SamFlag |= 0x8; // next segment unmapped
		}
		else if(read1.score > 0)
		{
			for (i = 0; i < read1.CanNum; i++)
			{
				if (read1.AlnReportArr[i].AlnScore > 0)
				{
					read1.AlnReportArr[i].SamFlag = 0x41; // read1 is the first read in a pair
					read1.AlnReportArr[i].SamFlag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);

					if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0) read1.AlnReportArr[i].SamFlag |= 0x2;// reads are mapped in a proper pair
					else read1.AlnReportArr[i].SamFlag |= 0x8; // next segment unmapped
				}
			}
		}
		else
		{
			read1.AlnReportArr[0].SamFlag = 0x41; // read1 is the first read in a pair
			read1.AlnReportArr[0].SamFlag |= 0x4;

			if (read2.score == 0) read1.AlnReportArr[0].SamFlag |= 0x8; // next segment unmapped
			else read1.AlnReportArr[0].SamFlag |= (read2.AlnReportArr[read2.iBestAlnCanIdx].coor.bDir ? 0x10 : 0x20);
		}

		if (read2.score > read2.sub_score) // unique mapping or bMultiHit=false
		{
			j = read2.iBestAlnCanIdx;
			read2.AlnReportArr[j].SamFlag = 0x81; // read2 is the second read in a pair
			read2.AlnReportArr[j].SamFlag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);

			if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0) read2.AlnReportArr[j].SamFlag |= 0x2;// reads are mapped in a proper pair
			else read2.AlnReportArr[j].SamFlag |= 0x8; // next segment unmapped
		}
		else if (read2.score > 0)
		{
			for (j = 0; j < read2.CanNum; j++)
			{
				if (read2.AlnReportArr[j].AlnScore > 0)
				{
					read2.AlnReportArr[j].SamFlag = 0x81; // read2 is the second read in a pair
					read2.AlnReportArr[j].SamFlag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);

					if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0) read2.AlnReportArr[j].SamFlag |= 0x2;// reads are mapped in a proper pair
					else read2.AlnReportArr[j].SamFlag |= 0x8; // next segment unmapped
				}
			}
		}
		else
		{
			read2.AlnReportArr[0].SamFlag = 0x81; // read2 is the second read in a pair
			read2.AlnReportArr[0].SamFlag |= 0x4; // segment unmapped
			if (read1.score == 0) read2.AlnReportArr[0].SamFlag |= 0x8; // next segment unmapped
			else read2.AlnReportArr[0].SamFlag |= (read1.AlnReportArr[read1.iBestAlnCanIdx].coor.bDir ? 0x10 : 0x20);
		}
	}
}

void EvaluateMAPQ(ReadItem_t& read)
{
	if (read.score == 0 || read.score == read.sub_score) read.mapq = 0;
	else
	{
		if (bPacBioData)
		{
			float fScale = 85.0*(int)(ceil(read.rlen / 100 + 0.5));
			if (fScale > 2000) fScale = 2000;
			read.mapq = (int)(Max_MAPQ * (read.score / fScale));
		}
		else if (read.sub_score == 0 || read.score - read.sub_score > 5) read.mapq = Max_MAPQ;
		else read.mapq = (int)(MAPQ_COEF * (1 - (float)(read.score - read.sub_score) / read.score)*log(read.score) + 0.4999);
		if (read.mapq > Max_MAPQ) read.mapq = Max_MAPQ;
	}
}

void OutputPairedAlignments(ReadItem_t& read1, ReadItem_t& read2, char* buffer, int& myUniqueMapping, int& myUnMapping, vector<string>& SamOutputVec)
{
	string rqual;
	char *seq, *rseq;
	int i, j, len, dist = 0;

	if (read1.score == 0)
	{
		myUnMapping++;
		len = sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0", read1.header, read1.AlnReportArr[0].SamFlag, read1.seq, (FastQFormat ? read1.qual : "*"));
		SamOutputVec.push_back(buffer);
	}
	else
	{
		if (read1.mapq == Max_MAPQ) myUniqueMapping++;

		seq = read1.seq; rseq = NULL;
		for (i = read1.iBestAlnCanIdx; i < read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore > 0)
			{
				if (read1.AlnReportArr[i].coor.bDir == false && rseq == NULL)
				{
					rseq = new char[read1.rlen + 1]; rseq[read1.rlen] = '\0'; GetComplementarySeq(read1.rlen, seq, rseq);
					if (FastQFormat)
					{
						rqual = read1.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0)
				{
					dist = (int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen));
					if (i == read1.iBestAlnCanIdx)
					{
						iPaired += 2;
						if (abs(dist) < 10000) iDistance += abs(dist);
					}
					len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read1.header, read1.AlnReportArr[i].SamFlag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), read2.AlnReportArr[j].coor.gPos, dist, (read1.AlnReportArr[i].coor.bDir ? seq : rseq), (FastQFormat ? (read1.AlnReportArr[i].coor.bDir ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.score, read1.score, read1.sub_score);
				}
				else
				{
					len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read1.header, read1.AlnReportArr[i].SamFlag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (read1.AlnReportArr[i].coor.bDir ? seq : rseq), (FastQFormat ? (read1.AlnReportArr[i].coor.bDir ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.score, read1.score, read1.sub_score);
				}
				SamOutputVec.push_back(buffer);
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
		myUnMapping++;
		len = sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0", read2.header, read2.AlnReportArr[0].SamFlag, read2.seq, (FastQFormat ? read2.qual : "*"));
		SamOutputVec.push_back(buffer);
	}
	else
	{
		if (read2.mapq == Max_MAPQ) myUniqueMapping++;

		rseq = read2.seq; seq = NULL;
		for (j = read2.iBestAlnCanIdx; j < read2.CanNum; j++)
		{
			if (read2.AlnReportArr[j].AlnScore > 0)
			{
				if (read2.AlnReportArr[j].coor.bDir == true && seq == NULL)
				{
					seq = new char[read2.rlen + 1]; seq[read2.rlen] = '\0'; GetComplementarySeq(read2.rlen, rseq, seq);
					if (FastQFormat)
					{
						rqual = read2.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0)
				{
					dist = 0 - ((int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen)));
					len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read2.header, read2.AlnReportArr[j].SamFlag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), read1.AlnReportArr[i].coor.gPos, dist, (read2.AlnReportArr[j].coor.bDir ? seq : rseq), (FastQFormat ? (read2.AlnReportArr[j].coor.bDir ? rqual.c_str() : read2.qual) : "*"), read2.rlen - read2.score, read2.score, read2.sub_score);
				}
				else len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read2.header, read2.AlnReportArr[j].SamFlag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (read2.AlnReportArr[j].coor.bDir ? seq : rseq), (FastQFormat ? (read2.AlnReportArr[j].coor.bDir ? rqual.c_str() : read2.qual) : "*"), read2.rlen - read2.score, read2.score, read2.sub_score);
				SamOutputVec.push_back(buffer);
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

void OutputSingledAlignments(ReadItem_t& read, char* buffer, int& myUniqueMapping, int& myUnMapping, vector<string>& SamOutputVec)
{
	int len;
	string rqual;

	if (read.score == 0)
	{
		myUnMapping++;
		len = sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0", read.header, read.AlnReportArr[0].SamFlag, read.seq, (FastQFormat ? read.qual : "*"));
		buffer[len] = '\0'; SamOutputVec.push_back(buffer);
	}
	else
	{
		int i;
		char *seq, *rseq;

		if (read.mapq == Max_MAPQ) myUniqueMapping++;

		seq = read.seq; rseq = NULL;
		for (i = read.iBestAlnCanIdx; i < read.CanNum; i++)
		{
			if (read.AlnReportArr[i].AlnScore == read.score)
			{
				if (read.AlnReportArr[i].coor.bDir == false && rseq == NULL)
				{
					rseq = new char[read.rlen + 1]; rseq[read.rlen] = '\0'; GetComplementarySeq(read.rlen, seq, rseq);
					if (FastQFormat)
					{
						rqual = read.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read.header, read.AlnReportArr[i].SamFlag, ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].name, read.AlnReportArr[i].coor.gPos, read.mapq, read.AlnReportArr[i].coor.CIGAR.c_str(), (read.AlnReportArr[i].coor.bDir? seq: rseq), (FastQFormat ? (read.AlnReportArr[i].coor.bDir ? read.qual : rqual.c_str()) : "*"), read.rlen - read.score, read.score, read.sub_score);
				buffer[len] = '\0'; SamOutputVec.push_back(buffer);

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
			//printf("#%d (s=%d) and #%d (s=%d) (dist=%lld / %d)\n", i+1, AlignmentVec1[i].Score, j+1, AlignmentVec2[j].Score, dist, EstiDistance), fflush(stdout);
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
	int i, j, s;

	//printf("BestIdx1=%d, BestIdx2=%d\n", read1.iBestAlnCanIdx + 1, read2.iBestAlnCanIdx + 1);
	if (read1.iBestAlnCanIdx != -1 && read2.iBestAlnCanIdx != -1) bMated = read1.AlnReportArr[read1.iBestAlnCanIdx].PairedAlnCanIdx == read2.iBestAlnCanIdx ? true : false;
	else bMated = false;

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
	char* buffer;
	bool bReadPairing;
	ReadItem_t* ReadArr = NULL;
	vector<string> SamOutputVec;
	vector<SeedPair_t> SeedPairVec1, SeedPairVec2;
	vector<AlignmentCandidate_t> AlignmentVec1, AlignmentVec2;
	int i, j, n, ReadNum, EstDistance, myUniqueMapping, myUnMapping;

	ReadArr = new ReadItem_t[ReadChunkSize];

	if (bPacBioData) buffer = new char[1024000];
	else buffer = new char[10240];

	while (true)
	{
		pthread_mutex_lock(&LibraryLock);
		if(gzCompressed) ReadNum = gzGetNextChunk(bSepLibrary, gzReadFileHandler1, gzReadFileHandler2, ReadArr);
		else ReadNum = GetNextChunk(bSepLibrary, ReadFileHandler1, ReadFileHandler2, ReadArr);
		fprintf(stdout, "\r%lld %s reads have been processed in %ld seconds...", (long long)iTotalReadNum, (bPairEnd ? "paired-end" : "singled-end"), (long)(time(NULL) - StartProcessTime)); fflush(stdout);
		pthread_mutex_unlock(&LibraryLock);
		
		if (ReadNum == 0) break;
		if (bPacBioData)
		{
			for (i = 0; i != ReadNum; i++)
			{
				if (bDebugMode) printf("\n\n\nMapping pacbio read#%d %s (len=%d):\n", i + 1, ReadArr[i].header, ReadArr[i].rlen);
				
				SeedPairVec1 = IdentifySeedPairs_SensitiveMode(ReadArr[i].rlen, ReadArr[i].EncodeSeq);
				AlignmentVec1 = GenerateAlignmentCandidateForPacBioSeq(ReadArr[i].rlen, SeedPairVec1);
				//if (bDebugMode) ShowAlignmentCandidateInfo(1, ReadArr[i].header, AlignmentVec1);
				RemoveRedundantCandidates(AlignmentVec1);
				if (bDebugMode) ShowAlignmentCandidateInfo(1, ReadArr[i].header, AlignmentVec1);
				GenMappingReport(true, ReadArr[i], AlignmentVec1);

				SetSingleAlignmentFlag(ReadArr[i]); EvaluateMAPQ(ReadArr[i]);
				//if (bDebugMode) printf("\nEnd of mapping for read#%s (len=%d)\n%s\n", ReadArr[i].header, ReadArr[i].rlen, string().assign(100, '=').c_str());
			}
		}
		else if (bPairEnd && ReadNum % 2 == 0)
		{
			pthread_mutex_lock(&DataLock);
			if (iPaired >= 1000)
			{
				EstDistance = (int)(iDistance / (iPaired >> 2));
				EstDistance = EstDistance + (EstDistance >> 1);
			}
			else  EstDistance = MaxInsertSize;
			pthread_mutex_unlock(&DataLock);
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
		myUniqueMapping = myUnMapping = 0; SamOutputVec.clear();
		if (bPairEnd && ReadNum % 2 == 0) for (i = 0; i != ReadNum; i += 2) OutputPairedAlignments(ReadArr[i], ReadArr[i+1], buffer, myUniqueMapping, myUnMapping, SamOutputVec);
		else for (i = 0; i != ReadNum; i++) OutputSingledAlignments(ReadArr[i], buffer, myUniqueMapping, myUnMapping, SamOutputVec);

		pthread_mutex_lock(&OutputLock);
		iTotalReadNum += ReadNum; iUniqueMapping += myUniqueMapping; iUnMapping += myUnMapping;
		if (OutputFileFormat == 0)
		{
			for (vector<string>::iterator iter = SamOutputVec.begin(); iter != SamOutputVec.end(); iter++)
			{
				fprintf(sam_out, "%s\n", iter->c_str()); fflush(sam_out);
			}
		}
		else
		{
			bam1_t *b = bam_init1();
			kstring_t str = { 0, 0, NULL };
			for (vector<string>::iterator iter = SamOutputVec.begin(); iter != SamOutputVec.end(); iter++)
			{
				//gzwrite(gzOutput, iter->c_str(), iter->length());
				str.s = (char*)iter->c_str(); str.l = iter->length();
				if (sam_parse1(&str, header, b) >= 0) sam_write1(bam_out, header, b);
			}
			bam_destroy1(b);
		}
		pthread_mutex_unlock(&OutputLock);

		for (i = 0; i != ReadNum; i++)
		{
			delete[] ReadArr[i].header;
			delete[] ReadArr[i].seq;
			if (FastQFormat) delete[] ReadArr[i].qual;
			delete[] ReadArr[i].EncodeSeq;
			if(ReadArr[i].CanNum > 0) delete[] ReadArr[i].AlnReportArr;
		}
	}
	delete[] buffer;
	delete[] ReadArr;
	//fprintf(stdout, "\niSimpleLen=%.2f, iNormalLen=%.2f\n", 1.0*iSimpleLength / (iSimpleLength + iNormalLength), 1.0*iNormalLength / (iSimpleLength + iNormalLength));

	return (void*)(1);
}

void Mapping()
{
	int i;
	vector<int> vec(iThreadNum); for (i = 0; i < iThreadNum; i++) vec[i] = i;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];

	for (MinSeedLength = 13; MinSeedLength < 16; MinSeedLength++) if (TwoGenomeSize < pow(4, MinSeedLength)) break;

	if (bDebugMode) iThreadNum = 1;
	else
	{
		int len;
		char buffer[1024];
		kstring_t str = { 0, 0, NULL };

		if (OutputFileFormat == 0) sam_out = fopen(OutputFileName, "w");
		else bam_out = sam_open_format(OutputFileName, "wb", NULL);

		len = sprintf(buffer, "@PG\tID:kart\tPN:Kart\tVN:%s\n", VersionStr);

		if (OutputFileFormat == 0) fprintf(sam_out, "%s", buffer);
		else kputsn(buffer, len, &str);

		for (i = 0; i < iChromsomeNum; i++)
		{
			len = sprintf(buffer, "@SQ\tSN:%s\tLN:%lld\n", ChromosomeVec[i].name, ChromosomeVec[i].len);

			if (OutputFileFormat == 0) fprintf(sam_out, "%s", buffer);
			else kputsn(buffer, len, &str);
		}
		if (OutputFileFormat == 1)
		{
			header = SamHdr2BamHdr(&str);
			sam_hdr_write(bam_out, header);
		}
	}

	StartProcessTime = time(NULL);
	for (int LibraryID = 0; LibraryID < (int)ReadFileNameVec1.size(); LibraryID++)
	{
		gzReadFileHandler1 = gzReadFileHandler2 = NULL; ReadFileHandler1 = ReadFileHandler2 = NULL;

		if (ReadFileNameVec1[LibraryID].substr(ReadFileNameVec1[LibraryID].find_last_of('.') + 1) == "gz") gzCompressed = true;
		else gzCompressed = false;

		FastQFormat = CheckReadFormat(ReadFileNameVec1[LibraryID].c_str());
		//fprintf(stdout, "gz=%s, format=%s\n", gzCompressed ? "Yes" : "No", FastQFormat ? "Fastq" : "Fasta");

		if (gzCompressed) gzReadFileHandler1 = gzopen(ReadFileNameVec1[LibraryID].c_str(), "rb");
		else ReadFileHandler1 = fopen(ReadFileNameVec1[LibraryID].c_str(), "r");

		if (ReadFileNameVec1.size() == ReadFileNameVec2.size())
		{
			bSepLibrary = bPairEnd = true;
			if (FastQFormat == CheckReadFormat(ReadFileNameVec2[LibraryID].c_str()))
			{
				if (gzCompressed) gzReadFileHandler2 = gzopen(ReadFileNameVec2[LibraryID].c_str(), "rb");
				else ReadFileHandler2 = fopen(ReadFileNameVec2[LibraryID].c_str(), "r");
			}
			else
			{
				fprintf(stdout, "Error! %s and %s are with different format...\n", (char*)ReadFileNameVec1[LibraryID].c_str(), (char*)ReadFileNameVec2[LibraryID].c_str());
				continue;
			}
		}
		else bSepLibrary = false;

		if (ReadFileHandler1 == NULL && gzReadFileHandler1 == NULL) continue;
		if (bSepLibrary && ReadFileHandler2 == NULL && gzReadFileHandler2 == NULL) continue;

		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, ReadMapping, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

		if (gzCompressed)
		{
			if (gzReadFileHandler1 != NULL) gzclose(gzReadFileHandler1);
			if (gzReadFileHandler2 != NULL) gzclose(gzReadFileHandler2);
		}
		else
		{
			if (ReadFileHandler1 != NULL) fclose(ReadFileHandler1);
			if (ReadFileHandler2 != NULL) fclose(ReadFileHandler2);
		}
	}
	fprintf(stdout, "\rAll the %lld %s reads have been processed in %lld seconds.\n", (long long)iTotalReadNum, (bPairEnd? "paired-end":"single-end"), (long long)(time(NULL) - StartProcessTime));
	delete[] ThreadArr;

	if (OutputFileFormat == 0) fclose(sam_out);
	else sam_close(bam_out);

	if(iTotalReadNum > 0)
	{
		if (bPairEnd) fprintf(stdout, "\t# of total mapped sequences = %lld (sensitivity = %.2f%%)\n\t# of paired sequences = %lld (%.2f%%), average insert size = %d\n", (long long)(iTotalReadNum - iUnMapping), (int)(10000 * (1.0*(iTotalReadNum - iUnMapping) / iTotalReadNum) + 0.5) / 100.0, (long long)iPaired, (int)(10000 * (1.0*iPaired / iTotalReadNum) + 0.5) / 100.0, (iPaired > 1 ? (int)(iDistance / (iPaired >> 1)) : 0));
		else fprintf(stdout, "\t# of total mapped sequences = %lld (sensitivity = %.2f%%)\n", (long long)(iTotalReadNum - iUnMapping), (int)(10000 * (1.0*(iTotalReadNum - iUnMapping) / iTotalReadNum) + 0.5) / 100.0);
		fprintf(stdout, "Alignment output: %s\n", OutputFileName);
	}
}
