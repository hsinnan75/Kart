#include "structure.h"

int IdentifyMaxAlignmentCandidateScore(vector<AlignmentCandidate_t>& AlignmentVec)
{
	int score = 0;
	vector<AlignmentCandidate_t>::iterator iter;

	for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
		if (iter->Score > score) score = iter->Score;

	return score;
}

int DetermineAnchorThreshold(vector<AlignmentCandidate_t>& AlignmentVec)
{
	int thr = 0;

	for (vector<AlignmentCandidate_t>::const_iterator iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
		if (iter->Score > thr) thr = iter->Score;

	if ((thr -= 30) < 50) thr = 50;

	return thr;
}

AlignmentCandidate_t IdnetifyRescueCandidate(int rlen, uint64_t gPos, vector<SeedPair_t>& vec)
{
	int i, j, num, s;
	vector<SeedPair_t> SeedPairVec;
	AlignmentCandidate_t AlignmentCandidate;

	//if (bDebugMode) printf("[RescueCandidate]Raw seeds from common kmers:\n"), ShowSeedInfo(vec);

	AlignmentCandidate.Score = s = 0; AlignmentCandidate.PairedAlnCanIdx = -1;
	for (num = (int)vec.size(), i = 0; i < num;)
	{
		vec[i].gPos += gPos; s = vec[i].rLen;
		SeedPairVec.resize(1); SeedPairVec[0] = vec[i];
		//if (bDebugMode) printf("\nMaster seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", vec[i].rPos, vec[i].rPos + vec[i].rLen - 1, vec[i].gPos, vec[i].gPos + vec[i].gLen - 1, vec[i].rLen, vec[i].PosDiff);
		for (j = i + 1; j < num; j++)
		{
			if (vec[j].PosDiff - vec[i].PosDiff < MaxGaps)
			{
				//if (bDebugMode) printf("add seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", vec[j].rPos, vec[j].rPos + vec[j].rLen - 1, vec[j].gPos, vec[j].gPos + vec[j].gLen - 1, vec[j].rLen, vec[i].PosDiff);
				vec[j].gPos += gPos; s += vec[j].rLen;
				SeedPairVec.push_back(vec[j]);
			}
			else break;
		}
		//s = CalculateSeedCoverage(rlen, SeedPairVec);
		if (s > AlignmentCandidate.Score)
		{
			AlignmentCandidate.Score = s;
			AlignmentCandidate.PosDiff = SeedPairVec[0].PosDiff + gPos;
			AlignmentCandidate.SeedVec = SeedPairVec;
		}
		i = j;
	}
	sort(AlignmentCandidate.SeedVec.begin(), AlignmentCandidate.SeedVec.end(), CompByGenomePos);

	for (vector<SeedPair_t>::iterator iter = AlignmentCandidate.SeedVec.begin(); iter != AlignmentCandidate.SeedVec.end(); iter++)
		iter->PosDiff += gPos;

	if (bDebugMode && AlignmentCandidate.Score > 0)
	{
		printf("\n\nCandidate score = %d\n", AlignmentCandidate.Score);
		ShowSeedLocationInfo(AlignmentCandidate.PosDiff);
		ShowSeedInfo(AlignmentCandidate.SeedVec);
	}
	return AlignmentCandidate;
}

bool RescueUnpairedAlignment(int EstDistance, ReadItem_t& r1, ReadItem_t& r2, vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	char* RefSeg;
	int iFixStrategy; // 1: rescue with r1, 2: rescue with r2, 3: rescue with r1 and r2, 4: give up
	bool bMated = false;
	vector<KmerItem_t> vec1, vec2;
	vector<KmerPair_t> KmerPairVec;
	vector<SeedPair_t> SimplePairVec;
	int64_t left_boundary, right_boundary;
	AlignmentCandidate_t AlignmentCandidate;
	int i, j, score1, score2, thr, num1, num2, slen;

	score1 = IdentifyMaxAlignmentCandidateScore(AlignmentVec1);
	score2 = IdentifyMaxAlignmentCandidateScore(AlignmentVec2);

	if (score1 == 0 && score2 == 0) return false;
	else if (score1 < (int)(r1.rlen*0.1) && score2 < (int)(r2.rlen*0.1)) iFixStrategy = 4;
	else if (score1 > score2 && score1 - score2 > 50) iFixStrategy = 1;
	else if (score2 > score1 && score2 - score1 > 50) iFixStrategy = 2;
	else iFixStrategy = 3;

	if (EstDistance > MaxInsertSize) EstDistance = MaxInsertSize;
	if (bDebugMode) printf("\n\nStart FixUnpairedAlignment with strategy %d (%d vs %d) and EsitDistance=%d\n\n", iFixStrategy, score1, score2, EstDistance), fflush(stdout);

	num1 = (int)AlignmentVec1.size(); num2 = (int)AlignmentVec2.size();
	if (iFixStrategy == 1 || iFixStrategy == 3)
	{
		// identify r2's alignment with r1's alignment region
		thr = DetermineAnchorThreshold(AlignmentVec1); vec1 = CreateKmerVecFromReadSeq(r2.rlen, r2.seq);
		for (j = num2, i = 0; i < num1; i++)
		{
			if (AlignmentVec1[i].Score < thr) continue;

			left_boundary  = AlignmentVec1[i].PosDiff;
			right_boundary = AlignmentVec1[i].PosDiff + EstDistance + r2.rlen;
			if (right_boundary > TwoGenomeSize) right_boundary = TwoGenomeSize - 1;

			if (ChrLocMap.lower_bound(left_boundary)->second != ChrLocMap.lower_bound(right_boundary)->second)
				right_boundary = ChrLocMap.lower_bound(left_boundary)->first;

			if ((slen = (int)(right_boundary - left_boundary)) < r2.rlen) continue;
			if (bDebugMode) printf("\n\nAnchor1-Candidate#%d (Score=%d) pos=%lld, Search region = [%lld - %lld], len = %d\n\n", i + 1, AlignmentVec1[i].Score, AlignmentVec1[i].PosDiff, left_boundary, right_boundary, slen), fflush(stdout);

			RefSeg  = RefSequence + left_boundary; vec2 = CreateKmerVecFromReadSeq(slen, RefSeg);
			KmerPairVec = IdentifyCommonKmers(slen, vec1, vec2);
			SimplePairVec = GenerateSimplePairsFromCommonKmers(10, KmerPairVec);
			AlignmentCandidate = IdnetifyRescueCandidate(r2.rlen, left_boundary, SimplePairVec);

			if (AlignmentCandidate.Score > score2)
			{
				bMated = true;
				AlignmentCandidate.PairedAlnCanIdx = i;
				AlignmentVec1[i].PairedAlnCanIdx = j++;
				AlignmentVec2.push_back(AlignmentCandidate);
			}
		}
	}
	if(iFixStrategy == 2 || iFixStrategy == 3)
	{
		// identify r1's alignment with r2's alignment region
		thr = DetermineAnchorThreshold(AlignmentVec2); vec1 = CreateKmerVecFromReadSeq(r1.rlen, r1.seq);
		for (i = num1, j = 0; j < num2; j++)
		{
			if (AlignmentVec2[j].Score < thr) continue;

			left_boundary = AlignmentVec2[j].PosDiff - EstDistance;
			right_boundary = AlignmentVec2[j].PosDiff + r2.rlen;

			if (right_boundary > TwoGenomeSize) right_boundary = TwoGenomeSize - 1;

			if (ChrLocMap.lower_bound(left_boundary)->second != ChrLocMap.lower_bound(right_boundary)->second)
				left_boundary = ChrLocMap.lower_bound(right_boundary)->first - ChromosomeVec[ChrLocMap.lower_bound(right_boundary)->second].len + 1;

			if ((slen = (int)(right_boundary - left_boundary)) < r1.rlen) continue;
			if (bDebugMode) printf("\n\nAnchor2-Candidate#%d (Score=%d) pos=%lld, Search region = [%lld - %lld], len = %d\n\n", i + 1, AlignmentVec2[i].Score, AlignmentVec2[i].PosDiff, left_boundary, right_boundary, slen);

			RefSeg = RefSequence + left_boundary; vec2 = CreateKmerVecFromReadSeq(slen, RefSeg);
			KmerPairVec = IdentifyCommonKmers(slen, vec1, vec2);
			SimplePairVec = GenerateSimplePairsFromCommonKmers(10, KmerPairVec);
			AlignmentCandidate = IdnetifyRescueCandidate(r2.rlen, left_boundary, SimplePairVec);

			if (AlignmentCandidate.Score > score1)
			{
				bMated = true;
				AlignmentCandidate.PairedAlnCanIdx = j;
				AlignmentVec2[j].PairedAlnCanIdx = i++;
				AlignmentVec1.push_back(AlignmentCandidate);
			}
		}
	}
	return bMated;
}
