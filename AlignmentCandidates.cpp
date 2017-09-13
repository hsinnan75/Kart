#include "structure.h"

void ShowFragmentPair(char* ReadSeq, SeedPair_t& sp)
{
	string frag1, frag2;
	frag1.resize(sp.rLen); strncpy((char*)frag1.c_str(), ReadSeq + sp.rPos, sp.rLen);
	frag2.resize(sp.gLen); strncpy((char*)frag2.c_str(), RefSequence + sp.gPos, sp.gLen);
	printf("FragmentPair:\n%s #read[%d-%d]=%d\n%s #chr[%lld-%lld]=%d\nscore = %d\n\n", frag1.c_str(), sp.rPos, sp.rPos + sp.rLen - 1, sp.rLen, frag2.c_str(), sp.gPos, sp.gPos + sp.gLen - 1, sp.gLen, sp.rLen), fflush(stdout);
}

bool CompByPosDiff(const SeedPair_t& p1, const SeedPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.rPos < p2.rPos;
	else return p1.PosDiff < p2.PosDiff;
}

bool CompByGenomePos(const SeedPair_t& p1, const SeedPair_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.rPos < p2.rPos;
	else return p1.gPos < p2.gPos;
}

bool CompByReadPos(const SeedPair_t& p1, const SeedPair_t& p2)
{
	return p1.rPos < p2.rPos;
}

bool CompByFirstInt(const pair<int, int>& p1, const pair<int, int>& p2)
{
	return p1.first < p2.first;
}

bool CheckCandidateValidity(vector<SeedPair_t>& SeedPairVec)
{
	int i, j, num;
	bool bValidity = true;

	for (num = (int)SeedPairVec.size(), i = 0, j = 1; j < num; i++, j++)
	{
		if ((SeedPairVec[i].rLen > 0 && SeedPairVec[j].rLen > 0 && SeedPairVec[i].rPos + SeedPairVec[i].rLen - SeedPairVec[j].rPos > 0) || (SeedPairVec[i].gLen > 0 && SeedPairVec[j].gLen > 0 && SeedPairVec[i].gPos + SeedPairVec[i].gLen - SeedPairVec[j].gPos > 0))
		{
			bValidity = false;
			break;
		}
	}
	return bValidity;
}

vector<SeedPair_t> IdentifySeedPairs_FastMode(int rlen, uint8_t* EncodeSeq)
{
	SeedPair_t SeedPair;
	int i, pos, end_pos;
	vector<SeedPair_t> SeedPairVec;
	bwtSearchResult_t bwtSearchResult;

	SeedPair.bSimple = true; pos = 0, end_pos = rlen - MinSeedLength;
	while (pos < end_pos)
	{
		if (EncodeSeq[pos] > 3) pos++;
		else
		{
			bwtSearchResult = BWT_Search(EncodeSeq, pos, rlen);
			//if (bDebugMode) printf("Pos=%d, Freq=%d, Len=%d\n", pos, bwtSearchResult.freq, bwtSearchResult.len);
			if (bwtSearchResult.freq > 0)
			{
				SeedPair.rPos = pos; SeedPair.rLen = SeedPair.gLen = bwtSearchResult.len;
				for (i = 0; i != bwtSearchResult.freq; i++)
				{
					SeedPair.PosDiff = (SeedPair.gPos = bwtSearchResult.LocArr[i]) - SeedPair.rPos;
					SeedPairVec.push_back(SeedPair);
				}
				delete[] bwtSearchResult.LocArr;
			}
			pos += (bwtSearchResult.len + 1);
		}
	}
	sort(SeedPairVec.begin(), SeedPairVec.end(), CompByPosDiff);

	return SeedPairVec;
}

vector<AlignmentCandidate_t> GenerateAlignmentCandidateForIlluminaSeq(int rlen, vector<SeedPair_t> SeedPairVec)
{
	int i, j, k, thr, num;
	AlignmentCandidate_t AlignmentCandidate;
	vector<AlignmentCandidate_t> AlignmentVec;

	if ((thr = (int)(rlen*0.2)) > 50) thr = 50;

	//if (bDebugMode) printf("\n\nRaw seeds:\n"), ShowSeedInfo(SeedPairVec);
	AlignmentCandidate.PairedAlnCanIdx = -1; num = (int)SeedPairVec.size();

	i = 0; while (i < num && SeedPairVec[i].PosDiff < 0) i++;
	for (; i < num;)
	{
		AlignmentCandidate.Score = SeedPairVec[i].rLen;
		AlignmentCandidate.SeedVec.resize(1); AlignmentCandidate.SeedVec[0] = SeedPairVec[i];
		//if (bDebugMode) printf("\nMaster seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", SeedPairVec[i].rPos, SeedPairVec[i].rPos + SeedPairVec[i].rLen - 1, SeedPairVec[i].gPos, SeedPairVec[i].gPos + SeedPairVec[i].gLen - 1, SeedPairVec[i].rLen, SeedPairVec[i].PosDiff);

		for (j = i, k = i + 1; k < num; k++)
		{
			if (SeedPairVec[k].PosDiff == SeedPairVec[j].PosDiff || SeedPairVec[k].PosDiff - SeedPairVec[j].PosDiff < MaxGaps)
			{
				//if (bDebugMode) printf("add seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", SeedPairVec[k].rPos, SeedPairVec[k].rPos + SeedPairVec[k].rLen - 1, SeedPairVec[k].gPos, SeedPairVec[k].gPos + SeedPairVec[k].gLen - 1, SeedPairVec[k].rLen, SeedPairVec[k].PosDiff);
				AlignmentCandidate.Score += SeedPairVec[k].rLen;
				AlignmentCandidate.SeedVec.push_back(SeedPairVec[k]);
				j = k;
			}
			else break;
		}
		if (AlignmentCandidate.Score > thr)
		{
			if (AlignmentCandidate.Score - 50 > thr) thr = AlignmentCandidate.Score - 50;
			AlignmentCandidate.PosDiff = AlignmentCandidate.SeedVec[0].PosDiff;
			if (AlignmentCandidate.PosDiff < 0) AlignmentCandidate.PosDiff = 0;
			sort(AlignmentCandidate.SeedVec.begin(), AlignmentCandidate.SeedVec.end(), CompByGenomePos);
			AlignmentVec.push_back(AlignmentCandidate);
			//if (bDebugMode)
			//{
			//	printf("Candidate score = %d\n", AlignmentCandidate.Score);
			//	ShowSeedLocationInfo(AlignmentCandidate.SeedVec[0]);
			//	ShowSeedInfo(AlignmentCandidate.SeedVec);
			//}
		}
		i = k;
	}
	return AlignmentVec;
}

vector<SeedPair_t> IdentifySeedPairs_SensitiveMode(int rlen, uint8_t* EncodeSeq)
{
	SeedPair_t SeedPair;
	int i, pos, stop_pos, end_pos;
	vector<SeedPair_t> SeedPairVec;
	bwtSearchResult_t bwtSearchResult;

	SeedPair.bSimple = true; pos = 0, stop_pos = 100; end_pos = rlen - MinSeedLength;
	while (pos < end_pos)
	{
		if (EncodeSeq[pos] > 3) pos++, stop_pos++;
		else
		{
			bwtSearchResult = BWT_Search(EncodeSeq, pos, stop_pos);
			//if (bDebugMode) printf("Pos=%d, Freq=%d, Len=%d\n", pos, bwtSearchResult.freq, bwtSearchResult.len);
			if (bwtSearchResult.freq > 0)
			{
				SeedPair.rPos = pos; SeedPair.rLen = SeedPair.gLen = bwtSearchResult.len;
				for (i = 0; i != bwtSearchResult.freq; i++)
				{
					SeedPair.PosDiff = (SeedPair.gPos = bwtSearchResult.LocArr[i]) - SeedPair.rPos;
					SeedPairVec.push_back(SeedPair);
				}
				delete[] bwtSearchResult.LocArr;

				pos += 30; stop_pos += 30;
				//pos += MinSeedLength; stop_pos += MinSeedLength;
			}
			else
			{
				pos += MinSeedLength; stop_pos += MinSeedLength;
			}
			if (stop_pos > rlen) stop_pos = rlen;
		}
	}
	sort(SeedPairVec.begin(), SeedPairVec.end(), CompByGenomePos);

	return SeedPairVec;
}

vector<AlignmentCandidate_t> GenerateAlignmentCandidateForPacBioSeq(int rlen, vector<SeedPair_t> SeedPairVec)
{
	bool *TakenArr;
	int i, j, k, thr, num;
	AlignmentCandidate_t AlignmentCandidate;
	vector<AlignmentCandidate_t> AlignmentVec;

	if ((num = (int)SeedPairVec.size()) > 0)
	{
		//thr = (int)(rlen*0.05);
		thr = 0;
		//printf("Raw seeds:\n"); ShowSeedInfo(SeedPairVec);
		AlignmentCandidate.PairedAlnCanIdx = -1; TakenArr = new bool[num]();
		i = 0; while (i < num && SeedPairVec[i].PosDiff < 0) i++;
		for (; i < num; i++)
		{
			if (TakenArr[i]) continue;
			AlignmentCandidate.Score = SeedPairVec[i].rLen; TakenArr[i] = true;
			AlignmentCandidate.SeedVec.clear(); AlignmentCandidate.SeedVec.push_back(SeedPairVec[i]);
			//if (bDebugMode) printf("Master seed: r[%d-%d] g[%lld-%lld], len=%d\n", SeedPairVec[i].rPos, SeedPairVec[i].rPos + SeedPairVec[i].rLen - 1, SeedPairVec[i].gPos, SeedPairVec[i].gPos + SeedPairVec[i].gLen - 1, SeedPairVec[i].rLen);
			for (j = i, k = i + 1; k < num; k++)
			{
				if (TakenArr[k]) continue;
				if (abs(SeedPairVec[k].PosDiff - SeedPairVec[j].PosDiff) < 300)
				{
					if (SeedPairVec[k].rPos > SeedPairVec[j].rPos)
					{
						//if (bDebugMode) printf("add seed: r[%d-%d] g[%lld-%lld], len=%d\n", SeedPairVec[k].rPos, SeedPairVec[k].rPos + SeedPairVec[k].rLen - 1, SeedPairVec[k].gPos, SeedPairVec[k].gPos + SeedPairVec[k].gLen - 1, SeedPairVec[k].rLen);
						AlignmentCandidate.Score += SeedPairVec[k].rLen;
						AlignmentCandidate.SeedVec.push_back(SeedPairVec[k]);
						TakenArr[(j = k)] = true;
					}
				}
				else if (SeedPairVec[k].gPos - SeedPairVec[j].gPos > 1000) break;
			}
			if (AlignmentCandidate.Score >= thr)
			{
				thr = AlignmentCandidate.Score;
				AlignmentCandidate.PosDiff = SeedPairVec[i].PosDiff;
				if (AlignmentCandidate.PosDiff < 0) AlignmentCandidate.PosDiff = 0;

				AlignmentVec.push_back(AlignmentCandidate);
				//if (bDebugMode)
				//{
				//	printf("\n\nCandidate score = %d\n", AlignmentCandidate.Score);
				//	ShowSeedLocationInfo(AlignmentCandidate.PosDiff);
				//	ShowSeedInfo(AlignmentCandidate.SeedVec);
				//}
			}
		}
		delete[] TakenArr;
	}
	return AlignmentVec;
}

void RemoveNullSeeds(vector<SeedPair_t>& SeedVec)
{
	for (vector<SeedPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end();)
	{
		if (iter->rLen == 0) iter = SeedVec.erase(iter);
		else iter++;
	}
}

void RemoveTandemRepeatSeeds(vector<SeedPair_t>& SeedVec)
{
	int i, j, k, num;
	bool bTandemRepeat = false;
	vector<pair<int, int> > vec;

	num = (int)SeedVec.size(); if (num < 2) return;
	vec.resize(num);
	for (i = 0; i < num; i++) vec[i] = make_pair(SeedVec[i].rPos, i);
	sort(vec.begin(), vec.end(), CompByFirstInt); //first:rPos, second: gPos rank

	for (i = 0; i < num;) // identify all repetitive rPos
	{
		j = i + 1; while (j < num && vec[j].first == vec[i].first) j++;
		if (j - i > 1)
		{
			bTandemRepeat = true;
			//printf("Tandem repeat found: rPos = %d\n", vec[i].first);
			for (k = i; k < j; k++) SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
		}
		i = j;
	}
	//vec.clear(); vector<pair<int, int> >().swap(vec);
	if (bTandemRepeat) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("After tandem repeat checking\n"), ShowSeedInfo(SeedVec);
}

int IdentifyTranslocationRange(int i, int num, vector<pair<int, int> >& vec, vector<SeedPair_t>& SeedVec)
{
	int j, max_idx = vec[i].second;

	for (j = i + 1; j <= max_idx; j++)
	{
		if (vec[j].second > max_idx) max_idx = vec[j].second;
	}
	return max_idx;
}

void RemoveTranslocatedSeeds(vector<SeedPair_t>& SeedVec)
{
	int i, j, k, s1, s2, num;
	bool bTranslocation = false;
	vector<pair<int, int> > vec;

	num = (int)SeedVec.size(); if (num < 2) return;
	vec.resize(num);
	for (i = 0; i < num; i++) vec[i] = make_pair(SeedVec[i].rPos, i);
	sort(vec.begin(), vec.end(), CompByFirstInt); //first:rPos, second: gPos rank

	// checking translocation
	for (i = 0; i < num; i++)
	{
		if (vec[i].first != SeedVec[i].rPos)
		{
			bTranslocation = true;
			j = IdentifyTranslocationRange(i, num, vec, SeedVec);
			s1 = 0; s2 = 0;
			for (k = i; k <= j; k++)
			{
				if (k < vec[k].second) s1 += SeedVec[vec[k].second].rLen;
				else s2 += SeedVec[vec[k].second].rLen;
			}
			if (s1 > s2)
			{
				for (k = i; k <= j; k++)
					if (k > vec[k].second)
					{
						//printf("set #%d seed to zero length\n", vec[k].second + 1);
						SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
					}
			}
			else
			{
				for (k = i; k <= j; k++)
					if (k < vec[k].second)
					{
						//printf("set #%d seed to zero length\n", vec[k].second + 1);
						SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
					}
			}
			i = j;
		}
	}
	//vec.clear(); vector<pair<int, int> >().swap(vec);
	if (bTranslocation) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("After translocation checking\n"), ShowSeedInfo(SeedVec);
}

bool CheckSeedOverlapping(SeedPair_t& p1, SeedPair_t& p2)
{
	int iOverlap;
	bool bMaster = true;

	//printf("[1]: p1: r[%d-%d]=%d g[%lld-%lld]=%d vs p2: r[%d-%d]=%d g[%lld-%lld]=%d\n", p1.rPos, p1.rPos + p1.rLen - 1, p1.rLen, p1.gPos, p1.gPos + p1.gLen - 1, p1.gLen, p2.rPos, p2.rPos + p2.rLen - 1, p2.rLen, p2.gPos, p2.gPos + p2.gLen - 1, p2.gLen);
	if ((iOverlap = p1.rPos + p1.rLen - p2.rPos) > 0)
	{
		if (p1.rLen < p2.rLen)
		{
			bMaster = false;
			if (p1.rLen > iOverlap)
			{
				p1.gLen = (p1.rLen -= iOverlap);
			}
			else p1.rLen = p1.gLen = 0;
		}
		else
		{
			if (p2.rLen > iOverlap)
			{
				p2.rPos += iOverlap; p2.gPos += iOverlap;
				p2.gLen = (p2.rLen -= iOverlap);
			}
			else p2.rLen = p2.gLen = 0;
		}
	}
	if ((p1.rLen > 0 && p2.rLen > 0) && (iOverlap = p1.gPos + p1.gLen - p2.gPos) > 0)
	{
		if (p1.gLen < p2.gLen)
		{
			bMaster = false;
			if (p1.rLen > iOverlap)
			{
				p1.gLen = (p1.rLen -= iOverlap);
			}
			else p1.rLen = p1.gLen = 0;
		}
		else
		{
			if (p2.rLen > iOverlap)
			{
				p2.rPos += iOverlap; p2.gPos += iOverlap;
				p2.gLen = (p2.rLen -= iOverlap);
			}
			else p2.rLen = p2.gLen = 0;
		}
	}
	//printf("[2]: p1: r[%d-%d]=%d g[%lld-%lld]=%d vs p2: r[%d-%d]=%d g[%lld-%lld]=%d\n\n", p1.rPos, p1.rPos + p1.rLen - 1, p1.rLen, p1.gPos, p1.gPos + p1.gLen - 1, p1.gLen, p2.rPos, p2.rPos + p2.rLen - 1, p2.rLen, p2.gPos, p2.gPos + p2.gLen - 1, p2.gLen);
	return bMaster;
}

int LocateThePreviousSeedIdx(int i, vector<SeedPair_t>& SeedVec)
{
	while (i > 0 && SeedVec[i].rLen == 0) i--;
	if (i < 0) return 0;
	else return i;
}

void CheckOverlappingSeeds(vector<SeedPair_t>& SeedVec)
{
	int64_t gEnd;
	int i, j, num, rEnd;
	bool bNullSeed = false;

	num = (int)SeedVec.size(); if (num < 2) return;
	for (i = 0; i < num;)
	{
		if (SeedVec[i].rLen > 0)
		{
			rEnd = SeedVec[i].rPos + SeedVec[i].rLen - 1;
			gEnd = SeedVec[i].gPos + SeedVec[i].gLen - 1;

			//printf("overlap check for seed#%d, rEnd=%d, gEned=%lld\n", i + 1, rEnd, gEnd);
			for (j = i + 1; j < num; j++)
			{
				if (SeedVec[j].rLen == 0) continue;
				if (rEnd < SeedVec[j].rPos && gEnd < SeedVec[j].gPos) break;
				//printf("\ttest seed#%d\n", j + 1);
				if (CheckSeedOverlapping(SeedVec[i], SeedVec[j]) == false) break;
			}
			if (SeedVec[i].rLen == 0)
			{
				bNullSeed = true;
				i = LocateThePreviousSeedIdx(i - 1, SeedVec);
			}
			else i++;
		}
		else
		{
			bNullSeed = true; i++;
		}
	}
	if (bNullSeed) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("after overlap checking\n"), ShowSeedInfo(SeedVec);
}

void IdentifyNormalPairs(int rlen, int glen, vector<SeedPair_t>& SeedVec)
{
	SeedPair_t SeedPair;
	int i, j, rGaps, gGaps, num;

	if (SeedVec.size() > 1)
	{
		RemoveTandemRepeatSeeds(SeedVec);
		RemoveTranslocatedSeeds(SeedVec);
		CheckOverlappingSeeds(SeedVec);

		SeedPair.bSimple = false;
		num = (int)SeedVec.size();

		for (i = 0, j = 1; j < num; i++, j++)
		{
			rGaps = SeedVec[j].rPos - (SeedVec[i].rPos + SeedVec[i].rLen); if (rGaps < 0) rGaps = 0;
			gGaps = SeedVec[j].gPos - (SeedVec[i].gPos + SeedVec[i].gLen); if (gGaps < 0) gGaps = 0;
			//printf("check %d and %d: rGaps=%d, gGaps=%d\n", i + 1, j + 1, rGaps, gGaps); fflush(stdout);
			if (rGaps > 0 || gGaps > 0)
			{
				SeedPair.rPos = SeedVec[i].rPos + SeedVec[i].rLen;
				SeedPair.gPos = SeedVec[i].gPos + SeedVec[i].gLen;
				SeedPair.PosDiff = SeedPair.gPos - SeedPair.rPos;
				SeedPair.rLen = rGaps; SeedPair.gLen = gGaps;
				//printf("Add normal pair:\nR1[%d-%d] G1[%lld-%lld]\nNR[%d-%d]=%d NG[%ld-%ld]=%d\nR2[%d-%d] G2[%lld-%lld]\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen, SeedVec[j].rPos, SeedVec[j].rPos + SeedVec[j].rLen - 1, SeedVec[j].gPos, SeedVec[j].gPos + SeedVec[j].gLen - 1); fflush(stdout);
				SeedVec.push_back(SeedPair);
			}
		}
		if ((int)SeedVec.size() > num) inplace_merge(SeedVec.begin(), SeedVec.begin()+num, SeedVec.end(), CompByGenomePos);
		//if(bDebugMode) printf("After filling gaps\n"), ShowSeedInfo(SeedVec);
	}
	// Check missing blocks at both ends
	if (SeedVec.size() > 0)
	{
		i = 0;
		//printf("Identify heading pair R[%d-%d]=%d, G[%lld-%lld]=%d\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].rLen, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedVec[i].gLen); fflush(stdout);
		rGaps = SeedVec[i].rPos > 0 ? SeedVec[i].rPos : 0;
		gGaps = glen > 0 ? SeedVec[i].gPos : rGaps;
		//printf("rGaps=%d, gGaps=%d\n", rGaps, gGaps); fflush(stdout);
		if (rGaps > 0 || gGaps > 0)
		{
			SeedPair.rPos = 0;
			SeedPair.gPos = SeedVec[i].gPos - gGaps;
			if (SeedPair.gPos < 0) { SeedPair.gPos = 0; gGaps += SeedPair.gPos; }
			SeedPair.PosDiff = SeedPair.gPos;
			SeedPair.bSimple = false;
			SeedPair.rLen = rGaps;
			SeedPair.gLen = gGaps;
			SeedVec.insert(SeedVec.begin(), SeedPair);

			//if (bDebugMode) printf("Add missing head: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
		}
		i = (int)SeedVec.size() - 1;
		rGaps = rlen - (SeedVec[i].rPos + SeedVec[i].rLen);
		gGaps = glen > 0 ? glen - (SeedVec[i].gPos + SeedVec[i].gLen) : rGaps;

		if (rGaps > 0 || gGaps > 0)
		{
			SeedPair.bSimple = false;
			SeedPair.rPos = SeedVec[i].rPos + SeedVec[i].rLen;
			SeedPair.gPos = SeedVec[i].gPos + SeedVec[i].gLen;
			SeedPair.rLen = rGaps;
			SeedPair.gLen = gGaps;
			SeedVec.push_back(SeedPair);

			//if (bDebugMode) printf("Add missing tail: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
		}
	}
	//if (bDebugMode) printf("after generating normal pairs\n"), ShowSeedInfo(SeedVec);
}

string GenerateCIGAR(vector<pair<int, char> >& cigar_vec)
{
	int i, num, c;
	char state, buf[8];
	string cigar_str;

	for (state = '\0', num = (int)cigar_vec.size(), c = 0, i = 0; i != num; i++)
	{
		if (cigar_vec[i].second != state)
		{
			if (c > 0) sprintf(buf, "%d%c", c, state), cigar_str += buf;

			c = cigar_vec[i].first; state = cigar_vec[i].second;
		}
		else c += cigar_vec[i].first;
	}
	if (c > 0) sprintf(buf, "%d%c", c, state), cigar_str += buf;

	if (bDebugMode) printf("CIGAR=%s\n\n\n", cigar_str.c_str());

	return cigar_str;
}

Coordinate_t GenCoordinateInfo(bool bFirstRead, int64_t gPos, int64_t end_gPos, vector<pair<int, char> >& cigar_vec)
{
	Coordinate_t coor;
	map<int64_t, int>::iterator iter;

	if (gPos < GenomeSize) // forward strand
	{
		if (bFirstRead) coor.bDir = true;
		else coor.bDir = false;

		if (iChromsomeNum == 1)
		{
			coor.ChromosomeIdx = 0;
			coor.gPos = gPos + 1;
		}
		else
		{
			iter = ChrLocMap.lower_bound(gPos);
			coor.ChromosomeIdx = iter->second;
			coor.gPos = gPos + 1 - ChromosomeVec[coor.ChromosomeIdx].FowardLocation;
			//if (bDebugMode) printf("matched chr=%s, loc=%lld, gPos: %lld -> %lld\n", ChromosomeVec[coor.ChromosomeIdx].name, ChromosomeVec[coor.ChromosomeIdx].FowardLocation, gPos, coor.gPos);
		}
	}
	else
	{
		if (bFirstRead) coor.bDir = false;
		else coor.bDir = true;

		reverse(cigar_vec.begin(), cigar_vec.end());

		if (iChromsomeNum == 1)
		{
			coor.ChromosomeIdx = 0;
			coor.gPos = TwoGenomeSize - end_gPos;
		}
		else
		{
			iter = ChrLocMap.lower_bound(gPos);
			coor.gPos = iter->first - end_gPos + 1; coor.ChromosomeIdx = iter->second;
			//if(bDebugMode) printf("matched chr=%s, loc=%lld, gPos: %lld -> %lld\n", ChromosomeVec[coor.ChromosomeIdx].name, ChromosomeVec[coor.ChromosomeIdx].ReverseLocation, gPos, coor.gPos);
		}
	}
	//if (bDebugMode) printf("gPos: %lld --> %lld %s\n", gPos, coor.gPos, (coor.bDir? "Forward":"Reverse"));

	coor.CIGAR = GenerateCIGAR(cigar_vec);

	return coor;
}

string ReverseCIGAR(string& CIGAR)
{
	string RevCIGAR;
	int len, pos1, pos2;

	len = (int)CIGAR.length();

	for (pos1 = 0, pos2 = 1; pos1<len; pos2++)
	{
		if (isalpha(CIGAR[pos2]))
		{
			RevCIGAR.insert(0, CIGAR.substr(pos1, pos2 - pos1 + 1));
			pos1 = pos2 + 1;
		}
	}
	return RevCIGAR;
}

bool CheckCoordinateValidity(vector<SeedPair_t>& SeedVec)
{
	bool bValid = true;
	int64_t gPos1=0, gPos2=TwoGenomeSize;

	for (vector<SeedPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end(); iter++)
	{
		if (iter->gLen > 0)
		{
			gPos1 = iter->gPos;
			break;
		}
	}
	for (vector<SeedPair_t>::reverse_iterator iter = SeedVec.rbegin(); iter != SeedVec.rend(); iter++)
	{
		if (iter->gLen > 0)
		{
			gPos2 = iter->gPos + iter->gLen - 1;
			break;
		}
	}	
	if ((gPos1 < GenomeSize && gPos2 >= GenomeSize) || (gPos1 >= GenomeSize && gPos2 < GenomeSize) || ChrLocMap.lower_bound(gPos1)->second != ChrLocMap.lower_bound(gPos2)->second)
	{
		bValid = false;
		//if (bDebugMode) fprintf(stderr, "%lld and %lld are not in the same chromosome!\n", gPos1, gPos2);
	}
	return bValid;
}

int GapPenalty(vector<pair<int, char> >& cigar_vec)
{
	int GP = 0;
	for (vector<pair<int, char> >::iterator iter = cigar_vec.begin(); iter != cigar_vec.end(); iter++)
	{
		if (iter->second == 'I' || iter->second == 'D') GP += iter->first;
	}
	//if (bDebugMode) printf("GapPenalty=%d\n", GP), fflush(stdout);

	return GP;
}

void GenMappingReport(bool bFirstRead, ReadItem_t& read, vector<AlignmentCandidate_t>& AlignmentVec)
{
	int i, j, num;
	map<int, int>::iterator iter;
	vector<pair<int, char> > cigar_vec;
	map<int64_t, int>::iterator ChrIter;

	//if (bDebugMode) printf("\n\n%s\nGenerate alignment for read %s (%d cans)\n", string().assign(100, '=').c_str(), read.header, (int)AlignmentVec.size()), fflush(stdout);
	
	read.score = read.iBestAlnCanIdx = 0;
	if ((read.CanNum = (int)AlignmentVec.size()) > 0)
	{
		read.AlnReportArr = new AlignmentReport_t[read.CanNum];
		for (i = 0; i != read.CanNum; i++)
		{
			read.AlnReportArr[i].AlnScore = 0;
			read.AlnReportArr[i].PairedAlnCanIdx = AlignmentVec[i].PairedAlnCanIdx;

			if (AlignmentVec[i].Score == 0 || (bPacBioData && read.score > 0)) continue;

			IdentifyNormalPairs(read.rlen, -1, AlignmentVec[i].SeedVec); // fill missing framgment pairs (normal pairs) between simple pairs
			//if (bDebugMode)
			//{
			//	printf("Process candidate#%d (Score = %d, SegmentPair#=%d): \n", i + 1, AlignmentVec[i].Score, (int)AlignmentVec[i].SeedVec.size());
			//	ShowSeedInfo(AlignmentVec[i].SeedVec);
			//}
			if (CheckCoordinateValidity(AlignmentVec[i].SeedVec) == false) continue;

			cigar_vec.clear();
			for (num = (int)AlignmentVec[i].SeedVec.size(), j = 0; j != num; j++)
			{
				if (AlignmentVec[i].SeedVec[j].rLen == 0 && AlignmentVec[i].SeedVec[j].gLen == 0) continue;
				else if (AlignmentVec[i].SeedVec[j].bSimple)
				{
					//if (bDebugMode) ShowFragmentPair(AlignmentVec[i].SeedVec[j]);
					cigar_vec.push_back(make_pair(AlignmentVec[i].SeedVec[j].rLen, 'M'));
					read.AlnReportArr[i].AlnScore += AlignmentVec[i].SeedVec[j].rLen;
				}
				else
				{
					//if (bDebugMode) printf("Check normal pair#%d: R[%d-%d]=%d G[%lld-%lld]=%d\n", j + 1, AlignmentVec[i].SeedVec[j].rPos, AlignmentVec[i].SeedVec[j].rPos + AlignmentVec[i].SeedVec[j].rLen - 1, AlignmentVec[i].SeedVec[j].rLen, AlignmentVec[i].SeedVec[j].gPos, AlignmentVec[i].SeedVec[j].gPos + AlignmentVec[i].SeedVec[j].gLen - 1, AlignmentVec[i].SeedVec[j].gLen);
					if (j == 0)
					{
						if (AlignmentVec[i].SeedVec[0].rLen > 5000)
						{
							cigar_vec.push_back(make_pair(AlignmentVec[i].SeedVec[0].rLen, 'S'));
							AlignmentVec[i].SeedVec[0].gPos = AlignmentVec[i].SeedVec[1].gPos;
							AlignmentVec[i].SeedVec[0].gLen = 0;
						}
						else read.AlnReportArr[i].AlnScore += ProcessHeadSequencePair(read.seq, AlignmentVec[i].SeedVec[0], cigar_vec);
					}
					else if (j == num - 1)
					{
						if (AlignmentVec[i].SeedVec[j].rLen > 5000)
						{
							cigar_vec.push_back(make_pair(AlignmentVec[i].SeedVec[j].rLen, 'S'));
							AlignmentVec[i].SeedVec[j].gPos = AlignmentVec[i].SeedVec[j-1].gPos + AlignmentVec[i].SeedVec[j - 1].gLen;
							AlignmentVec[i].SeedVec[j].gLen = 0;
						}
						else read.AlnReportArr[i].AlnScore += ProcessTailSequencePair(read.seq, AlignmentVec[i].SeedVec[j], cigar_vec);
					}
					else
					{
						read.AlnReportArr[i].AlnScore += ProcessNormalSequencePair(read.seq, AlignmentVec[i].SeedVec[j], cigar_vec);
					}
				}
			}
			//if (bDebugMode) printf("Alignment score = %d (rlen=%d) \n", read.AlnReportArr[i].AlnScore, read.rlen), fflush(stdout);

			if (!bPacBioData && cigar_vec.size() > 1)
			{
				read.AlnReportArr[i].AlnScore -= GapPenalty(cigar_vec);
				if (read.AlnReportArr[i].AlnScore <= 0)
				{
					read.AlnReportArr[i].AlnScore = 0;
					continue;
				}
			}
			if((read.AlnReportArr[i].coor = GenCoordinateInfo(bFirstRead, AlignmentVec[i].SeedVec[0].gPos, (AlignmentVec[i].SeedVec[num - 1].gPos + AlignmentVec[i].SeedVec[num - 1].gLen - 1), cigar_vec)).gPos <= 0) read.AlnReportArr[i].AlnScore = 0;

			if (read.AlnReportArr[i].AlnScore > read.score)
			{
				read.iBestAlnCanIdx = i;
				read.sub_score = read.score;
				read.score = read.AlnReportArr[i].AlnScore;
			}
			else if (read.AlnReportArr[i].AlnScore == read.score)
			{
				read.sub_score = read.score;
				if (!bMultiHit && ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].len > ChromosomeVec[read.AlnReportArr[read.iBestAlnCanIdx].coor.ChromosomeIdx].len) read.iBestAlnCanIdx = i;
			}
		}
	}
	else
	{
		read.CanNum = 1; read.iBestAlnCanIdx = 0;
		read.AlnReportArr = new AlignmentReport_t[1];
		read.AlnReportArr[0].AlnScore = 0;
		read.AlnReportArr[0].PairedAlnCanIdx = -1;
	}
}
