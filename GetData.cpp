#include "structure.h"

int iChromsomeNum;
map<int64_t, int> ChrLocMap;
int64_t GenomeSize, TwoGenomeSize;
vector<Chromosome_t> ChromosomeVec;

int IdentifyHeaderBoundary(string& str)
{
	int i = 0;
	for (string::iterator iter = str.begin(); iter != str.end(); iter++)
	{
		if (*iter == ' ' || *iter == '\t')
		{
			i = iter - str.begin();
			break;
		}
	}
	if (i == 0) i = (int)str.length();

	return i-1;
}

bool GetNextEntry(fstream& file, string& header, string& seq)
{
	char ch;
	string str;

	getline(file, header);
	if (header == "" || file.eof()) return false;
	else
	{
		if (FastQFormat)
		{
			getline(file, seq); getline(file, str), getline(file, str);
		}
		else
		{
			seq.clear();
			while (!file.eof())
			{
				getline(file, str); seq += str; file.get(ch);
				if (file.eof()) break;
				else file.unget();
				if (ch == '>') break;
			}
		}
		return true;
	}
}

int GetNextChunk(bool bSepLibrary, fstream& file, fstream& file2, ReadItem_t* ReadArr)
{
	int i, p, len, iCount = 0;
	string header, seq, quality;

	while (!file.eof())
	{
		if (!GetNextEntry(file, header, seq)) break;
		//header
		p = IdentifyHeaderBoundary(header);
		ReadArr[iCount].header = new char[p + 1]; ReadArr[iCount].header[p] = '\0';
		strncpy(ReadArr[iCount].header, (char*)header.c_str() + 1, p);

		//sequence
		ReadArr[iCount].rlen = len = (int)seq.length();
		ReadArr[iCount].seq = new char[len + 1];
		strcpy(ReadArr[iCount].seq, seq.c_str());

		ReadArr[iCount].EncodeSeq = new uint8_t[len];
		for (i = 0; i != len; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];

		iCount++;

		if (bSepLibrary) GetNextEntry(file2, header, seq);
		else if (!GetNextEntry(file, header, seq)) break;
		//header
		p = IdentifyHeaderBoundary(header);
		ReadArr[iCount].header = new char[p + 1]; ReadArr[iCount].header[p] = '\0';
		strncpy(ReadArr[iCount].header, (char*)header.c_str() + 1, p);
		//sequence
		ReadArr[iCount].rlen = len = (int)seq.length();
		ReadArr[iCount].seq = new char[len + 1];

		if (bPairEnd) GetComplementarySeq(len, (char*)seq.c_str(), ReadArr[iCount].seq);
		else strcpy(ReadArr[iCount].seq, seq.c_str());

		ReadArr[iCount].EncodeSeq = new uint8_t[len];
		for (i = 0; i != len; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];

		iCount++;

		if ((bPacBioData && iCount == 8) || iCount == ReadChunkSize) break;
	}
	return iCount;
}

bool CheckBWAIndexFiles(string IndexPrefix)
{
	fstream file;
	string filename;
	bool bChecked=true;

	filename = IndexPrefix + ".ann"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();
	
	filename = IndexPrefix + ".amb"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();

	filename = IndexPrefix + ".pac"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();

	return bChecked;
}
