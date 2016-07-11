#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>

#define iShift 30

using namespace std;

bool bShowWrongCase=false;
bool bShowUnmapped=true;
bool bShowBadAlign=false;

int TotalQuery=0;
map<string, int> QueryMap;

bool CheckPosConsistency(int Length, int TrueLocation, int PredictedLocation)
{
	if(TrueLocation >= PredictedLocation && TrueLocation - PredictedLocation < (iShift+ Length)) return true;
	else if (TrueLocation < PredictedLocation && PredictedLocation - TrueLocation < (iShift + Length)) return true;
	else return false;
}

void Evaluation(string SamFileName)
{
	fstream file;
	stringstream ss;
	bool bCorLocation, bMapped;
	unsigned int pos, TrueLocation;
	long long iTotalBase, iTotalCorBase;
	int iFlag, MAPQ, iReadNum=0, iUnmapped=0, iCorLocation=0, iBadMAPQ=0;
	string str, Prevheader, Readheader, PosInfoStr, ChName, tmp, seq, region, CIGAR;

	iTotalBase = iTotalCorBase = 0;

	pos = 0; Prevheader = "";
	bMapped=bCorLocation=false;
	file.open(SamFileName.c_str(), ios_base::in);

	if(file.is_open())
	{
		while(!file.eof())
		{
			getline(file, str); if(str=="") break;
			if(str[0]=='@') continue;

			ss.clear(); ss.str(str);
			ss >> Readheader >> iFlag >> ChName >> PosInfoStr >> MAPQ >> tmp >> tmp >> tmp >> tmp >> seq;

			if (Readheader == Prevheader) continue;
			Prevheader = Readheader; iReadNum++;

			bMapped=bCorLocation=false;
			if(ChName!="*" && MAPQ==0) iBadMAPQ++;

			if(PosInfoStr!="*") pos = (unsigned int)atoi(PosInfoStr.c_str()); pos-=1;
			if(Readheader.find_first_of('=')!=string::npos) PosInfoStr = Readheader.substr(Readheader.find_first_of('=')+1);
			else PosInfoStr = Readheader.substr(0, Readheader.find_first_of('_'));
			TrueLocation = (unsigned int)atoi(PosInfoStr.c_str());

			if(ChName=="*")
			{
				iUnmapped++;
				if(bShowUnmapped && bShowWrongCase) cout << Readheader << endl;
			}
			else
			{
				if(bCorLocation==false && CheckPosConsistency((int)seq.length(), TrueLocation, pos))
					bCorLocation=true;
				bMapped=true;

				if(bCorLocation==false)
				{
					if(MAPQ > 0 && bShowWrongCase) cout << Readheader << endl;
				}
				else iCorLocation++;
			}
		}
		file.close();
	}
	if(TotalQuery==0) TotalQuery=iReadNum;
	else if(iReadNum > (int)(TotalQuery*0.995)) TotalQuery=iReadNum;

	cerr << endl << endl << "filename=" << SamFileName << endl;
	cerr << "# of reads= " << TotalQuery << endl;
	if(iReadNum>0) cerr << "# of mapped reads= " << (iReadNum-iUnmapped) << " (" << (int)(10000*(1.0*(iReadNum-iUnmapped)/TotalQuery)+0.5)/100.0 << "%)" << endl;
	if(iReadNum>0) cerr << "# of mapq_0=" << iBadMAPQ << " (" << (int)(10000*(1.0*iBadMAPQ/iReadNum)+0.5)/100.0 << "%)" << endl;
	if(iReadNum>0) cerr << "precision= " << iCorLocation << " (" << (int)(10000*(1.0*iCorLocation/(iReadNum-iUnmapped))+0.5)/100.0 << "%)" << endl;
	if(iReadNum>0) cerr << "recall= " << iCorLocation << " (" << (int)(10000*(1.0*iCorLocation/TotalQuery)+0.5)/100.0 << "%)" << endl;
	cerr << endl << endl;

	//// precision, recall
	//FILE *rfile = fopen("exp.csv", "a");
	//if (iReadNum > 0) fprintf(rfile, "%s: %.1f\t%.1f\n", SamFileName.c_str(), (int)(1000 * (1.0*iCorLocation / (iReadNum - iUnmapped)) + 0.5) / 10.0, (int)(1000 * (1.0*iCorLocation / TotalQuery) + 0.5) / 10.0);
	//else fprintf(rfile, "%s: 0\t0\n", SamFileName.c_str());

	//fclose(rfile);
}

int main(int argc, char* argv[])
{
	string str;

	if (argc == 1)
	{
		printf("usage: %s SamFile\n\n", argv[0]);
	}
	else
	{
		for (int i = 2; i < argc; i++)
		{
			if ((str = argv[i]) == "-d") bShowWrongCase = true;
		}
		Evaluation(argv[1]);
	}
	return 0;
}
