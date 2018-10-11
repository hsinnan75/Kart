#include <sys/stat.h>
#include "structure.h"

bwt_t *Refbwt;
bwaidx_t *RefIdx;
const char* VersionStr = "2.4.6";
vector<string> ReadFileNameVec1, ReadFileNameVec2;
char *RefSequence, *IndexFileName, *OutputFileName;
int iThreadNum, MaxInsertSize, MaxGaps, MinSeedLength, OutputFileFormat;
bool bDebugMode, bPairEnd, bPacBioData, bMultiHit, gzCompressed, FastQFormat;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "kart v%s (Hsin-Nan Lin & Wen-Lian Hsu)\n\n", VersionStr);
	fprintf(stderr, "Usage: %s -i Index_Prefix -f <ReadFile_A1 ReadFile_B1 ...> [-f2 <ReadFile_A2 ReadFile_B2 ...>] > out.sam\n\n", program);
	fprintf(stderr, "Options: -t INT        number of threads [4]\n");
	fprintf(stderr, "         -f            files with #1 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stderr, "         -f2           files with #2 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stderr, "         -o            output filename [stdout]\n");
	fprintf(stderr, "         -m            output multiple alignments\n");
	fprintf(stderr, "         -g INT        max gaps (indels) [5]\n");
	fprintf(stderr, "         -p            paired-end reads are interlaced in the same file\n");
	fprintf(stderr, "         -pacbio       pacbio data\n");
	fprintf(stderr, "         -v			version\n");
	fprintf(stderr, "\n");
}

bool CheckOutputFileName()
{
	bool bRet = true;

	if (OutputFileName != NULL)
	{
		struct stat s;
		string filename, FileExt;

		filename = OutputFileName; FileExt = filename.substr(filename.find_last_of('.') + 1);
		if (FileExt == "gz") OutputFileFormat = 1;
		if (stat(OutputFileName, &s) == 0)
		{
			if (s.st_mode & S_IFDIR)
			{
				bRet = false;
				fprintf(stderr, "Warning: %s is a directory!\n", OutputFileName);
			}
			else if (s.st_mode & S_IFREG)
			{
			}
			else
			{
				bRet = false;
				fprintf(stderr, "Warning: %s is not a regular file!\n", OutputFileName);
			}
		}
	}
	return bRet;
}

bool CheckInputFiles()
{
	struct stat s;
	bool bRet = true;

	for (vector<string>::iterator iter = ReadFileNameVec1.begin(); iter != ReadFileNameVec1.end(); iter++)
	{
		if (stat(iter->c_str(), &s) == -1)
		{
			bRet = false;
			fprintf(stderr, "Cannot access file:[%s]\n", (char*)iter->c_str());
		}
	}
	for (vector<string>::iterator iter = ReadFileNameVec2.begin(); iter != ReadFileNameVec2.end(); iter++)
	{
		if (stat(iter->c_str(), &s) == -1)
		{
			bRet = false;
			fprintf(stderr, "Cannot access file:[%s]\n", (char*)iter->c_str());
		}
	}
	return bRet;
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	MaxGaps = 5;
	iThreadNum = 4;
	bPairEnd = false;
	bDebugMode = false;
	MaxInsertSize = 1500;
	MinSeedLength = 0;
	bPacBioData = false;
	bMultiHit = false;
	FastQFormat = true;
	OutputFileFormat = 0; // 0:sam 1:sam.gz
	RefSequence = IndexFileName = OutputFileName = NULL;

	if (argc == 1 || strcmp(argv[1], "-h") == 0) ShowProgramUsage(argv[0]);
	else if (strcmp(argv[1], "update") == 0)
	{
		system("git fetch; git merge origin/master master;make");
		exit(0);
	}
	else
	{
		for (i = 1; i < argc; i++)
		{
			parameter = argv[i];

			if (parameter == "-i") IndexFileName = argv[++i];
			else if (parameter == "-f")
			{
				while (++i < argc && argv[i][0] != '-') ReadFileNameVec1.push_back(argv[i]);
				i--;
			}
			else if (parameter == "-f2")
			{
				while (++i < argc && argv[i][0] != '-') ReadFileNameVec2.push_back(argv[i]);
				i--;
			}
			else if (parameter == "-t")
			{
				if ((iThreadNum = atoi(argv[++i])) > 40)
				{
					fprintf(stderr, "Warning! Thread number is limited to 40!\n");
					iThreadNum = 40;
				}
			}
			else if (parameter == "-g")
			{
				if ((MaxGaps = atoi(argv[++i])) < 0) MaxGaps = 0;
			}
			else if (parameter == "-o") OutputFileName = argv[++i];
			else if (parameter == "-pacbio") bPacBioData = true;
			else if (parameter == "-m") bMultiHit = true;
			else if (parameter == "-pair" || parameter == "-p") bPairEnd = true;
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else if (parameter == "-v" || parameter == "--version")
			{
				fprintf(stderr, "kart v%s\n\n", VersionStr);
				exit(0);
			}
			else
			{
				fprintf(stderr, "Error! Unknow parameter: %s\n", argv[i]);
				ShowProgramUsage(argv[0]);
				exit(1);
			}
		}

		if (ReadFileNameVec1.size() == 0)
		{
			fprintf(stderr, "Error! Please specify a valid read input!\n");
			ShowProgramUsage(argv[0]);
			exit(1);
		}
		if (ReadFileNameVec2.size() > 0 && ReadFileNameVec1.size() != ReadFileNameVec2.size())
		{
			fprintf(stderr, "Error! Paired-end reads input numbers do not match!\n");
			fprintf(stderr, "Read1:\n"); for (vector<string>::iterator iter = ReadFileNameVec1.begin(); iter != ReadFileNameVec1.end(); iter++) fprintf(stderr, "\t%s\n", (char*)iter->c_str());
			fprintf(stderr, "Read2:\n"); for (vector<string>::iterator iter = ReadFileNameVec2.begin(); iter != ReadFileNameVec2.end(); iter++) fprintf(stderr, "\t%s\n", (char*)iter->c_str());
			exit(1);
		}
		if (CheckInputFiles() == false || CheckOutputFileName() == false) exit(0);
		if (IndexFileName != NULL && CheckBWAIndexFiles(IndexFileName)) RefIdx = bwa_idx_load(IndexFileName);
		else
		{
			fprintf(stderr, "Error! Please specify a valid reference index!\n");
			ShowProgramUsage(argv[0]);
			exit(1);
		}
		if (RefIdx == 0)
		{
			fprintf(stderr, "\n\nError! Index files are corrupt!\n");
			exit(1);
		}
		else
		{
			Refbwt = RefIdx->bwt;
			RestoreReferenceInfo();
			Mapping();
			bwa_idx_destroy(RefIdx);
			if (RefSequence != NULL) delete[] RefSequence;
		}
	}
	return 0;
}
