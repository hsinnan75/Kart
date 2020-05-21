#include <sys/stat.h>
#include "structure.h"
extern "C"
{
	int bwa_idx_build(const char *fa, const char *prefix);
}

bwt_t *Refbwt;
bwaidx_t *RefIdx;
const char* VersionStr = "2.5.6";
vector<string> ReadFileNameVec1, ReadFileNameVec2;
char *RefSequence, *IndexFileName, *OutputFileName;
int iThreadNum, MaxInsertSize, MaxGaps, MinSeedLength, OutputFileFormat;
bool bDebugMode, bPairEnd, bPacBioData, bMultiHit, gzCompressed, FastQFormat, bSilent;

void ShowProgramUsage(const char* program)
{
	fprintf(stdout, "kart v%s (Hsin-Nan Lin & Wen-Lian Hsu)\n\n", VersionStr);
	fprintf(stdout, "Usage: %s -i Index_Prefix -f <ReadFile_A1 ReadFile_B1 ...> [-f2 <ReadFile_A2 ReadFile_B2 ...>] -o Output\n\n", program);
	fprintf(stdout, "Options: -t INT        number of threads [4]\n");
	fprintf(stdout, "         -f            files with #1 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stdout, "         -f2           files with #2 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stdout, "         -o            alignment filename in SAM format [output.sam]\n");
	fprintf(stdout, "         -bo           alignment filename in BAM format\n");
	fprintf(stdout, "         -m            output multiple alignments\n");
	fprintf(stdout, "         -g INT        max gaps (indels) [5]\n");
	fprintf(stdout, "         -p            paired-end reads are interlaced in the same file\n");
	fprintf(stdout, "         -pacbio       pacbio data\n");
	fprintf(stdout, "         -v            version\n");
	fprintf(stdout, "\n");
}

bool CheckOutputFileName()
{
	bool bRet = true;

	if (strcmp(OutputFileName, "output.sam")!=0)
	{
		struct stat s;
		if (stat(OutputFileName, &s) == 0)
		{
			if (s.st_mode & S_IFDIR)
			{
				bRet = false;
				fprintf(stdout, "Warning: %s is a directory!\n", OutputFileName);
			}
		}
		int i, len = strlen(OutputFileName);
		for (i = 0; i < len; i++)
		{
			if (isalnum(OutputFileName[i]) || OutputFileName[i] == '/' || OutputFileName[i] == '.' || OutputFileName[i] == '-');
			else
			{
				bRet = false;
				fprintf(stdout, "Warning: [%s] is not a valid filename!\n", OutputFileName);
				break;
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
			fprintf(stdout, "Cannot access file:[%s]\n", (char*)iter->c_str());
		}
	}
	for (vector<string>::iterator iter = ReadFileNameVec2.begin(); iter != ReadFileNameVec2.end(); iter++)
	{
		if (stat(iter->c_str(), &s) == -1)
		{
			bRet = false;
			fprintf(stdout, "Cannot access file:[%s]\n", (char*)iter->c_str());
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
	bSilent = false;
	FastQFormat = true;
	OutputFileName = (char*)"output.sam";
	OutputFileFormat = 0; // 0:sam 1:bam
	RefSequence = IndexFileName = NULL;

	if (argc == 1 || strcmp(argv[1], "-h") == 0) ShowProgramUsage(argv[0]);
	else if (strcmp(argv[1], "update") == 0)
	{
		i = system("git fetch; git merge origin/master master;make");
		exit(0);
	}
	else if (strcmp(argv[1], "index") == 0)
	{
		if (argc == 4) bwa_idx_build(argv[2], argv[3]);
		else
		{
			fprintf(stderr, "usage: %s index ref.fa prefix\n", argv[0]);
		}
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
			else if (parameter == "-t" && i + 1 < argc)
			{
				if ((iThreadNum = atoi(argv[++i])) <= 0)
				{
					fprintf(stdout, "Warning! Thread number should be a positive number!\n");
					iThreadNum = 4;
				}
			}
			else if (parameter == "-g")
			{
				if ((MaxGaps = atoi(argv[++i])) < 0) MaxGaps = 0;
			}
			else if (parameter == "-o")
			{
				OutputFileFormat = 0;
				OutputFileName = argv[++i];
			}
			else if (parameter == "-bo")
			{
				OutputFileFormat = 1;
				OutputFileName = argv[++i];
			}
			else if (parameter == "-silent") bSilent = true;
			else if (parameter == "-pacbio") bPacBioData = true;
			else if (parameter == "-m") bMultiHit = true;
			else if (parameter == "-pair" || parameter == "-p") bPairEnd = true;
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else if (parameter == "-v" || parameter == "--version")
			{
				fprintf(stdout, "kart v%s\n\n", VersionStr);
				exit(0);
			}
			else
			{
				fprintf(stdout, "Error! Unknown parameter: %s\n", argv[i]);
				ShowProgramUsage(argv[0]);
				exit(1);
			}
		}

		if (ReadFileNameVec1.size() == 0)
		{
			fprintf(stdout, "Error! Please specify a valid read input!\n");
			ShowProgramUsage(argv[0]);
			exit(1);
		}
		if (ReadFileNameVec2.size() > 0 && ReadFileNameVec1.size() != ReadFileNameVec2.size())
		{
			fprintf(stdout, "Error! Paired-end reads input numbers do not match!\n");
			fprintf(stdout, "Read1:\n"); for (vector<string>::iterator iter = ReadFileNameVec1.begin(); iter != ReadFileNameVec1.end(); iter++) fprintf(stdout, "\t%s\n", (char*)iter->c_str());
			fprintf(stdout, "Read2:\n"); for (vector<string>::iterator iter = ReadFileNameVec2.begin(); iter != ReadFileNameVec2.end(); iter++) fprintf(stdout, "\t%s\n", (char*)iter->c_str());
			exit(1);
		}
		if (CheckInputFiles() == false || CheckOutputFileName() == false) exit(0);
		if (IndexFileName != NULL && CheckBWAIndexFiles(IndexFileName)) RefIdx = bwa_idx_load(IndexFileName);
		else
		{
			fprintf(stdout, "Error! Please specify a valid reference index!\n");
			ShowProgramUsage(argv[0]);
			exit(1);
		}
		if (RefIdx == 0)
		{
			fprintf(stdout, "\n\nError! Index files are corrupt!\n");
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
