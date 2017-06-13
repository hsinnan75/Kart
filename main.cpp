#include <sys/stat.h>
#include "structure.h"

bwt_t *Refbwt;
bwaidx_t *RefIdx;
const char* VersionStr = "2.2.0";
int iThreadNum, MaxInsertSize, MaxGaps, MinSeedLength, OutputFileFormat;
bool bDebugMode, bPairEnd, bPacBioData, bMultiHit, gzCompressed, FastQFormat;
char *RefSequence, *IndexFileName, *ReadFileName, *ReadFileName2, *OutputFileName;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "kart v%s\n", VersionStr);
	fprintf(stderr, "Usage: %s -i Index_Prefix -f ReadFile [-f2 ReadFile2] > out.sam\n\n", program);
	fprintf(stderr, "Options: -t INT        number of threads [16]\n");
	fprintf(stderr, "         -f            files with #1 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stderr, "         -f2           files with #2 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stderr, "         -o            output filename [stdout]\n");
	fprintf(stderr, "         -m            output multiple alignments\n");
	fprintf(stderr, "         -g INT        max gaps (indels) [5]\n");
	fprintf(stderr, "         -p            paired-end reads are interlaced in the same file\n");
	fprintf(stderr, "         -pacbio       pacbio data\n");
	fprintf(stderr, "\n");
}

bool CheckReadFile(char* filename, bool& bReadFormat)
{
	fstream file;
	string header;
	bool bCheck = true;

	file.open(filename, ios_base::in);
	if (!file.is_open()) return false;
	else
	{
		getline(file, header);
		if (header == "") return false;
		else
		{
			if (header[0] == '@') bReadFormat = true;
			else bReadFormat = false;
		}
	}
	file.close();

	return bCheck;
}

bool CheckCompressedFile(char* filename)
{
	gzFile file;
	bool bCheck = true;

	if ((file = gzopen(filename, "rb")) == Z_NULL) bCheck = false;
	gzclose(file);

	return bCheck;
}

bool Check_gzInputFormat()
{
	string filename, FileExt;
	bool bCompressed = false;

	filename = ReadFileName; FileExt = filename.substr(filename.find_last_of('.') + 1);
	if (FileExt == "gz") bCompressed = true;

	return bCompressed;
}

bool CheckOutputFileName()
{
	struct stat s;
	bool bRet = true;
	string filename, FileExt;

	filename = OutputFileName; FileExt = filename.substr(filename.find_last_of('.') + 1);

	if (FileExt == "gz") OutputFileFormat = 1;
	else if (FileExt == "bam") OutputFileFormat = 2;

	//if (OutputFileFormat == 0) fprintf(stderr, "OutputFile = %s [format=sam]\n", OutputFileName);
	//else if (OutputFileFormat == 1) fprintf(stderr, "OutputFile = %s [format=sam.gz]\n", OutputFileName);
	//else fprintf(stderr, "OutputFile = %s [format=bam]\n", OutputFileName);

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
	return bRet;
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	MaxGaps = 5;
	iThreadNum = 16;
	bPairEnd = false;
	bDebugMode = false;
	MaxInsertSize = 1500;
	MinSeedLength = 0;
	bPacBioData = false;
	bMultiHit = false;
	FastQFormat = true;
	OutputFileFormat = 0; // 0:sam 1:sam.gz, 2:bam
	RefSequence = IndexFileName = ReadFileName = ReadFileName2 = OutputFileName = NULL;

	if (argc == 1 || strcmp(argv[1], "-h") == 0) ShowProgramUsage(argv[0]);
	else
	{
		for (i = 1; i < argc; i++)
		{
			parameter = argv[i];

			if (parameter == "-i") IndexFileName = argv[++i];
			else if (parameter == "-f" || parameter == "-q") ReadFileName = argv[++i];
			else if (parameter == "-f2" || parameter =="-q2") ReadFileName2 = argv[++i];
			else if (parameter == "-t")
			{
				if ((iThreadNum = atoi(argv[++i])) > 16)
				{
					fprintf(stderr, "Warning! Thread number is limited to 16!\n");
					iThreadNum = 16;
				}
			}
			else if (parameter == "-g")
			{
				if((MaxGaps = atoi(argv[++i])) < 0) MaxGaps = 0;
			}
			else if (parameter == "-o") OutputFileName = argv[++i];
			else if (parameter == "-pacbio") bPacBioData = true;
			else if (parameter == "-m") bMultiHit = true;
			else if (parameter == "-pair" || parameter == "-p") bPairEnd = true;
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else
			{
				fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
				ShowProgramUsage(argv[0]);
				exit(0);
			}
		}

		if (IndexFileName == NULL || ReadFileName == NULL)
		{
			fprintf(stderr, "Warning! Please specify a valid index prefix and read files!\n");
			ShowProgramUsage(argv[0]);
			exit(0);
		}
		gzCompressed = Check_gzInputFormat();

		if ((gzCompressed && CheckCompressedFile(ReadFileName) == false) ||
			(!gzCompressed && CheckReadFile(ReadFileName, FastQFormat) == false))
			fprintf(stderr, "Cannot open the read file: %s\n", ReadFileName), exit(0);

		if (ReadFileName2 != NULL)
		{
			bool FastQFormat2 = true;
			if ((gzCompressed && CheckCompressedFile(ReadFileName2) == false) ||
				(!gzCompressed && CheckReadFile(ReadFileName2, FastQFormat2) == false))
				fprintf(stderr, "Cannot open the read file: %s\n", ReadFileName2), exit(0);

			if(FastQFormat2 != FastQFormat) fprintf(stderr, "The input files are not with the same format!\n"), exit(0);
		}

		if (OutputFileName != NULL && CheckOutputFileName() == false) exit(0);

		if (CheckBWAIndexFiles(IndexFileName)) RefIdx = bwa_idx_load(IndexFileName);
		else RefIdx = 0;

		if (RefIdx == 0) fprintf(stderr, "\n\nError! Index files are corrupt!\n");
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
