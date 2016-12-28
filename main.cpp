#include "structure.h"

bwt_t *Refbwt;
bwaidx_t *RefIdx;
const char* VersionStr = "2.1.0";
int iThreadNum, MaxInsertSize, MaxGaps, MinSeedLength;
bool bDebugMode, bPairEnd, FastQFormat, bPacBioData, bMultiHit;
char *RefSequence, *IndexFileName, *ReadFileName, *ReadFileName2, *SamFileName;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\nKart v%s [Developers: Hsin-Nan Lin and Wen-Lian Hsu]\n\n", VersionStr);
	fprintf(stderr, "Usage: %s aln|index\n\n", program);
	fprintf(stderr, "Command: index		index the reference sequences with FASTA format\n");
	fprintf(stderr, "         aln		read alignment\n");
	fprintf(stderr, "\n");
}

bool CheckReadFile(char* filename, bool bFirst)
{
	char ch;
	fstream file;
	bool bCheck = true;
	string str, header, seq;
	int iCount = 0, iBaseCount = 0;

	file.open(filename, ios_base::in);
	if (!file.is_open()) return false;
	else
	{
		getline(file, header);
		if (header == "") return false;
		else
		{
			if (bFirst)
			{
				if (header[0] == '>') FastQFormat = false;
				else FastQFormat = true;
			}
			else
			{
				if ((FastQFormat && header[0] != '@') || (!FastQFormat && header[0] != '>')) return false;
			}

			//if (FastQFormat) getline(file, seq);
			//else {
			//	while (!file.eof())
			//	{
			//		getline(file, str); seq += str;
			//		file.get(ch); file.unget();
			//		if (ch == '>')
			//		{
			//			iBaseCount += (int)seq.length();
			//			seq.clear();
			//			if (++iCount == 100) break;
			//		}
			//	}
			//}
			//if (seq.length() > 0)
			//{
			//	iCount++;
			//	iBaseCount += (int)seq.length();
			//}
			//if (!bPacBioData && iCount > 0 && iBaseCount / iCount > 1000)
			//{
			//	fprintf(stderr, "Warning! The input library contains long reads, sensitive mode (-pacbio) is switched on!\n");
			//	bPacBioData = true; bPairEnd = false;
			//}
		}
	}
	file.close();

	return bCheck;
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	MaxGaps = 5;
	iThreadNum = 16;
	bPairEnd = false;
	bDebugMode = false;
	FastQFormat = true; // fastq:true, fasta:false
	MaxInsertSize = 1500;
	MinSeedLength = 0;
	bPacBioData = false;
	bMultiHit = false;
	RefSequence = IndexFileName = ReadFileName = ReadFileName2 = SamFileName = NULL;

	if (argc == 1) ShowProgramUsage(argv[0]);
	else if (strcmp(argv[1], "index") == 0)
	{
		string cmd(argv[0]);
		cmd = cmd.substr(0, cmd.find_last_of('/')+1) + "bwt_index";

		for (i = 2; i < argc; i++) cmd += " " + (string)argv[i];
		system((char*)cmd.c_str());
	}
	else if (strcmp(argv[1], "aln") == 0)
	{
		for (i = 2; i < argc; i++)
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
			else if (parameter == "-o") SamFileName = argv[++i];
			else if (parameter == "-pacbio") bPacBioData = true;
			else if (parameter == "-m") bMultiHit = true;
			else if (parameter == "-pair" || parameter == "-p") bPairEnd = true;
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
		}

		if (IndexFileName == NULL || ReadFileName == NULL)
		{
			fprintf(stderr, "\n");
			fprintf(stderr, "Kart v%s\n", VersionStr);
			fprintf(stderr, "Usage: %s aln [-i IndexFile Prefix] -f ReadFile [-f2 ReadFile2] > out.sam\n\n", argv[0]);
			fprintf(stderr, "Options: -t INT        number of threads [16]\n");
			fprintf(stderr, "         -f            files with #1 mates reads\n");
			fprintf(stderr, "         -f2           files with #2 mates reads\n");
			fprintf(stderr, "         -o            sam filename for output\n");
			fprintf(stderr, "         -m            output multiple alignments\n");
			fprintf(stderr, "         -g INT        max gaps (indels)\n");
			fprintf(stderr, "         -p            paired-end reads are interlaced in the same file\n");
			fprintf(stderr, "         -pacbio       pacbio data\n");
			fprintf(stderr, "\n");
			exit(0);
		}
		if (CheckReadFile(ReadFileName, true) == false) fprintf(stderr, "Cannot open the read file: %s\n", ReadFileName), exit(0);
		if (ReadFileName2 != NULL && CheckReadFile(ReadFileName2, false) == false) fprintf(stderr, "Read file: %s cannot be accessed or with incompatible format!\n", ReadFileName2), exit(0);

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
	else ShowProgramUsage(argv[0]);

	return 0;
}
