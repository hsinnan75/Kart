# Changes
- 2.5.6: Fix a bug when the read data contains unexpected characters.
- 2.5.5: Fix several compiler warnings.
- 2.5.4: Integrated bwt_index into Kart.
- 2.5.3: Fix a bug when running with multi-threads on Mac PCs
- 2.5.2: Fix a bug in parameter processing.
- 2.5.1: Dicard the thread limit.
- 2.5.0: Add BAM formatted output.
- 2.4.9: Fix a bug in reading wrong formatted read header
- 2.4.8: Fix a bug when reading long read data in a gzip file.
- 2.4.7: Fix a bug on paired-end mapping.
- 2.4.6: Replaced 0 exit codes with 1 and the corresponding 'Warning' with 'Error' for cases where program termination is not the expected result (revised by Rad Suchecki).
- 2.4.5: Add a '-v' option to show the version of kart. 
- 2.4.4: Fix a bug on resuce mode
- 2.4.3: Fix a bug on reading paired-end inputs.
- 2.4.2: Modify the MAPQ measure for the PacBio read mapping. 
- 2.4.1: Fix a bug on reporting alignment coordinate with soft clipping.
- 2.4.0: Use ksw algorithm for PacBio long read alignment.
- 2.3.9: Fix a bug in SAM output.
- 2.3.8: Fix the alignment of segment pairs with poor sequence identity.
- 2.3.7: Add an update command.
- 2.3.6: Fix the alignment of segment pairs with multiple mismatches.
- 2.3.5: Fix the alignment when DNA sequences are shown in lower case.
- 2.3.4: Fix the alignment of poor-quality sequences at both ends.
- 2.3.3: Remove edlib and ksw algorithms.
- 2.3.2: Allow multiple read files as the input.
- 2.3.1: Fix the buf when read number exceeds 2^32.
- 2.3.0: Add ksw2 and edlib alignment method to replace the Needleman-Wunsch algorithm.
- 2.2.6: Fix a bug in the alignment report.
- 2.2.5: Fix a bug in the alignment report.
