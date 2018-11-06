#include "structure.h"

const float MaxPenalty = -65536;
const float OPEN_GAP = -1;
const float EXTEND_GAP = -0.5;
const float NEW_GAP = -1.5;

double max(float x, float y)
{
	return x > y ? x : y;
}

double max(float x, float y, float z)
{
	return x > y ? max(x, z) : max(y, z);
}

void nw_alignment(int m, string& s1, int n, string& s2)
{
	int i, j;

	m = m + 1, n = n + 1;

	float** r = new float*[m];
	float** t = new float*[m];
	float** s = new float*[m];

	for (i = 0; i < m; i++)
	{
		r[i] = new float[n];
		t[i] = new float[n];
		s[i] = new float[n];
	}

	// initialization
	r[0][0] = t[0][0] = s[0][0] = 0;
	for (i = 1; i < m; i++)
	{
		r[i][0] = MaxPenalty;
		s[i][0] = t[i][0] = OPEN_GAP + i*EXTEND_GAP;
	}

	for (j = 1; j < n; j++)
	{
		t[0][j] = MaxPenalty;
		s[0][j] = r[0][j] = OPEN_GAP + j*EXTEND_GAP;
	}

	for (i = 1; i < m; i++)
	{
		for (j = 1; j < n; j++)
		{
			r[i][j] = max(r[i][j - 1] + EXTEND_GAP, s[i][j - 1] + NEW_GAP);
			t[i][j] = max(t[i - 1][j] + EXTEND_GAP, s[i - 1][j] + NEW_GAP);
			s[i][j] = max(s[i - 1][j - 1] + (nst_nt4_table[s1[i - 1]] == nst_nt4_table[s2[j - 1]] ? 1.5 : -1.5), r[i][j], t[i][j]);
		}
	}
	// back tracking
	i = m - 1, j = n - 1;
	while (i > 0 || j > 0) {
		if (s[i][j] == r[i][j]) {
			s1.insert(i, 1, '-');
			j--;
		}
		else if (s[i][j] == t[i][j]) {
			s2.insert(j, 1, '-');
			i--;
		}
		else {
			i--, j--;
		}
	}
	for (i = 0; i < m; i++)
	{
		delete[] r[i]; 
		delete[] t[i]; 
		delete[] s[i];
	}
	delete[] r; delete[] t; delete[] s;
}
