#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
int * initial(int N, int K, int *sol)
{
	int i;
	for (i = 0; i<N; i++)
	{
		sol[i] = rand() % K;
	}
	return sol;
}
int tabu_search(int * * G, int N, int K, int maxIter, int *sol)
{
	/*adjacent_color_table  tabu_table*/
	int **act = new int*[N], **tabu = new int *[N];
	/*moves_table  location_in_moves_table numbers  bucket delt_of_move*/
	int ***m = new int**[9], ***loc = new int**[N], n[9], buc, delt;

	int f = 0, best_f = 0;
	int iter = 0;
	int i, j, k;

	/*initial*/
	for (i = 0; i < 9; i++)
	{
		m[i] = new int*[N*K];
		n[i] = 0;
		for (j = 0; j < N*K; j++)
		{
			m[i][j] = new int[2];
		}
	}

	for (i = 0; i<N; i++)
	{
		act[i] = new int[K];
		tabu[i] = new int[K];
		loc[i] = new int*[K];
		for (j = 0; j<K; j++)
		{
			act[i][j] = 0;
			tabu[i][j] = 0;
			loc[i][j] = new int[9];
			for (k = 0; k < 9; k++)
			{
				loc[i][j][k] = -1;
			}
		}
	}
	for (i = 0; i<N; i++)
	{
		for (j = 0; G[i][j] != -1; j++)
		{
			act[i][sol[G[i][j]]] += 1;
		}
	}
	for (i = 0; i<N; i++)
	{
		f += act[i][sol[i]];
		if (act[i][sol[i]])
		{
			for (j = 0; j < K; j++)
			{
				if (j == sol[i]) continue;
				delt = act[i][j] - act[i][sol[i]];
				if (delt < -4) buc = 0;
				else if (delt > 4) buc = 8;
				else buc = delt + 4;
				if (loc[i][j][buc] == -1)
				{
					m[buc][n[buc]][0] = i;
					m[buc][n[buc]][1] = j;
					loc[i][j][buc] = n[buc];
					n[buc]++;
				}
			}
		}
	}
	f /= 2;
	best_f = f;

	/*best_nontabu_move_vetex color delt_f flag count*/
	int bmv, bmc, bmd, bmf, bmcnt;
	/*best_tabu_move_vertex color delt_f flag count*/
	int btmv, btmc, btmd, btmf;
	/*temp_vetex color old_color old_color_of_bmv neighbors act loc*/
	int v, c, tc, tbc, *ngb, *t_act, **t_loc;
	/*random number*/
	int r;

	for (iter = 0; iter<maxIter && f>0; iter++)
	{

		/*findMove*/
		bmf = 0;
		btmf = 0;
		bmd = N;
		btmd = N;
		bmcnt = 0;
		for (i = 0; !bmf &&i < 9; i++)
		{
			for (j = 0; !bmf&&j < n[i]; j++)
			{
				r = rand() % n[i];
				v = m[i][r][0];
				c = m[i][r][1];
				delt = act[v][c] - act[v][sol[v]];
				if (delt < -4) buc = 0;
				else if (delt > 4) buc = 8;
				else buc = delt + 4;
				if (buc != i || c == sol[v] || act[v][sol[v]] == 0)
				{
					m[i][r][0] = m[i][n[i] - 1][0];
					m[i][r][1] = m[i][n[i] - 1][1];
					loc[m[i][r][0]][m[i][r][1]][i] = r;
					loc[v][c][i] = -1;
					n[i]--;
					j--;
					continue;
				}
				if (tabu[v][c] <= iter)
				{
					bmv = v;
					bmc = c;
					bmd = delt;
					bmf = 1;
				}
				else if (!btmf)
				{
					btmv = v;
					btmc = c;
					btmd = delt;
					btmf = 1;
				}
			}
			for (j = 0; j < n[i]; j++)
			{
				v = m[i][j][0];
				c = m[i][j][1];
				if (tabu[v][c] > iter) continue;
				delt = act[v][c] - act[v][sol[v]];
				if (delt < -4) buc = 0;
				else if (delt > 4) buc = 8;
				else buc = delt + 4;
				if (buc != i || c == sol[v] || act[v][sol[v]] == 0)
				{
					m[i][j][0] = m[i][n[i] - 1][0];
					m[i][j][1] = m[i][n[i] - 1][1];
					loc[m[i][j][0]][m[i][j][1]][i] = j;
					loc[v][c][i] = -1;
					n[i]--;
					j--;
					continue;
				}
				bmcnt++;
				if (rand() % bmcnt == 0)
				{
					bmv = v;
					bmc = c;
					bmd = delt;
					bmf = 1;
				}
			}
		}
		/*aspiration*/
		if (btmd < bmd && f + btmd < best_f)
		{
			bmd = btmd;
			bmv = btmv;
			bmc = btmc;
		}

		/*makeMove*/
		ngb = G[bmv];
		tbc = sol[bmv];
		for (i = 0; (v = ngb[i]) != -1; i++)
		{
			t_act = act[v];
			tc = sol[v];
			t_act[tbc]--;
			t_act[bmc]++;
			if (t_act[tc] == 0) continue;
			t_loc = loc[v];
			if (tc == bmc || tc == tbc)
			{
				for (j = 0; j < K; j++)
				{
					if (j == tc) continue;
					delt = t_act[j] - t_act[tc];
					if (delt < -4) buc = 0;
					else if (delt > 4) buc = 8;
					else buc = delt + 4;
					if (t_loc[j][buc] == -1)
					{
						m[buc][n[buc]][0] = v;
						m[buc][n[buc]][1] = j;
						t_loc[j][buc] = n[buc];
						n[buc]++;
					}
				}
			}
			else
			{
				delt = t_act[bmc] - t_act[tc];
				if (delt < -4) buc = 0;
				else if (delt > 4) buc = 8;
				else buc = delt + 4;
				if (t_loc[bmc][buc] == -1)
				{
					m[buc][n[buc]][0] = v;
					m[buc][n[buc]][1] = bmc;
					t_loc[bmc][buc] = n[buc];
					n[buc]++;
				}

				delt = t_act[tbc] - t_act[tc];
				if (delt < -4) buc = 0;
				else if (delt > 4) buc = 8;
				else buc = delt + 4;
				if (t_loc[tbc][buc] == -1)
				{
					m[buc][n[buc]][0] = v;
					m[buc][n[buc]][1] = tbc;
					t_loc[tbc][buc] = n[buc];
					n[buc]++;
				}
			}
		}
		t_act = act[bmv];
		t_loc = loc[bmv];
		if (t_act[bmc]) for (i = 0; i < K; i++)
		{
			if (i == bmc) continue;
			delt = t_act[i] - t_act[bmc];
			if (delt < -4) buc = 0;
			else if (delt > 4) buc = 8;
			else buc = delt + 4;
			if (t_loc[i][buc] == -1)
			{
				m[buc][n[buc]][0] = bmv;
				m[buc][n[buc]][1] = i;
				t_loc[i][buc] = n[buc];
				n[buc]++;
			}
		}

		tabu[bmv][tbc] = iter + rand() % 10 + f + bmd;
		sol[bmv] = bmc;
		f += bmd;
		if (f<best_f)
		{
			best_f = f;
		}
	}

	for (i = 0; i < N; i++)
	{
		delete(act[i]);
		delete(tabu[i]);
		for (j = 0; j < K; j++)
		{
			delete loc[i][j];
		}
	}
	delete act;
	delete tabu;
	for (i = 0; i < 9; i++)
	{
		for (j = 0; j < N*K; j++)
		{
			delete (m[i][j]);
		}
		delete(m[i]);
	}
	return best_f;
}
int * * input(char * fileName, int *pn)
{
	int m, u, v, i, j;
	char s[100];
	int * * G = 0;
	FILE *pf;
	pf = fopen(fileName, "r");
	for (i = 0; i<12; i++)
		fgets(s, 100, pf);
	fscanf(pf, "p edge %d %d", pn, &m);
	int *n = new int[*pn];
	G = new int*[*pn];
	for (i = 0; i < *pn; i++)
		n[i] = 0;
	for (i = 0; i<m; i++)
	{
		fgets(s, 100, pf);
		fscanf(pf, "e %d %d", &u, &v);
		u--;
		v--;
		n[u]++;
		n[v]++;
	}
	fclose(pf);
	for (i = 0; i<*pn; i++)
	{
		G[i] = new int[n[i] + 1];
		G[i][n[i]] = -1;
		n[i] = 0;
	}
	pf = fopen(fileName, "r");
	for (i = 0; i<12; i++)
		fgets(s, 100, pf);
	for (i = 0; i<m; i++)
	{
		fgets(s, 100, pf);
		fscanf(pf, "e %d %d", &u, &v);
		u--;
		v--;
		G[u][n[u]] = v;
		n[u]++;
		G[v][n[v]] = u;
		n[v]++;
	}
	fclose(pf);
	return G;
}
void end(int ** G, int n)
{
	int i;
	for (i = 0; i<n; i++)
	{
		delete(G[i]);
	}
	delete G;
}
int *crossover(int N, int K, int *p1, int *p2, int *child)
{
	int *t[2] = { new int[N],new int[N] }, *cnt[2] = { new int[K],new int[K] };
	int max = 0;
	int i, j;
	for (i = 0; i<K; i++)
	{
		cnt[0][i] = 0;
		cnt[1][i] = 0;
	}
	for (i = 0; i<N; i++)
	{
		t[0][i] = p1[i];
		t[1][i] = p2[i];
		cnt[0][p1[i]]++;
		cnt[1][p2[i]]++;
	}
	for (i = 0; i<K; i++)
	{
		for (j = 0; j<K; j++)
			if (cnt[i % 2][j] >= cnt[i % 2][max])
				max = j;
		cnt[i % 2][max] = 0;
		for (j = 0; j<N; j++)
		{
			if (t[i % 2][j] == max)
			{
				child[j] = i;
				cnt[(i + 1) % 2][t[(i + 1) % 2][j]]--;
				t[0][j] = -1;
				t[1][j] = -1;
			}
		}
	}
	for (i = 0; i<N; i++)
	{
		if (t[0][i] != -1)
			child[i] = rand() % K;
	}
	return child;
}
int findPath(int v1, int **weight, int *l1, int* l2, int * mat1, int *mat2, int *path1, int *path2, int K)
{
	int i;
	path1[v1] = 1;
	for (i = 0; i<K; i++)
	{
		if (path2[i] == 0
			&& l1[v1] + l2[i] == weight[v1][i])
		{
			path2[i] = 1;
			if (mat2[i] == -1 || findPath(mat2[i], weight, l1, l2, mat1, mat2, path1, path2, K))
			{
				mat2[i] = v1;
				mat1[v1] = i;
				return 1;
			}
		}
	}
	return 0;
}
int distance(int*r1, int*r2, int N, int K)
{
	int **weight = new int*[K];
	int *l1 = new int[K], *l2 = new int[K];
	int *mat1 = new int[K], *mat2 = new int[K];
	int *path1 = new int[K], *path2 = new int[K];
	int sum = 0;
	int i, j, k, delt;
	for (i = 0; i<K; i++)
	{
		l1[i] = 0;
		l2[i] = 0;
		mat1[i] = -1;
		mat2[i] = -1;
		path1[i] = 0;
		path2[i] = 0;
		weight[i] = new int[K];
		for (j = 0; j<K; j++)
		{
			weight[i][j] = 0;
		}
	}
	for (i = 0; i<N; i++)
	{
		weight[r1[i]][r2[i]]++;
	}
	for (i = 0; i<K; i++)
	{
		for (j = 0; j<K; j++)
		{
			if (weight[i][j]>l1[i])
				l1[i] = weight[i][j];
		}
	}

	for (i = 0; i<K; i++)
	{
		if (mat1[i] != -1) continue;
		while (1)
		{
			for (j = 0; j<K; j++)
				path1[j] = path2[j] = 0;
			if (findPath(i, weight, l1, l2, mat1, mat2, path1, path2, K))
				break;
			delt = N;
			for (j = 0; j<K; j++)
			{
				if (path1[j])
					for (k = 0; k<K; k++)
					{
						if (!path2[k] && delt>l1[j] + l2[k] - weight[j][k])
							delt = l1[j] + l2[k] - weight[j][k];
					}
			}
			for (j = 0; j<K; j++)
			{
				if (path1[j])
					l1[j] -= delt;
				if (path2[j])
					l2[j] += delt;
			}
		}
	}
	for (i = 0; i<K; i++)
	{
		sum += weight[i][mat1[i]];
	}

	return N - sum;
}
int HEA(int ** G, int N, int K, int *result, FILE *pout)
{
	int Ncycle = 10;
	int *p1 = new int[N], *p2 = new int[N], *c1, *c2, *elite1 = new int[N], *elite2 = new int[N], *best,*t;
	int fp1, fp2, fe1, fe2, fbest;
	int cycle=0,gnr = 1, maxGnr = 4000, maxIter = 50000;
	int i, j;

	time_t t1, t2;
	t1 = time(NULL);
	srand(t1);
	fprintf(pout, "%lld,%d,%d,", t1, Ncycle, maxIter);

	fp1 = tabu_search(G, N, K, maxIter, initial(N, K, p1));
	fbest = fp1,best=p1;
	fp2 = tabu_search(G, N, K, maxIter, initial(N, K, p2));
	if (fp2 < fbest) fbest = fp2,best=p2;
	fe1 = tabu_search(G, N, K, maxIter, initial(N, K, elite1));
	if (fe1 < fbest) fbest = fe1,best=elite1;
	fe2 = tabu_search(G, N, K, maxIter, initial(N, K, elite2));
	if (fe2 < fbest) fbest = fe2,best=elite2;
	
	for (; fbest>0 && gnr<maxGnr; gnr++)
	{
		t = NULL;
		c1 = new int[N];
		c2 = new int[N];
		fp1= tabu_search(G, N, K, maxIter, crossover(N, K, p1, p2, c1));
		fp2= tabu_search(G, N, K, maxIter, crossover(N, K, p2, p1, c2));
		if (fp1 < fe1)
		{
			t = c1;
			fe1 = fp1;
		}
		if (fp2 < fe1)
		{
			t = c2;
			fe1 = fp2;
		}
		if(t)
			for (i = 0; i < N; i++)
			{
				elite1[i] = t[i];
			}
		delete p1;
		delete p2;
		p1 = c1;
		p2 = c2;
		printf("generation %d:f=%d  %d\n", gnr, fp1,fp2);

		if (!(gnr%Ncycle))
		{
			printf("cycle %d:felite1=%d\n", cycle, fe1);
			if (fe1 < fbest) fbest = fe1,best=elite1;
			if (fp1 > fp2)
			{
				delete p1;
				p1 = elite2;
			}
			else
			{
				delete p2;
				p2 = elite2;
			}
			elite2 = elite1;
			elite1 = new int[N];
			fe1 = tabu_search(G, N, K, maxIter, initial(N, K, elite1));
			if (fe1 < fbest) fbest = fe1,best=elite1;
			cycle++;
		}
	}
	
	t2 = time(NULL);
	fprintf(pout, "%d,%f,%d,", gnr, difftime(t2, t1), fbest);
	for (i = 0; i<N; i++)
	{
		result[i] = best[i];
		fprintf(pout, "%d ", result[i]);
	}
	fprintf(pout, "\n");
	return fbest;
}
void check(int *result, int **G, int N, int K)
{
	int i, j, flag = 1;
	for (i = 0; i<N; i++)
	{
		for (j = 0; G[i][j] != -1; j++)
		{
			if (result[i] == result[G[i][j]])
				flag = 0;
		}
	}
	printf("check:\n");
	if (flag)
		printf("correct\n");
	else printf("wrong\n");

	int *use = new int[K], count = 0;
	for (i = 0; i<K; i++)
		use[K] = 0;
	for (i = 0; i<N; i++)
		use[result[i]] = 1;
	for (i = 0; i<K; i++)
		if (use[i]) count++;
	printf("use %d colors\n", count);
	return;
}
int main()
{
	char filename[100] = "DSJC500.5.col";
	int N, K = 48, f;
	int i;

	FILE *pout;
	pout = fopen("record.csv", "at");
	fprintf(pout, "%s,%d,", filename, K);

	int** G;
	G = input(filename, &N);

	int *result = new int[N];
	f = HEA(G, N, K, result, pout);

	check(result, G, N, K);
	printf("f=%d\n", f);
	end(G, N);
	fclose(pout);
	return 0;
}
