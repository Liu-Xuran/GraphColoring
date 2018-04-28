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
	/*conflict_vetexs  location_in_arrey numbers*/
	int *v = new int[N], *loc = new int[N], n = 0;
	int f = 0, best_f = 0;
	int iter = 0;
	int i, j;

	/*initial*/
	for (i = 0; i<N; i++)
	{
		loc[i] = -1;
		act[i] = new int[K];
		tabu[i] = new int[K];
		for (j = 0; j<K; j++)
		{
			act[i][j] = 0;
			tabu[i][j] = 0;
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
			v[n] = i;
			loc[i] = n;
			n++;
		}
	}
	f /= 2;
	best_f = f;

	/*best_nontabu_move best_tabu_move*/
	int bmd, bmv, bmc, btmd, btmv, btmc;
	/*temp_vertex temp_color temp_act temp_tabu temp_neighbor*/
	int tv, tc, *t_act, *t_tabu, *t_n;
	int delt, same, samet;
	int c;
	for (iter = 0; iter<maxIter && f>0; iter++)
	{

		/*findMove*/
		bmd = N;
		btmd = N;
		for (i = 0; i<n; i++)
		{
			tv = v[i];
			tc = sol[tv];
			t_act = act[tv];
			t_tabu = tabu[tv];
			for (j = 0; j<K; j++)
			{
				if (j == tc) continue;
				delt = t_act[j] - t_act[tc];
				if (t_tabu[j] <= iter)
				{
					if (delt<bmd)
					{
						bmv = tv;
						bmc = j;
						bmd = delt;
						same = 1;
					}
					else if (delt == bmd)
					{
						same++;
						if (rand() % same == 0)
						{
							bmv = tv;
							bmc = j;
							bmd = delt;
						}
					}
				}
				else
				{
					if (delt<btmd)
					{
						btmv = tv;
						btmc = j;
						btmd = delt;
						samet = 1;
					}
					else if (delt == btmd)
					{
						samet++;
						if (rand() % samet == 0)
						{
							btmv = tv;
							btmc = j;
							btmd = delt;
						}
					}
				}
			}
		}
		/*aspiration*/
		if (btmd<bmd && f + btmd < best_f)
		{
			bmd = btmd;
			bmv = btmv;
			bmc = btmc;
		}

		/*makeMove*/
		t_n = G[bmv];
		tc = sol[bmv];
		for (i = 0; (tv = t_n[i]) != -1; i++)
		{
			act[tv][tc]--;
			act[tv][bmc]++;
			c = sol[tv];
			if (c == tc && act[tv][c] == 0)
			{
				v[loc[tv]] = v[n - 1];
				loc[v[n - 1]] = loc[tv];
				loc[tv] = -1;
				n--;
			}
			else if (c == bmc&&act[tv][c] == 1)
			{
				v[n] = tv;
				loc[tv] = n;
				n++;
			}
		}
		if (act[bmv][bmc] == 0)
		{
			v[loc[bmv]] = v[n - 1];
			loc[v[n - 1]] = loc[bmv];
			loc[bmv] = -1;
			n--;
		}

		tabu[bmv][tc] = iter + rand() % 10 + f + bmd;
		sol[bmv] = bmc;
		f += bmd;
		if (f<best_f)
		{
			best_f = f;
		}
	}

	for (i = 0; i < N; i++)
	{
		delete act[i];
		delete tabu[i];
	}
	delete act;
	delete tabu;
	delete v;
	delete loc;

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
	int **weight=new int*[K];
	int *l1=new int[K], *l2 = new int[K];
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
	int population = 10;
	int **pop=new int*[population], *f=new int[population], **d = new int*[population], *sumD=new int[population];
	int best = 0, worst = 0;
	int *newpop, newf, *newD=new int[population], newSum;
	int gnr = 0, maxGnr = 4000, maxIter = 50000;
	int p1, p2;
	int i, j, t;

	time_t t1, t2;
	t1 = time(NULL);
	srand(t1);
	fprintf(pout, "%lld,%d,%d,", t1, population, maxIter);

	for (i = 0; i<population; i++)
	{
		sumD[i] = 0;
		pop[i] = new int [N];
		d[i] = new int[population];
		f[i] = tabu_search(G, N, K, maxIter, initial(N, K, pop[i]));
		printf("initial %d:f=%d\n", i, f[i]);
		if (f[i] <= f[best]) best = i;
		if (f[i] == 0) goto finish;
	}
	for (i = 0; i<population; i++)
	{
		for (j = i; j<population; j++)
		{
			d[i][j] = d[j][i] = distance(pop[i], pop[j], N, K);
			sumD[i] += d[i][j];
			sumD[j] += d[i][j];
		}
	}


	for (; f[best]>0 && gnr<maxGnr; gnr++)
	{
		/*choose parent*/
		p1 = rand() % population;
		p2 = (p1 + 1 + rand() % (population - 1)) % population;
		/*generate new population*/
		newpop = (int *)malloc(N * sizeof(int));
		newf = tabu_search(G, N, K, maxIter, crossover(N, K, pop[p1], pop[p2], newpop));
		newSum = 0;
		for (i = 0; i<population; i++)
		{
			newD[i] = distance(newpop, pop[i], N, K);
			newSum += newD[i];
			if (newD[i] == 0) goto cancel;
		}
		/*update pop*/
		for (i = 0, worst = 0; i<population; i++)
		{
			if (f[i]>f[worst] || (f[i] == f[worst] && sumD[i]<sumD[worst]))
				worst = i;
		}
		if (newf<f[worst] || (newf == f[worst] && newSum >= sumD[worst]))
		{
			free(pop[worst]);
			pop[worst] = newpop;
			f[worst] = newf;
			if (newf<f[best])
			{
				best = worst;
			}
			for (i = 0; i<population; i++)
			{
				sumD[i] = sumD[i] + newD[i] - d[i][worst];
				d[i][worst] = d[worst][i] = newD[i];
			}
			sumD[worst] = newSum - newD[worst];
			d[worst][worst] = 0;
		}

		else
			cancel:
		free(newpop);
		printf("generation %d:f=%d\n", gnr, newf);
	}
	for (i = best + 1; i<population; i++)
	{
		free(pop[i]);
	}
finish:
	t2 = time(NULL);
	fprintf(pout, "%d,%f,%d,", gnr, difftime(t2, t1), f[best]);
	printf("generation=%d\n", gnr);
	for (i = 0; i<N; i++)
	{
		result[i] = pop[best][i];
		fprintf(pout, "%d ", result[i]);
	}
	fprintf(pout, "\n");
	for (i = 0; i <= best; i++)
	{
		free(pop[i]);
	}
	return f[best];
}
void check(int *result, int **G, int N, int K)
{
	int i, j, flag = 1;
	for (i = 0; i<N; i++)
	{
		for (j = 0; G[i][j]!=-1; j++)
		{
			if (result[i] == result[G[i][j]])
				flag = 0;
		}
	}
	printf("check:\n");
	if (flag)
		printf("correct\n");
	else printf("wrong\n");

	int *use=new int[K], count = 0;
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

	int *result=new int[N];
	f = HEA(G, N, K, result, pout);

	check(result, G, N, K);
	printf("f=%d\n", f);
	end(G, N);
	fclose(pout);
	return 0;
}
