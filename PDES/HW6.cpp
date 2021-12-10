#define _CRT_SECURE_NO_WARNINGS//跳过安全检查
#include <graphics.h>//图形界面头文件
#include <conio.h>//标准输入输出头文件
#include<Windows.h>
#include <iostream>
#include <vector>
#include"math.h"
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cmath>
double PI = 3.1415926535;
using namespace std;
double A(int j, int N, int k)
{
	double A_L = 4.0, A_R = 1.0;

	if (j == (N - 1) / 2)
		if (k == 1)
			return double((A_L + A_R) / 2);
		else
			return double(2 * A_L * A_R / (A_L + A_R));
	else if (j < (N - 1) / 2)
		return double(A_L);
	else
		return double(A_R);
}
double U(double t, double x)
{
	if (x <= 0)
		return sin(x) * exp(-4 * t) / 2;
	else
		return sin(2 * x) * exp(-4 * t);
}
int main()
{
	int K = 5;
	int type = 2;//k=1用的是算术平均值，k=2用的是调和平均值；
	double mu = 0.1;
	vector<int>N(K);
	vector<double>LL2(K), LLL(K);
	for (int i = 0; i < K; i++)
	{
		N[i] = 10 * pow(2, i + 1) + 1;
		double dx = 2 * PI / N[i];
		double dt = mu * dx * dx;
		int n = 1 / dt;
		double T = n * dt;
		vector<vector<double>>V(n + 1, vector<double>(N[i] + 1));
		for (int j = 0; j <= N[i]; j++)
		{
			V[0][j] = U(0, j * dx - PI);
		}
		for (int j = 1; j <= n; j++)
		{
			V[j][0] = V[j][N[i]] = 0;
		}
		for (int k = 0; k < n; k++)
		{
			for (int j = 1; j < N[i]; j++)
			{
				V[k + 1][j] = V[k][j] + mu * (A(j, N[i], type) * (V[k][j + 1] - V[k][j]) - A(j - 1, N[i], type) * (V[k][j] - V[k][j - 1]));
			}
		}
		//计算L2范数误差,L无穷误差
		vector<double>L0(N[i] + 1);
		double L2 = 0, LL = 0;
		for (int j = 0; j <= N[i]; j++)
		{
			L0[j] = U(T, j * dx - PI) - V[n][j];
			L2 = L2 + L0[j] * L0[j];
			if (L0[j] < 0)
			{
				L0[j] = -L0[j];
			}
			if (L0[j] > LL)
			{
				LL = L0[j];
			}
		}
		LL2[i] = L2; LLL[i] = LL;
	}
	vector<double>ord_L2(K), ord_LL(K);
	ord_L2[0] = 0; ord_LL[0] = 0;
	for (int i = 1; i < K; i++)
	{
		ord_L2[i] = log2(LL2[i - 1] / LL2[i]);
		ord_LL[i] = log2(LLL[i - 1] / LLL[i]);
	}
	cout << "次数" << "\t" << "L2误差" << "\t" << "L2_ord" << "\t" << "L无穷误差" << "\t" << "LL_ord" << endl;
	for (int i = 0; i < K; i++)
	{
		cout << N[i] << "\t" << LL2[i] << "\t" << ord_L2[i] << "\t" << LLL[i] << "\t" << ord_LL[i] << endl;
	}
}