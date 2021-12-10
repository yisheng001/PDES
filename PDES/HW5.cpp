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
double FTCS(double a, double b, double c, double dt, double h)
{
	return b - dt * (a - c) / (2 * h);
}
double Lax_Friedrich(double a, double b, double c, double dt, double h)
{
	return (a + c) / 2 + dt * (a - c) / (2 * h);
}
double Lax_Wendroff(double a, double b, double c, double dt, double h)
{
	return b + dt * (a - c) / (2 * h) + dt * dt * (a + c - 2 * b) / (2 * h * h);
}
double FTBS(double a, double b, double c, double dt, double h)
{
	return b - dt * (b - c) / h;
}
double FTFS(double a, double b, double c, double dt, double h)
{
	return b + dt * (a - b) / h;
}
double CTCS(double a, double b, double c, double dt, double h)
{
	return c + dt * (a - c) / h;
}
double FD(double a, double b, double c, double dt, double h, int K)
{
	if (K == 1)
	{
		return FTCS(a, b, c, dt, h);
	}
	else if (K == 2)
	{
		return Lax_Friedrich(a, b, c, dt, h);
	}
	else if (K == 3)
	{
		return Lax_Wendroff(a, b, c, dt, h);
	}
	else if (K == 4)
	{
		return FTBS(a, b, c, dt, h);
	}
	else if (K == 5)
	{
		return FTFS(a, b, c, dt, h);
	}
	else
	{
		return CTCS(a, b, c, dt, h);
	}
}
void problem(double r)
{
	double dx = 0.05;
	double dt = dx * r;
	int N = 1.0 / dx;
	int T = N/r;
	int T1 = T, T2 = 2 * T, T3 = 5 * T;
	int N0 = N + 2 * T3;
	vector<vector<double>>V1(T3 + 1, vector<double>(N0)), V2(T3 + 1, vector<double>(N0));
	for (int i = 0; i < N0; i++)
	{
		if (i >= T3+8 && i <= T3 + N-8)
		{
			V2[0][i] = 1;
		}
		else
		{
			V2[0][i] = 0;
		}
	}
	for (int i = 1; i <= T3; i++)
	{
		for (int j = 0; j < N0; j++)
		{
			if (j == 0)
			{
				//V1[i][j] = FD(V1[i - 1][j + 1], V1[i - 1][j], 0, dt, dx, 1);
				V2[i][j] = FD(V2[i - 1][j + 1], V2[i - 1][j], 0, dt, dx, 4);
			}
			else if (j == N0 - 1)
			{
				//V1[i][j] = FD(0, V1[i - 1][j], V1[i - 1][j - 1], dt, dx, 1);
				V2[i][j] = FD(0, V2[i - 1][j], V2[i - 1][j - 1], dt, dx, 4);
			}
			else
			{
				//V1[i][j] = FD(V1[i - 1][j + 1], V1[i - 1][j], V1[i - 1][j - 1], dt, dx, 1);
				V2[i][j] = FD(V2[i - 1][j + 1], V2[i - 1][j], V2[i - 1][j - 1], dt, dx, 4);
			}
		}
	}

	int beishu_x = 100, beishu_y = 400;
	int graphic_x = 1600, graphic_y = 1000;//窗口大小

	initgraph(graphic_x, graphic_y);//初始化绘图界面

	setbkcolor(WHITE);
	cleardevice();
	//画坐标
	setlinecolor(LIGHTBLUE);//设置画线颜色

	line(0, graphic_y / 2, graphic_x, graphic_y / 2);
	line(graphic_x / 2, 0, graphic_x / 2, graphic_y);

	settextcolor(RED);//设置字体颜色
	outtextxy(graphic_x - 20, graphic_y / 2 + 5, 'x');
	outtextxy(graphic_x / 2 + 5, 0, 'y');

	setaspectratio(1, -1);//这个函数用于设置当前缩放因子。
	setorigin(graphic_x / 2, graphic_y / 2);// 设置坐标原点到屏幕中央点

	POINT* p0 = new POINT[N0];
	for (int j = 0; j <= N0 - 1; j++)
	{
		p0[j].x = beishu_x * (j * dx- 0.5 - 5.0 / r);
		p0[j].y = beishu_y * V2[T1][j];
	}
	setlinecolor(GREEN);
	polyline(p0, N0);

	POINT* p1 = new POINT[N0];
	for (int j = 0; j <= N0 - 1; j++)
	{
		p1[j].x = beishu_x * (j * dx - 0.5 - 5.0 / r);
		p1[j].y = beishu_y * V2[T2][j];
	}
	setlinecolor(GREEN);
	polyline(p1, N0);

	POINT* p2 = new POINT[N0];
	for (int j = 0; j <= N0 - 1; j++)
	{
		p2[j].x = beishu_x * (j * dx - 0.5-5.0/r);
		p2[j].y = beishu_y * V2[T3][j];
	}
	setlinecolor(GREEN);
	polyline(p2, N0);

	POINT* p3 = new POINT[N0];
	for (int j = 0; j <= N0 - 1; j++)
	{
		p3[j].x = beishu_x * (j * dx - 0.5 - 5.0 / r +1.0);
		p3[j].y = beishu_y * V2[0][j];
	}
	setlinecolor(RED);
	polyline(p3, N0);

	POINT* p4 = new POINT[N0];
	for (int j = 0; j <= N0 - 1; j++)
	{
		p4[j].x = beishu_x * (j * dx - 0.5 - 5.0 / r +2.0);
		p4[j].y = beishu_y * V2[0][j];
	}
	setlinecolor(RED);
	polyline(p4, N0);

	POINT* p5 = new POINT[N0];
	for (int j = 0; j <= N0 - 1; j++)
	{
		p5[j].x = beishu_x * (j * dx - 0.5 - 5.0 / r +5.0);
		p5[j].y = beishu_y * V2[0][j];
	}
	setlinecolor(RED);
	polyline(p5, N0);

	
	_getch();
	closegraph();//关闭图形界面

	delete[]p0,p1,p2,p3,p4,p5;
}
int main()
{
	problem(0.8);
	return 0;
}