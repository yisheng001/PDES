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
	return b + dt * (a - c) / (2 * h);
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
	return b + dt * (b - c) / h;
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
void Problem_3(int k)
{
	int K = k;//K=1时为FTCS，K=2时为Lax_Friedrich,K=3时为Lax_Wendroff
	int J = 80;
	double h = 1.0 / J;
	double lamda = 0.5;
	double dt = lamda * h;
	int N = 1.0 / dt;
	double t1 = 0.2, t2 = 0.5;
	int n1 = t1 / dt, n2 = t2 / dt;
	vector<vector<double>>v(N + 1, vector<double>(J + 1));
	for (int j = 0; j <= J; j++)
	{
		v[0][j] = sin(2 * PI * j * h);
	}
	for (int i = 1; i <= N; i++)
	{
		for (int j = 0; j <= J; j++)
		{
			if (j == 0 || j == J)
			{
				if (j == 0)
				{
					v[i][j] = FD(v[i - 1][j + 1], v[i - 1][j], v[i - 1][j - 1 + J], dt, h, K);
				}
				else
				{
					v[i][j] = FD(v[i - 1][j + 1 - J], v[i - 1][j], v[i - 1][j - 1], dt, h, K);
				}
			}
			else
			{
				v[i][j] = FD(v[i - 1][j + 1], v[i - 1][j], v[i - 1][j - 1], dt, h, K);
			}
		}
	}

	long beishu_x = 100, beishu_y = 300;
	int graphic_x = 1600, graphic_y = 1000;//窗口大小

	initgraph(graphic_x, graphic_y);//初始化绘图界面

	setbkcolor(WHITE);
	cleardevice();
	//画坐标
	setlinecolor(LIGHTBLUE);//设置画线颜色

	line(0, graphic_y / 2, graphic_x, graphic_y / 2);
	line(graphic_x / 2, 0, graphic_x / 2, graphic_y);

	settextcolor(RED);//设置字体颜色
	outtextxy(graphic_x - 20, graphic_y / 2 + 5, 'y');
	outtextxy(graphic_x / 2 + 5, 0, 'x');

	settextcolor(LIGHTGREEN);
	outtextxy(graphic_x / 10, graphic_x / 10, L"T1红色 T2绿色 T3蓝色 T4黑色 的图像");

	//setaspectratio(0.005, -1);//这个函数用于设置当前缩放因子。
	setorigin(graphic_x / 2, graphic_y / 2);// 设置坐标原点到屏幕中央点

	POINT* p0 = new POINT[J];
	for (int j = 0; j <= J - 1; j++)
	{
		p0[j].x = beishu_x * j * h;
		p0[j].y = beishu_y * v[n1][j];
	}
	setlinecolor(BLUE);
	polyline(p0, J);

	POINT* p1 = new POINT[J];
	for (int j = 0; j <= J - 1; j++)
	{
		p1[j].x = beishu_x * j * h;
		p1[j].y = beishu_y * v[n2][j];
	}
	setlinecolor(GREEN);
	polyline(p1, J);

	POINT* p2 = new POINT[J];
	for (int j = 0; j <= J - 1; j++)
	{
		p2[j].x = beishu_x * j * h;
		p2[j].y = beishu_y * sin(2 * PI * (j * h + t1));
	}
	setlinecolor(RED);
	polyline(p2, J);

	POINT* p3 = new POINT[J];
	for (int j = 0; j <= J - 1; j++)
	{
		p3[j].x = beishu_x * j * h;
		p3[j].y = beishu_y * sin(2 * PI * (j * h + t2));
	}
	setlinecolor(BLACK);
	polyline(p3, J);

	/*vector<double>vmax(N);
	double Max = 0;
	for (int i = 0; i <= N - 1; i++)
	{
		for (int j = 0; j <= J; j++)
		{
			v[i][j] = v[i][j] - sin(2 * PI * (j * h + i * dt));
			if (v[i][j] < 0)
			{
				v[i][j] = -v[i][j];
			}
		}
		Max = 0;
		for (int j = 0; j <= J; j++)
		{
			if (v[i][j] >= Max)
			{
				Max = v[i][j];
			}
		}
		vmax[i] = Max;
	}

	POINT* p4 = new POINT[N];
	for (int j = 0; j <= N - 1; j++)
	{
		p4[j].x = beishu_x * j * dt;
		p4[j].y = beishu_y * vmax[j];
	}
	setlinecolor(LIGHTBLUE);
	polyline(p4, N);*/

	_getch();
	closegraph();//关闭图形界面

	delete[]p0, p1, p2, p3;
}
void Problem_1(int J, double lamda)
{
	int K = 3;
	double h = 1.0 / J;
	double dt = lamda * h;
	int N = 1.0 / dt;
	double T = 1.0;
	int n = T / dt;
	vector<vector<double>>v(N + 1, vector<double>(J + 1));
	for (int j = 0; j <= J; j++)
	{
		v[0][j] = sin(2 * PI * j * h);
	}
	for (int i = 1; i <= N; i++)
	{
		if (i == 1)
		{
			K = 5;
			for (int j = 0; j <= J; j++)
			{
				if (j == 0 || j == J)
				{
					if (j == 0)
					{
						v[i][j] = FD(v[i - 1][j + 1], v[i - 1][j], v[i - 1][j - 1 + J], dt, h, K);
					}
					else
					{
						v[i][j] = FD(v[i - 1][j + 1 - J], v[i - 1][j], v[i - 1][j - 1], dt, h, K);
					}
				}
				else
				{
					v[i][j] = FD(v[i - 1][j + 1], v[i - 1][j], v[i - 1][j - 1], dt, h, K);
				}
			}
		}
		else
		{
			K = 6;
			for (int j = 0; j <= J; j++)
			{
				if (j == 0 || j == J)
				{
					if (j == 0)
					{
						v[i][j] = FD(v[i - 1][j + 1], v[i - 1][j], v[i - 1][j - 1 + J], dt, h, K);
					}
					else
					{
						v[i][j] = FD(v[i - 1][j + 1 - J], v[i - 1][j], v[i - 1][j - 1], dt, h, K);
					}
				}
				else
				{
					v[i][j] = FD(v[i - 1][j + 1], v[i - 1][j], v[i - 1][j - 1], dt, h, K);
				}
			}
		}
	}

	long beishu_x = 400, beishu_y = 300;
	int graphic_x = 1600, graphic_y = 1000;//窗口大小

	initgraph(graphic_x, graphic_y);//初始化绘图界面

	setbkcolor(WHITE);
	cleardevice();
	//画坐标
	setlinecolor(LIGHTBLUE);//设置画线颜色

	line(0, graphic_y / 2, graphic_x, graphic_y / 2);
	line(graphic_x / 2, 0, graphic_x / 2, graphic_y);

	settextcolor(RED);//设置字体颜色
	outtextxy(graphic_x - 20, graphic_y / 2 + 5, 'y');
	outtextxy(graphic_x / 2 + 5, 0, 'x');

	settextcolor(LIGHTGREEN);
	outtextxy(graphic_x / 10, graphic_x / 10, L"T1红色 T2绿色 T3蓝色 T4黑色 的图像");

	//setaspectratio(0.005, -1);//这个函数用于设置当前缩放因子。
	setorigin(graphic_x / 2, graphic_y / 2);// 设置坐标原点到屏幕中央点

	POINT* p0 = new POINT[J];
	for (int j = 0; j <= J - 1; j++)
	{
		p0[j].x = beishu_x * j * h;
		p0[j].y = beishu_y * v[n][j];
	}
	setlinecolor(GREEN);
	polyline(p0, J);

	POINT* p1 = new POINT[160];
	for (int j = 0; j <= 160 - 1; j++)
	{
		p1[j].x = beishu_x * j * 1.0 / 160;
		p1[j].y = beishu_y * sin(2 * PI * (j * 1.0 / 160));
	}
	setlinecolor(RED);
	polyline(p1, 160);

	_getch();
	closegraph();//关闭图形界面

	delete[]p0, p1;
}
int main()
{
	//Problem_1(160, 0.5);//PG4_2_1 andPG4_2_2
	Problem_3(4);//PG4_2_3,K=4时为FTBS格式
}