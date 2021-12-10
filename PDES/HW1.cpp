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
void problem(double dx, double dt)
{
	int N, J;
	J = (int)1 / dx;
	N = (int)1 / dt;
	double t0 = 0.3;
	int n;
	n = (int)(t0 / dt);
	vector<vector<double>>v(J + 1, vector<double>(n + 1));
	for (int j = 0; j <= J; j++)
	{
		v[j][0] = sin(2 * PI * j * dx);
	}
	for (int i = 1; i <= n; i++)
	{
		for (int j = 0; j <= J; j++)
		{
			if (j != J)
			{
				v[j][i] = v[j][i - 1] + dt / dx * (v[j + 1][i - 1] - v[j][i - 1]);
			}
			else
			{
				v[j][i] = v[j][i - 1] + dt / dx * (v[1][i - 1] - v[j][i - 1]);
			}
		}
	}

	long beishu_x = 200, beishu_y = 100;
	int graphic_x = 800, graphic_y = 600;//窗口大小

	initgraph(graphic_x, graphic_y);//初始化绘图界面

 //画坐标
	setlinecolor(BLUE);//设置画线颜色
	line(0, graphic_y / 2, graphic_x, graphic_y / 2);
	line(graphic_x / 2, 0, graphic_x / 2, graphic_y);

	settextcolor(RED);//设置字体颜色
	outtextxy(graphic_x - 20, graphic_y / 2 + 5, 'y');
	outtextxy(graphic_x / 2 + 5, 0, 'x');

	settextcolor(LIGHTGREEN);
	outtextxy(graphic_x / 10, graphic_x / 10, L"函数 t=0.3 的近似解(绿），精确解（红）  X属于(0,1) 的图像");

	//setaspectratio(0.005, -1);//这个函数用于设置当前缩放因子。
	setorigin(graphic_x / 2, graphic_y / 2);// 设置坐标原点到屏幕中央点

	POINT* p1 = new POINT[J + 1];

	for (int j = 0; j <= J; j++)
	{
		p1[j].x = beishu_x * j * dx;
		p1[j].y = beishu_y * v[j][n];
	}
	setlinecolor(GREEN);
	polyline(p1, J + 1);

	POINT* p0 = new POINT[J + 1];
	for (int j = 0; j <= J; j++)
	{
		p0[j].x = beishu_x * j * dx;
		p0[j].y = beishu_y * sin(2 * PI * (j * dx + t0));
	}
	setlinecolor(RED);
	polyline(p0, J + 1);

	_getch();
	closegraph();//关闭图形界面

	delete[]p0, p1;
}
int main()
{
	//problem(0.02, 0.01);
	problem(0.02, 0.03);
}