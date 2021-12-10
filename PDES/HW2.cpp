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
double f_x(double x)
{
	return (PI - x) / 2;
}
double fN_x(double x, int N)
{
	double S = 0;
	for (int i = 1; i <= N; i++)
	{
		S = S + sin(i * x) / i;
	}
	return S;
}
double FN_x(double x, int N)
{
	double S = 0;
	for (int i = 1; i <= N; i++)
	{
		S = S + sin(i * PI / N) / (i * PI / N) * sin(i * x) / i;
	}
	return S;
}
void problem(int M, int N)
{
	double dx = 2 * PI / M;

	long beishu_x = 200, beishu_y = 100;
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
	outtextxy(graphic_x / 10, graphic_x / 10, L"f(x)红色 fN(x)绿色 FN(x)蓝色 f(x)-FN(x)黑色 的图像");

	//setaspectratio(0.005, -1);//这个函数用于设置当前缩放因子。
	setorigin(graphic_x / 2, graphic_y / 2);// 设置坐标原点到屏幕中央点

	POINT* p0 = new POINT[M];
	for (int j = 0; j <= M - 1; j++)
	{
		p0[j].x = beishu_x * j * dx;
		p0[j].y = beishu_y * f_x(j * dx);
	}
	setlinecolor(RED);
	polyline(p0, M);

	POINT* p1 = new POINT[M];
	for (int j = 0; j <= M - 1; j++)
	{
		p1[j].x = beishu_x * j * dx;
		p1[j].y = beishu_y * fN_x(j * dx, N);
	}
	setlinecolor(GREEN);
	polyline(p1, M);

	POINT* p2 = new POINT[M];
	for (int j = 0; j <= M - 1; j++)
	{
		p2[j].x = beishu_x * j * dx;
		p2[j].y = beishu_y * FN_x(j * dx, N);
	}
	setlinecolor(BLUE);
	polyline(p2, M);

	POINT* p3 = new POINT[M];
	for (int j = 0; j <= M - 1; j++)
	{
		p3[j].x = beishu_x * j * dx;
		p3[j].y = beishu_y * (f_x(j * dx) - FN_x(j * dx, N));
	}
	setlinecolor(BLACK);
	polyline(p3, M);

	POINT* p4 = new POINT[M];
	for (int j = 0; j <= M - 1; j++)
	{
		p3[j].x = beishu_x * j * dx;
		p3[j].y = beishu_y * (f_x(j * dx) - fN_x(j * dx, N));
	}
	setlinecolor(LIGHTGREEN);
	polyline(p4, M);

	_getch();
	closegraph();//关闭图形界面

	delete[]p0, p1, p2, p3;
}
int main()
{
	//problem(20, 10);
	//problem(20, 100);
	//problem(160, 10);
	problem(160, 100);
}