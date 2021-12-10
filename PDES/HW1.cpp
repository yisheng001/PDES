#define _CRT_SECURE_NO_WARNINGS//������ȫ���
#include <graphics.h>//ͼ�ν���ͷ�ļ�
#include <conio.h>//��׼�������ͷ�ļ�
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
	int graphic_x = 800, graphic_y = 600;//���ڴ�С

	initgraph(graphic_x, graphic_y);//��ʼ����ͼ����

 //������
	setlinecolor(BLUE);//���û�����ɫ
	line(0, graphic_y / 2, graphic_x, graphic_y / 2);
	line(graphic_x / 2, 0, graphic_x / 2, graphic_y);

	settextcolor(RED);//����������ɫ
	outtextxy(graphic_x - 20, graphic_y / 2 + 5, 'y');
	outtextxy(graphic_x / 2 + 5, 0, 'x');

	settextcolor(LIGHTGREEN);
	outtextxy(graphic_x / 10, graphic_x / 10, L"���� t=0.3 �Ľ��ƽ�(�̣�����ȷ�⣨�죩  X����(0,1) ��ͼ��");

	//setaspectratio(0.005, -1);//��������������õ�ǰ�������ӡ�
	setorigin(graphic_x / 2, graphic_y / 2);// ��������ԭ�㵽��Ļ�����

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
	closegraph();//�ر�ͼ�ν���

	delete[]p0, p1;
}
int main()
{
	//problem(0.02, 0.01);
	problem(0.02, 0.03);
}