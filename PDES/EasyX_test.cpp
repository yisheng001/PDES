//������ƣ�����ͼ��������
//�ļ�����Beauty_Of_Formula.cpp
//�汾��2019_12_20
//���˵�����򵥵Ļ�����ѧ����ͼ�񣬱����û�ͼ���ߣ���ͼ���ڣ���ʵ�ֻ��ƺ���ͼ��
//���ߣ�A������Abr��
//ʱ�䣺2019��12��20��20:04��
#define _CRT_SECURE_NO_WARNINGS//������ȫ���
#include <graphics.h>//ͼ�ν���ͷ�ļ�
#include <conio.h>//��׼�������ͷ�ļ�
#include <math.h>//��ѧ����ͷ�ļ�
#define PI 3.1415926//���ַ����������

int main()//������
{
	int z = 0;
	double* temp = NULL;
	double* hanshu_y = (double*)malloc(sizeof(double) * 100000);//�����ڴ�ռ�

	long beishu_x = 1, beishu_y = 50;
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
	outtextxy(graphic_x / 10, graphic_x / 10, L"��ѧ��ɫ��   ���� y=2sin(2x)  X����(-PI -- PI) ��ͼ��");

	setaspectratio(0.005, -1);//��������������õ�ǰ�������ӡ�
	setorigin(graphic_x / 2, graphic_y / 2);// ��������ԭ�㵽��Ļ�����

	int f = 0;
	//���㺯��ֵ
	for (double x = -PI; x < PI; x += PI / 50000)// -PI~~~~PI
		hanshu_y[z++] = 2 * sin(2 * x);//���㺯��y=2sin(2x)�ģ�-PI~~~~PI���ĺ���ֵ

 //����һ���Ϳ�ʼ������ͼ��
	for (int i = -z / 2; i < z / 2; i++)
	{
		putpixel(beishu_x * i, beishu_y * hanshu_y[f++], LIGHTMAGENTA);
	}

	_getch();
	closegraph();//�ر�ͼ�ν���
	return 0;
}