//软件名称：函数图像生成器
//文件名：Beauty_Of_Formula.cpp
//版本：2019_12_20
//软件说明：简单的画出数学函数图像，本例用绘图工具（绘图窗口）来实现绘制函数图像
//作者：A贝尔（Abr）
//时间：2019年12月20日20:04分
#define _CRT_SECURE_NO_WARNINGS//跳过安全检查
#include <graphics.h>//图形界面头文件
#include <conio.h>//标准输入输出头文件
#include <math.h>//数学运算头文件
#define PI 3.1415926//用字符常量定义π

int main()//主函数
{
	int z = 0;
	double* temp = NULL;
	double* hanshu_y = (double*)malloc(sizeof(double) * 100000);//申请内存空间

	long beishu_x = 1, beishu_y = 50;
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
	outtextxy(graphic_x / 10, graphic_x / 10, L"数学的色彩   函数 y=2sin(2x)  X属于(-PI -- PI) 的图像");

	setaspectratio(0.005, -1);//这个函数用于设置当前缩放因子。
	setorigin(graphic_x / 2, graphic_y / 2);// 设置坐标原点到屏幕中央点

	int f = 0;
	//计算函数值
	for (double x = -PI; x < PI; x += PI / 50000)// -PI~~~~PI
		hanshu_y[z++] = 2 * sin(2 * x);//计算函数y=2sin(2x)的（-PI~~~~PI）的函数值

 //到这一步就开始画函数图像
	for (int i = -z / 2; i < z / 2; i++)
	{
		putpixel(beishu_x * i, beishu_y * hanshu_y[f++], LIGHTMAGENTA);
	}

	_getch();
	closegraph();//关闭图形界面
	return 0;
}