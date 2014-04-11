package com.ethon.plugin;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.lang.StrictMath;
import java.text.*;
import Jama.*;

import com.ethon.DataBase;
import com.ethon.model.DataPoint;
import com.ethon.plugin.LUVtest.Coord;
import com.ethon.tools.ImageGenerator;
import com.ethon.ui.DataField;
import java.math.*;

/**
 * FHS 快速两步图像分割算法
 * 
 * @author liwei
 * 
 */

public class FHS {

	// 图片的行宽和列宽
	int line_num, row_num;

	// 小网格的起点和终点
	public double xmin, ymin, zmin, xmax, ymax, zmax;

	// FHS算法开始运行时间
	private long start;

	// FHS算法运行结束时间
	private long end;

	// 小网格类型的矩阵
	public Cells[][][] cells;

	// 网格系统中所有非空小网格集合
	private ArrayList<Cells> notEmptyCells = new ArrayList<Cells>();

	// 网格系统中非空网格的个数
	private int notEmptyCellsNum = 0;

	// 初始密度函数（每个小网格对应一个初始密度函数）
	public double[][][] initial;

	// 自适应密度函数（每个小网格对应一个自适应密度函数）
	public double[][][] adaptive;

	// 比例因子(用于计算小网格的边长)
	public double[][][] varFactor;

	// 可变带宽（用于计算自适应密度的一个参数值）
	public double[][][] varBandwidth;

	// h(x)的d次方
	public double[][][] varHd;

	// 累乘积（用于计算初始密度函数的几何平均值）
	private BigDecimal multiply;
	
	// 累乘积（用于计算自适应密度函数的几何平均值）
	private BigDecimal multiply2;
	
	// 初始密度函数的几何平均值
	private double geoMean = 0;

	// 自适应性密度函数的几何平均值
	private double geoMean2 = 0;

	// 噪声水平
	private double noiseLevel = 0;

	// 网格系统每边包含小网格的数量
	private int K;

	// 固定带宽h0(由用户输入)
	private int h0;

	// ro因子的初始值是多少？？？
	private double ro = 1.0;

	// 固定窗宽(也就是每个小网格的边长)
	private double width;

	// 所有原始图像数据点
	DataPoint[] RGBPoints;

	// 所有经过转化为LUV坐标系的数据点
	DataPoint[] LUVPoints;

	// 该矩阵中，所有数据点按其在原图像中的坐标(x,y)记录于对应的位置，原图像为line_num*row_num
	DataPoint[][] allPix;

	// 记录所有数据的点
	private DataPoint[] points1 = null;

	// 原始图像数据点的数量
	private int allPointsNum = 0;

	// 初始聚类数目
	private int initalClusterNum;

	// 中间聚类数目
	private int middleClusterNum;
	
	// 最终聚类数目
	private int finalClusterNum;

	// 初始分块数目
	private int initSegmentNum;

	// 最终分块数目
	private int finalSegmentNum;

	// 平滑阈值
	private double smoothValve;

	// 存放平滑后图像中所有的块,也是图像按块划分的最终结果
	LinkedList<Segment> segments = new LinkedList<Segment>();

	// 块所属分组的编号
	private int segmentGroupLabelIndex = 0;

	/**
	 * 构造函数FHS， 接受用户输入参数h0，将所有原始图像数据点映射入网格系统，并确保初始化全部小网格对象 Parameters:
	 * h0:固定带宽（用户输入） line_num:图片的行宽 row_num:图片的列宽 smoothValve:平滑阈值
	 */
	public FHS(int h0, int line_num, int row_num) {

		// 初始化
		start = System.currentTimeMillis();

		this.h0 = h0;
		this.line_num = line_num;
		this.row_num = row_num;

		this.smoothValve = 40;

		allPix = new DataPoint[line_num + 1][row_num + 1];

		// 初始化数据点（原始图像数据点的集合 ）
		RGBPoints = DataBase.getInstance().getPoints();
		points1 = RGBPoints;

		// 公共的数据池
		DataBase db = new DataBase();
		// 保存全局变量points1
		LUVPoints = db.new_public_points(points1);

		// RGB值转换为LUV，结果存到LUVPoint数组中
		LUVPoints = RGBToLUV(LUVPoints);

		// 统计LUV空间像素点的个数（即计算LUVPoints数组长度）
		allPointsNum = LUVPoints.length;

		// 为xmin、xmin、zmin，xmax、ymax、zmax初始化
		xmin = Double.MAX_VALUE;
		ymin = Double.MAX_VALUE;
		zmin = Double.MAX_VALUE;

		xmax = Double.MIN_VALUE;
		ymax = Double.MIN_VALUE;
		zmax = Double.MIN_VALUE;

		// 遍历每个数据点
		for (DataPoint p : LUVPoints) {
			// 数据点在LUV空间的x、y、z坐标值
			double x = p.getCoord_X();
			double y = p.getCoord_Y();
			double z = p.getCoord_Z();
			
			//4月23号新添加内容
			//p.setClusterLabel(-1);
			
			// 计算三维网格系统各维的起点（xmin，ymin，zmin）
			if (xmin > x)
				xmin = x;
			if (ymin > y)
				ymin = y;
			if (zmin > z)
				zmin = z;

			// 计算三维网格系统各维的终点（xmax，ymax，zmax）
			if (xmax < x)
				xmax = x;
			if (ymax < y)
				ymax = y;
			if (zmax < z)
				zmax = z;
		}// foreach循环结束

		// 计算窗宽（也就是小网格的边长）
		this.width = h0 / 3.0;

		System.out.println("小网格的宽度width  = " + width);
		double dx = xmax - xmin;
		double dy = ymax - ymin;
		double dz = zmax - zmin;
		// 比较dx、dy、dz，取其中最大的那个
		double minusMax = 0;
		double temp = dx > dy ? dx :dy;
		minusMax = temp > dz ? temp : dz;
		// 计算网格系统每边包含小网格的数量
		/**
		 * 函数public static double ceil(double a) 返回最小的（最接近负无穷大）double
		 * 值，该值大于等于参数，并等于某个整数。
		 */
		//this.K = (int) Math.ceil((xmax - xmin) / width);
		this.K = (int) Math.ceil(minusMax / width);
		System.out.println(" 每个方向包含的小网格的个数K = " + K);
		
		// 小网格类型矩阵初始化
		cells = new Cells[K][K][K];

		// 将每个数据点分配到相应的小网格中，可能有一部分小网格为空
		for (DataPoint p : LUVPoints) {
			// 问题:每一个数据点的坐标相同？？
			// 数据点的坐标值
			double x = p.getCoord_X();
			double y = p.getCoord_Y();
			double z = p.getCoord_Z();

			/**
			 * public static double floor(double a) 返回最大的（最接近正无穷大）double
			 * 值，该值小于等于参数，并等于某个整数。 cells[xl][yl][zl] 为数据点（x,y,z）对应的网格
			 */
			// 数据点在网格的坐标
			int xl = (int) Math.floor((x - xmin) / width);
			int yl = (int) Math.floor((y - ymin) / width);
			int zl = (int) Math.floor((z - zmin) / width);
			//问题：竟然出现了zl > 11的情况
			System.out.println("cells[" + xl + "][" + yl + "][" + zl +"]");

			if (xl == K)
				xl = K - 1;
			if (yl == K)
				yl = K - 1;
			if (zl == K)
				zl = K - 1;

			// 将数据点添加到小网格矩阵之中
			if (cells[xl][yl][zl] == null)
				//cells[xl][yl][zl] = new Cells();
				cells[xl][yl][zl] = new Cells(xl, yl, zl);

			cells[xl][yl][zl].add(p);

			// 记录所有非空小网格
			/**
			 * public boolean contains(Object o) 如果此列表中包含指定的元素，则返回 true。
			 */

			/**
			 * public boolean add(E e) 将指定的元素添加到此列表的尾部。
			 */
			if (!notEmptyCells.contains(cells[xl][yl][zl]))
				notEmptyCells.add(cells[xl][yl][zl]);

		}// foreach循环结束
		System.out.println("数据点分配到相应小网格完成！");
		
		// 确保所有小网格均被初始化
		// 统计小网格的数量
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					// 小网格的几何中心坐标
					double xloc = xmin + width / 2 + width * x;
					double yloc = ymin + width / 2 + width * y;
					double zloc = zmin + width / 2 + width * z;

					Cells cell = cells[x][y][z];
					// 初始化那些空的小网格
					if (cell == null) {						
						cell = new Cells();
						// 小空网格的中心设为其几何中心
						cell.setCenter(xloc, yloc, zloc);
						cells[x][y][z] = cell;
					}
				}
			}
		}// for循环结束
		System.out.println("小网格初始化完成！");		
		
		// 初始化初始密度函数与自适应密度函数
		initial = new double[K][K][K];
		adaptive = new double[K][K][K];
		varFactor = new double[K][K][K];
		varBandwidth = new double[K][K][K];
		varHd = new double[K][K][K];
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					varFactor[x][y][z] = 0.0;
					initial[x][y][z] = 0.0;
					adaptive[x][y][z] = 0.0;
					varBandwidth[x][y][z] = 0.0;
					varHd[x][y][z] = 0.0;
				}
			}
		}
		//初始化multiply
		multiply = new BigDecimal(Double.toString(1.0));
		System.out.println("初始化的 multiply = " + multiply);
		
		/**
		 * 第3步 初始密度估计
		 */
		// 初始密度估计
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {

					// 初始密度估计,每个小网格有一个初始密度
					// 每个小网格的初始密度 = 每个小网格内包含的数据点个数/LUV空间数据点个数；
					System.out.println("该小网格内数据点的个数："
							+ cells[x][y][z].getpNum());
					System.out.println("LUV颜色空间数据点总个数：" + LUVPoints.length);
					initial[x][y][z] = (double) (cells[x][y][z].getpNum() )/ LUVPoints.length;
					System.out.println("初始密度值：" + initial[x][y][z]);

					// 初始密度累乘
					// 注意：去掉那些空的小网格，只统计那些非空的小网格数量
					//if (!(Math.abs(initial[x][y][z] - 0.0) < 1e-6))
					if((initial[x][y][z] - 0.0) != 0)
					{
						// 每遇到一个初始密度值不为0的网格，非空网格数量就+1
						notEmptyCellsNum++;// 统计非空网格的个数
						System.out.println("非空的网格个数： "+notEmptyCellsNum);						
					}
				}
			}
		}// for循环结束
		
		System.out.println("非空小网格的个数：" + notEmptyCellsNum);
		
		for(int x = 0; x < K; x++){
			for(int y = 0; y < K; y++){
				for(int z = 0; z < K; z++){
					if((initial[x][y][z] - 0.0) != 0){
						
						//double类型必须现转化为String类型，然后通过BigDecimal的构造函数转为BigDecimal类型，进行计算
						/*
						 * BigDecimal bd = new BigDecimal(string);
						 * 详见BigDecimal用法文档
						 */
						double v1 = initial[x][y][z];
						double multiplyDou = 1.0;
						BigDecimal initial = new BigDecimal(Double.toString(v1));
						//先计算每一个初始密度的notEmptyCellsNum次方的值
						multiplyDou = StrictMath.pow(initial.doubleValue(), 1.0/notEmptyCellsNum);
						//累乘计算,得到几何平均值
						multiply = new BigDecimal(Double.toString(multiplyDou));
						multiply = multiply.multiply(initial);
						
						System.out.println("每一次累乘的结果multiply  = " + multiply.doubleValue());
					}//if结束
					
				}
			}
		}//for循环结束	

		System.out.println("multiply  = " + multiply.doubleValue());
		
		// 计算初始密度函数的几何平均值
		geoMean = multiply.doubleValue();
		System.out.println("初始密度函数的集合平均值是：geoMean = " + geoMean);
		
		/**
		 * 第4步 自适应性密度估计
		 */
		// 自适应性密度估计
		//统计个数
		int counters = 1;
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					counters++;
					double plus = 0;// 累加和

					//只计算当前小网格周围的26个小网格的情况
					double Z1 = disLeftCell(x, y, z);
					System.out.println("Z1 = " + Z1);
					
					double Z2 = disRightCell(x, y, z);
					System.out.println("Z2 = " + Z2);
					
					double Z3 = disUpCell(x, y, z);
					System.out.println("Z3 = " + Z3);
					
					double Z4 = disDownCell(x, y, z);
					System.out.println("Z4 = " + Z4);
					
					double Z5 = disFrontCell(x, y, z);
					System.out.println("Z5 = " + Z5);
					
					double Z6 = disBackCell(x, y, z);
					System.out.println("Z6 = " + Z6);
					
					double Z7 = disLeftUpCell(x, y, z);
					System.out.println("Z7 = " + Z7);
					
					double Z8 = disLeftDownCell(x, y, z);
					System.out.println("Z8 = " + Z8);
					
					double Z9 = disLeftFrontCell(x, y, z);
					System.out.println("Z9 = " + Z9);
					
					double Z10 = disLeftBackCell(x, y, z);
					System.out.println("Z10 = " + Z10);
					
					double Z11= disRightUpCell(x, y, z);
					System.out.println("Z11 = " + Z11);
					
					double Z12 = disRightDownCell(x, y, z);
					System.out.println("Z12 = " + Z12);
					
					double Z13 = disRightFrontCell(x, y, z);
					System.out.println("Z13 = " + Z13);
					
					double Z14 = disRightBackCell(x, y, z);
					System.out.println("Z14 = " + Z14);
					
					double Z15 = disUpFrontCell(x, y, z);
					System.out.println("Z15 = " + Z15);
					
					double Z16 = disUpBackCell(x, y, z);
					System.out.println("Z16 = " + Z16);
					
					double Z17 = disDownFrontCell(x, y, z);
					System.out.println("Z17 = " + Z17);
					
					double Z18 = disDownBackCell(x, y, z);
					System.out.println("Z18 = " + Z18);
					
					double Z19 = disLeftUpFrontCell(x, y, z);
					System.out.println("Z19 = " + Z19);
					
					double Z20 = disLeftUpBackCell(x, y, z);
					System.out.println("Z20 = " + Z20);
					
					double Z21 = disLeftDownFrontCell(x, y, z);
					System.out.println("Z21 = " + Z21);
					
					double Z22 = disLeftDownBackCell(x, y, z);
					System.out.println("Z22 = " + Z22);
					
					double Z23 = disRightUpFrontCell(x, y, z);
					System.out.println("Z23 = " + Z23);
					
					double Z24 = disRightUpBackCell(x, y, z);
					System.out.println("Z24 = " + Z24);
					
					double Z25 = disRightDownFrontCell(x, y, z);
					System.out.println("Z25 = " + Z25);
					
					double Z26 = disRightDownBackCell(x, y, z);
					System.out.println("Z26 = " + Z26);
					
					plus = Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10 + Z11 + Z12 + Z13 + Z14 +Z15 +Z16
							+ Z17 + Z18 + Z19 + Z20 + Z21 + Z22 + Z23 + Z24 + Z25 + Z26;
					System.out.println("plus = " + plus);
					
					// 计算以x为基准点的x的自适应性密度函数值(有问题：是这个位置是所有网格的个数还是非空网格的个数？？)					
					adaptive[x][y][z] = plus / notEmptyCellsNum;
					System.out.println("自适应密度函数: " + adaptive[x][y][z]);
				}
			}
		}// 外层的for循环结束

		
		//初始化multiply2
		multiply2 = new BigDecimal(Double.toString(1.0));
		System.out.println("初始化的 multiply2 = " + multiply2);
		
		// 计算自适应密度函数值
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					//只统计非空网格
					if ((adaptive[x][y][z] - 0.0) != 0) {
						
						double v2 = adaptive[x][y][z];
						
						BigDecimal adaptive = new BigDecimal(Double.toString(v2));
						
						double multiplyDou = 1.0;
						//先计算每一个初始密度的notEmptyCellsNum次方的值				
						multiplyDou = StrictMath.pow(adaptive.doubleValue(), 1.0 / notEmptyCellsNum);
						
						//累乘计算,得到几何平均值
						multiply2 = new BigDecimal(Double.toString(multiplyDou));
						multiply2 = multiply2.multiply(adaptive);
						
						System.out.println("每一次累乘的结果multiply2  = " + multiply2.doubleValue());												
					}// if结束
				}
			}
		}// for循环结束

		// 自适应密度函数的几何平均值
		geoMean2 = multiply2.doubleValue();
		System.out.println("自适应密度函数的几何平均值 :" + geoMean2);

	}// 构造函数结束
	
	
	

	/**
	 * 类Cells
	 * 定义空间小网格对象
	 * 
	 */
	public class Cells {

		// X，Y，Z轴坐标值
		private int x;

		private int y;

		private int z;

		// 各坐标轴的数值统计，累加和
		private double xall = 0;

		private double yall = 0;

		private double zall = 0;

		// 对应的图像的坐标值累加和
		private double spacialXall = 0;

		private double spacialYall = 0;

		// 小网格对象的x、y、z轴的坐标值（网格内数据点坐标的平均值）
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// 对应的图像的坐标平均值(数据点对应的图像的x、y坐标值的算数平均值)
		private double spacialX = -1;

		private double spacialY = -1;

		// x、y、z 三个方向的三个偏导势
		private double xdp = -1;

		private double ydp = -1;

		private double zdp = -1;

		// 小网格中包含的数据点的个数
		private int pNum;

		// 单个网格所属的聚类编号
		private int clusterNum;

		// 小网格的初始密度值
		private double density1 = -1;

		// 小网格的自适应密度值
		private double density2 = -1;

		// 小网格的势值
		private double potential = -1;

		// 小网格包含的所有数据点集
		private ArrayList<DataPoint> list = new ArrayList<DataPoint>();

		// 小网格是否被搜索过
		private int alreadysearched = 0;

		// 小网格是否被合并过
		private int alreadyMerged = 0;
		
		// 带参数的构造函数
		public Cells(int x, int y, int z) {
			this.setX(x);
			this.setY(y);
			this.setZ(z);
		}

		// 不带参数的构造函数
		public Cells() {
			
		}

		// 把数据点添加入相应的小网格之中
		public void add(DataPoint p) {
			// 相应的小网格内数据点的个数+1
			pNum++;

			xall += p.getCoord_X();
			yall += p.getCoord_Y();
			zall += p.getCoord_Z();

			// 数据点在图片中的坐标(x,y)=[line,row]
			spacialXall += p.getLine();
			spacialYall += p.getRow();

			// 把数据点加入到网格中
			list.add(p);
		}

		
		// 设置小网格对象的坐标(xloc, yloc, zloc)
		public void setCenter(double xloc, double yloc, double zloc) {
			this.xloc = xloc;
			this.yloc = yloc;
			this.zloc = zloc;
		}

		// getter and setter 方法
		
		public double getDensity1() {
			return density1;
		}

		public int getAlreadyMerged() {
			return alreadyMerged;
		}

		public void setAlreadyMerged(int alreadyMerged) {
			this.alreadyMerged = alreadyMerged;
		}

		public double getPotential() {
			return potential;
		}

		public void setPotential(double potential) {
			this.potential = potential;
		}

		public void setDensity1(double density1) {
			this.density1 = density1;
		}

		public double getDensity2() {
			return density2;
		}

		public void setDensity2(double density2) {
			this.density2 = density2;
		}

		public int getX() {
			return x;
		}

		public void setX(int x) {
			this.x = x;
		}

		public int getY() {
			return y;
		}

		public void setY(int y) {
			this.y = y;
		}

		public int getZ() {
			return z;
		}

		public void setZ(int z) {
			this.z = z;
		}

		public double getXall() {
			return xall;
		}

		public void setXall(double xall) {
			this.xall = xall;
		}

		public double getYall() {
			return yall;
		}

		public void setYall(double yall) {
			this.yall = yall;
		}

		public double getZall() {
			return zall;
		}

		public void setZall(double zall) {
			this.zall = zall;
		}

		public double getSpacialXall() {
			return spacialXall;
		}

		public void setSpacialXall(double spacialXall) {
			this.spacialXall = spacialXall;
		}

		public double getSpacialYall() {
			return spacialYall;
		}

		public void setSpacialYall(double spacialYall) {
			this.spacialYall = spacialYall;
		}

		//小网格对象的坐标值,取几何中心的坐标值作为该对象的坐标值
		public double getXloc() {

			if(xloc == -1){
				if(this.pNum != 0){
					xloc = this.xall / this.pNum;
				}
			}			
			return xloc;
		}

		public void setXloc(double xloc) {
			this.xloc = xloc;
		}

		public double getYloc() {
			if(yloc == -1){
				if(this.pNum != 0){
					yloc = this.yall / this.pNum;
				}
			}
			
			return yloc;
		}

		public void setYloc(double yloc) {
			this.yloc = yloc;
		}

		public double getZloc() {
			if(zloc == -1){
				if(this.pNum != 0){
					zloc = this.zall / this.pNum;
				}
			}
			
			return zloc;
		}

		public void setZloc(double zloc) {
			this.zloc = zloc;
		}

		public double getSpacialX() {
			if (spacialX == -1) {
				if (pNum != 0)
					spacialX = spacialXall / pNum;
				else
					spacialX = 0;
			}
			return spacialX;
		}

		public void setSpacialX(double spacialX) {
			this.spacialX = spacialX;
		}

		public double getSpacialY() {
			if (spacialY == -1) {
				if (pNum != 0)
					spacialY = spacialYall / pNum;
				else
					spacialY = 0;
			}
			return spacialY;
		}

		public void setSpacialY(double spacialY) {
			this.spacialY = spacialY;
		}

		public int getpNum() {
			return pNum;
		}

		public void setpNum(int pNum) {
			this.pNum = pNum;
		}

		public int getClusterNum() {
			return clusterNum;
		}

		public void setClusterNum(int clusterNum) {
			this.clusterNum = clusterNum;
		}

		public ArrayList<DataPoint> getList() {
			return list;
		}

		public void setList(ArrayList<DataPoint> list) {
			this.list = list;
		}

		public double getXdp() {
			return xdp;
		}

		public void setXdp(double xdp) {
			this.xdp = xdp;
		}

		public double getYdp() {
			return ydp;
		}

		public void setYdp(double ydp) {
			this.ydp = ydp;
		}

		public double getZdp() {
			return zdp;
		}

		public void setZdp(double zdp) {
			this.zdp = zdp;
		}

		public int getAlreadysearched() {
			return alreadysearched;
		}

		public void setAlreadysearched(int alreadysearched) {
			this.alreadysearched = alreadysearched;
		}

		// 以下方法不知道有什么用
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result + getX();
			result = prime * result + getY();
			result = prime * result + getZ();
			return result;
		}

		@Override
		// 判断两个网格坐标是否相等
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Cells other = (Cells) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (getX() != other.getX())
				return false;
			if (getY() != other.getY())
				return false;
			if (getZ() != other.getZ())
				return false;
			return true;
		}

		public String toString() {
			String info = new String();
			info = "cells [" + this.getX() + "]" + "[" + this.getY() + "]"
					+ "[" + this.getZ() + "]";
			return info;
		}

		private FHS getOuterType() {
			return FHS.this;
		}

	}

	/**
	 * 类Cluster
	 * 定义网格聚类对象 一个聚类代表若干个小网格
	 * 
	 */
	public class Cluster {
		// 聚类对象的编号
		private int clusterLabel;

		// 聚类对象的x、y、z三轴坐标值
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// 不同聚类对象内的多个像素点的z、y、z坐标轴的累加和
		private double xAll = 0;

		private double yAll = 0;

		private double zAll = 0;

		// 聚类范围域内包含的数据点总个数
		private int pointsAll = 0;

		// 聚类范围域网格坐标集
		LinkedList<Cells> contentList;

		// 聚类几何中心网格（该小网格坐标值作为该聚类的坐标值）
		Cells clusterCenter;

		// 初始聚类中心小网格
		Cells initialClusterCenter;

		//多余的变量
		private double distanceSQ;
		
		// 以下是三个构造函数

		// 寻找初始聚类中心小网格周围的六邻域（相邻的六个小网格）
		public Cluster(int x, int y, int z, int clusterLabel) {
			// 初始聚类中心的小网格是cell[x][y][z]
			this.initialClusterCenter = new Cells(x, y, z);

			// 聚类对象的编号
			this.clusterLabel = clusterLabel;

			// 初始聚类范围域内小网格坐标集
			contentList = new LinkedList<Cells>();

			// 将初始聚类中心小网格添加到初始聚类范围域
			// 初始聚类范围域只有一个元素（即一个作为基准的小网格）
			contentList.add(initialClusterCenter);

		}// 构造函数结束

		//无参数构造函数
		public Cluster() {

		}
		
		//带一个参数的构造函数
		public Cluster(int clusterLabel) {
			this.clusterLabel = clusterLabel;
			contentList = new LinkedList<Cells>();
		}
		
		public Cluster(int x, int y, int z) {
			// 初始聚类中心的小网格是cell[x][y][z]
			this.initialClusterCenter = new Cells(x, y, z);

			// 初始聚类范围域内小网格坐标集
			contentList = new LinkedList<Cells>();

			// 将初始聚类中心小网格添加到初始聚类范围域
			// 初始聚类范围域只有一个元素（即一个作为基准的小网格）
			contentList.add(initialClusterCenter);

		}// 构造函数结束
		
		
		
		// getter and setter方法
		
		public Cells getClusterCenter() {
			return clusterCenter;
		}

		public double getDistanceSQ() {
			return distanceSQ;
		}

		//已经作出修改
		public void setDistanceSQ(Cluster c) {
			this.distanceSQ = computeDistance(c);
		}

		public void setClusterCenter(Cells clusterCenter) {
			this.clusterCenter = clusterCenter;
		}

		public void setPointsAll(int pointsAll) {
			this.pointsAll = pointsAll;
		}

		public int getClusterLabel() {
			return clusterLabel;
		}

		public void setClusterLabel(int clusterLabel) {
			this.clusterLabel = clusterLabel;
		}

		public double getXloc() {
			return xloc;
		}

		public void setXloc(double xloc) {
			this.xloc = xloc;
		}

		public double getYloc() {
			return yloc;
		}

		public void setYloc(double yloc) {
			this.yloc = yloc;
		}

		public double getZloc() {
			return zloc;
		}

		public void setZloc(double zloc) {
			this.zloc = zloc;
		}

		// 统计聚类范围域内包含的所有小网格内数据点的x坐标值的累加和
		public double getxAll() {
			xAll = 0;
			for (Cells c : contentList) {
				xAll += cells[c.getX()][c.getY()][c.getZ()].getXall();
			}
			System.out.println("xAll = " + xAll);
			return xAll;			
		}

		public void setxAll(double xAll) {
			this.xAll = xAll;
		}

		// 统计聚类范围域内包含的所有小网格内数据点的y坐标值的累加和
		public double getyAll() {
			yAll = 0;
			for (Cells c : contentList) {
				yAll += cells[c.getX()][c.getY()][c.getZ()].getYall();
			}
			System.out.println("yAll = " + yAll);
			return yAll;
		}

		public void setyAll(double yAll) {
			this.yAll = yAll;
		}

		// 统计聚类范围域内包含的所有小网格内数据点的z坐标值的累加和
		public double getzAll() {
			zAll = 0;
			for (Cells c : contentList) {
				zAll += cells[c.getX()][c.getY()][c.getZ()].getZall();
			}
			System.out.println("zAll = " + zAll);
			return zAll;
		}

		public void setzAll(double zAll) {
			this.zAll = zAll;
		}

		public int getPointsAll() {
			return pointsAll;
		}

		// setPointsAll()方法做修改，统计聚类范围域内数据点的个数（不是小网格的个数）
		public void setPointsAll() {
			for (Cells c : contentList) {
				// 首先寻找到聚类范围域的每个小网格，通过小网格对象的getpNum()方法获得每个小网格内的数据点的个数
				pointsAll += cells[c.getX()][c.getY()][c.getZ()].getpNum();
			}
		}		

		//聚类范围域内所有的小网格
		public LinkedList<Cells> getContentList() {
			return contentList;
		}

		public void setContentList(LinkedList<Cells> contentList) {
			this.contentList = contentList;
		}

		public Cells getInitialClusterCenter() {
			return initialClusterCenter;
		}

		public void setInitialClusterCenter(Cells initialClusterCenter) {
			this.initialClusterCenter = initialClusterCenter;
		}

		// 聚类之后计算聚类范围域的重点
		public void adjustClusterCenter() {
			this.setPointsAll();			
			
			/*
			 * public static double floor(double a)
			 * 返回最大的（最接近正无穷大）double 值，该值小于等于参数，并等于某个整数
			 * return :最大（最接近正无穷大）浮点值，该值小于等于该参数，并等于某个整数
			 */
			System.out.println("this.getxAll() = " + this.getxAll());
			System.out.println("this.getPointsAll() = " + this.getPointsAll());
			System.out.println("this.getyAll() = " + this.getyAll());
			System.out.println("this.getzAll() = " + this.getzAll());
			System.out.println("xmin = " + xmin);
			System.out.println("width = " + width);
			
			int x = (int) Math
					.floor((this.getxAll() / this.getPointsAll() - xmin)
							/ width);
			if (x == K)
				x = K - 1;

			int y = (int) Math
					.floor((this.getyAll() / this.getPointsAll() - ymin)
							/ width);
			if (y == K)
				y = K - 1;

			int z = (int) Math
					.floor((this.getzAll() / this.getPointsAll() - zmin)
							/ width);
			if (z == K)
				z = K - 1;
			
			this.setClusterCenter(new Cells(x, y, z));
		}
		
		//多余的函数
		public double computeDistance(Cluster c) {

			distanceSQ = (this.getxAll() / this.getPointsAll() - c.getxAll()
					/ c.getPointsAll())
					* (this.getxAll() / this.getPointsAll() - c.getxAll()
							/ c.getPointsAll())
					+ (this.getyAll() / this.getPointsAll() - c.getyAll()
							/ c.getPointsAll())
					* (this.getyAll() / this.getPointsAll() - c.getyAll()
							/ c.getPointsAll())
					+ (this.getzAll() / this.getPointsAll() - c.getzAll()
							/ c.getPointsAll())
					* (this.getzAll() / this.getPointsAll() - c.getzAll()
							/ c.getPointsAll());
			return distanceSQ;
			//this.setDistanceSQ(distanceSQ);
		}
		
		//当参数为单独一个小网格的时候,如何把该小网格添加到聚类范围域内
		public void addCell(Cells c){
			//如果小网格为空，添加到聚类范围域（也就是说，即使是空的小网格也属于某个聚类范围域内）
			if(xloc == -1 && yloc == -1 && zloc ==-1){
				xloc = cells[c.getX()][c.getY()][c.getZ()].getXloc();
				yloc = cells[c.getX()][c.getY()][c.getZ()].getYloc();
				zloc = cells[c.getX()][c.getY()][c.getZ()].getZloc();
			}
			contentList.add(c);
		}
		
		//当参数为多个小网格的时候，如何把该小网格添加到聚类范围域内
		public void addCells(LinkedList<Cells> cs){
			for(Cells c : cs){
				contentList.add(c);//一次把每一个小网格添加到聚类范围域网格系统中
			}
		}
		
	}

	/**
	 * 类Segment
	 * 定义图像分块对象 每一块代表若干个数据点
	 * 
	 */
	private class Segment {
		// 块编号
		private int segmentLabel;

		// 该块所在组的编号
		private int groupLabel = -1;

		// 块中心特征空间三轴坐标
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// 特征空间各坐标值统计（累加和）
		private double xAll = 0;

		private double yAll = 0;

		private double zAll = 0;;

		// 块中包含的数据点总数
		private double numAll = 0;

		// 块内数据点集合
		LinkedList<DataPoint> contentList;

		// 块内边界点集合
		LinkedList<DataPoint> marginList;

		// 相邻块的编号集合
		LinkedList<Integer> neighborSegLabelList;

		//getter and setter方法
		
		public void setSegmentLabel(int segmentLabel) {
			this.segmentLabel = segmentLabel;
		}

		public void setXloc(double xloc) {
			this.xloc = xloc;
		}

		public void setYloc(double yloc) {
			this.yloc = yloc;
		}

		public void setZloc(double zloc) {
			this.zloc = zloc;
		}

		public void setxAll(double xAll) {
			this.xAll = xAll;
		}

		public void setyAll(double yAll) {
			this.yAll = yAll;
		}

		public void setzAll(double zAll) {
			this.zAll = zAll;
		}

		public void setNumAll(double numAll) {
			this.numAll = numAll;
		}

		public int getGroupLabel() {
			return groupLabel;
		}

		public void setGroupLabel(int groupLabel) {
			this.groupLabel = groupLabel;
		}

		public double getXloc() {
			xloc = this.getxAll() / this.getNumAll();
			return xloc;
		}

		public double getYloc() {
			yloc = this.getyAll() / this.getNumAll();
			return yloc;
		}

		public double getZloc() {
			zloc = this.getzAll() / this.getNumAll();
			return zloc;
		}

		public double getxAll() {
			xAll = 0;
			for (DataPoint p : contentList) {
				xAll += p.getCoord_X();
			}
			return xAll;
		}

		public double getyAll() {
			yAll = 0;
			for (DataPoint p : contentList) {
				yAll += p.getCoord_Y();
			}
			return yAll;
		}

		public double getzAll() {
			zAll = 0;
			for (DataPoint p : contentList) {
				zAll += p.getCoord_Z();
			}
			return zAll;
		}

		public double getNumAll() {
			return numAll;
		}

		public void setNumAll() {
			this.numAll = this.contentList.size();
		}

		public LinkedList<DataPoint> getContentList() {
			return contentList;
		}

		public LinkedList<DataPoint> getMarginList() {
			return marginList;
		}

		public LinkedList<Integer> getNeighborSegLabelList() {
			return neighborSegLabelList;
		}

		public int getSegmentLabel() {
			return segmentLabel;
		}

		public void setContentList(LinkedList<DataPoint> contentList) {
			this.contentList = contentList;
		}

		public void setMarginList(LinkedList<DataPoint> marginList) {
			this.marginList = marginList;
		}

		public void setNeighborSegLabelList(
				LinkedList<Integer> neighborSegLabelList) {
			this.neighborSegLabelList = neighborSegLabelList;
		}

		// 带参数的构造函数
		public Segment(int segmentLabel, int groupLabel) {
			this.segmentLabel = segmentLabel;
			this.groupLabel = groupLabel;
			contentList = new LinkedList<DataPoint>();
			marginList = new LinkedList<DataPoint>();
			neighborSegLabelList = new LinkedList<Integer>();
		}
		// 不参与分组的构造函数
		public Segment(int segmentLabel) {
			this.segmentLabel = segmentLabel;
			contentList = new LinkedList<DataPoint>();
			marginList = new LinkedList<DataPoint>();
			neighborSegLabelList = new LinkedList<Integer>();
		}

	}

	/**
	 * 类Group
	 * 定义图像组对象,及块的的集合 每一组包含很多的图像快
	 * 
	 */
	private class Group {
		// 组编号
		private int groupLabel;

		// 住中心特征空间三轴坐标
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// 特征空间各坐标值统计
		private double xAll = 0;

		private double yAll = 0;

		private double zAll = 0;

		// 组中包含的数据点总数
		private double numAll = 0;

		// 组内数据点集合
		LinkedList<DataPoint> contentList;

		//getter and setter 方法
		
		public int getGroupLabel() {
			return groupLabel;
		}

		public double getXloc() {
			xloc = this.getxAll() / this.getNumAll();
			return xloc;
		}

		public double getYloc() {
			yloc = this.getyAll() / this.getNumAll();
			return yloc;
		}

		public double getZloc() {
			zloc = this.getzAll() / this.getNumAll();
			return zloc;
		}

		public double getxAll() {
			xAll = 0;
			for (DataPoint p : contentList) {
				xAll += p.getCoord_X();
			}
			return xAll;
		}

		public double getyAll() {
			yAll = 0;
			for (DataPoint p : contentList) {
				yAll += p.getCoord_Y();
			}
			return yAll;
		}

		public double getzAll() {
			zAll = 0;
			for (DataPoint p : contentList) {
				zAll += p.getCoord_Z();
			}
			return zAll;
		}

		public double getNumAll() {
			return numAll;
		}

		public void setNumAll() {
			this.numAll = this.contentList.size();
		}

		public LinkedList<DataPoint> getContentList() {
			return contentList;
		}

		public void setContentList(LinkedList<DataPoint> contentList) {
			this.contentList = contentList;
		}

		// 带参数的构造函数
		public Group(int groupLabel) {
			this.groupLabel = groupLabel;
			contentList = new LinkedList<DataPoint>();
		}
		
		
		//以下方法貌似没用到
		public void setGroupLabel(int groupLabel) {
			this.groupLabel = groupLabel;
		}

		public void setXloc(double xloc) {
			this.xloc = xloc;
		}

		public void setYloc(double yloc) {
			this.yloc = yloc;
		}

		public void setZloc(double zloc) {
			this.zloc = zloc;
		}

		public void setxAll(double xAll) {
			this.xAll = xAll;
		}

		public void setyAll(double yAll) {
			this.yAll = yAll;
		}

		public void setzAll(double zAll) {
			this.zAll = zAll;
		}

		public void setNumAll(double numAll) {
			this.numAll = numAll;
		}

	}

	
	
	/**
	 * 函数RGBToLUV，
	 * 接收数组RGBPoints,将所有原始图像数据点的特征三维坐标值(R,G,B)按公式转换为(L,U,V),用于聚类和分割操作
	 * Parameters: RGBPoints 所有RGB数据点 
	 * return DataPoint[]
	 */
	public DataPoint[] RGBToLUV(DataPoint[] RGBPoints) {

		double L = 0;
		double U = 0;
		double V = 0;

		// 二维数据存放转换矩阵
		double[][] tansformMatrixRGBToLUV = new double[][] {
				{ 0.412453, 0.357580, 0.180423 },
				{ 0.212671, 0.715160, 0.072169 },
				{ 0.019334, 0.119193, 0.950227 } };

		double Xn = 0.950456;
		double Yn = 1;
		double Zn = 1.088754;

		double un = 4 * Xn / (Xn + 15 * Yn + 3 * Zn);
		double vn = 9 * Yn / (Xn + 15 * Yn + 3 * Zn);

		for (int i = 0; i < RGBPoints.length; i++) {
			double LUVx = RGBPoints[i].getCoord_X() / 255
					* tansformMatrixRGBToLUV[0][0] + RGBPoints[i].getCoord_Y()
					/ 255 * tansformMatrixRGBToLUV[0][1]
					+ RGBPoints[i].getCoord_Z() / 255
					* tansformMatrixRGBToLUV[0][2];
			double LUVy = RGBPoints[i].getCoord_X() / 255
					* tansformMatrixRGBToLUV[1][0] + RGBPoints[i].getCoord_Y()
					/ 255 * tansformMatrixRGBToLUV[1][1]
					+ RGBPoints[i].getCoord_Z() / 255
					* tansformMatrixRGBToLUV[1][2];
			double LUVz = RGBPoints[i].getCoord_X() / 255
					* tansformMatrixRGBToLUV[2][0] + RGBPoints[i].getCoord_Y()
					/ 255 * tansformMatrixRGBToLUV[2][1]
					+ RGBPoints[i].getCoord_Z() / 255
					* tansformMatrixRGBToLUV[2][2];

			double u = 0;
			double v = 0;
			if (LUVx + 15 * LUVy + 3 * LUVz != 0) {
				u = 4 * LUVx / (LUVx + 15 * LUVy + 3 * LUVz);
				v = 9 * LUVy / (LUVx + 15 * LUVy + 3 * LUVz);
			}

			if (LUVy / Yn > 0.008856) {
				L = 116 * StrictMath.pow(LUVy / Yn, 1.0 / 3) - 16;
			} else {
				L = 116 * (7.787 * LUVy / Yn + 16.0 / 116) - 16;
			}

			U = 13 * L * (u - un);
			V = 13 * L * (v - vn);
			// 将新得到的LUV值存到原像素点对象
			RGBPoints[i].setCoord_X(L);
			RGBPoints[i].setCoord_Y(U);
			RGBPoints[i].setCoord_Z(V);

		}
		return RGBPoints;
	}

	
	/**
	 * 函数LUVToRGB，接收数组LUVPoints,将处理所得的图像数据点的特征三维坐标值(L,U,V)转换为(R,G,B),用于图像绘制
	 * Parameters: LUVPoints 所有LUV数据点 
	 * return DataPoint[]
	 */
	public DataPoint[] LUVToRGB(DataPoint[] LUVPoints) {

		int R = 0;
		int G = 0;
		int B = 0;

		double[][] tansformMatrixLUVToRGB = new double[][] {
				{ 3.240479, -1.537150, -0.498535 },
				{ -0.969256, 1.875992, 0.041556 },
				{ 0.055648, -0.204043, 1.057311 } };

		double Xn = 0.950456;
		double Yn = 1;
		double Zn = 1.088754;

		double un = 4 * Xn / (Xn + 15 * Yn + 3 * Zn);
		double vn = 9 * Yn / (Xn + 15 * Yn + 3 * Zn);

		for (int i = 0; i < LUVPoints.length; i++) {
			double u = LUVPoints[i].getCoord_Y()
					/ (13 * LUVPoints[i].getCoord_X()) + un;

			double v = LUVPoints[i].getCoord_Z()
					/ (13 * LUVPoints[i].getCoord_X()) + vn;

			double RGBy = 0;
			if ((LUVPoints[i].getCoord_X() + 16) / 116 > 0.206893) {
				RGBy = StrictMath.pow((LUVPoints[i].getCoord_X() + 16) / 116,
						3.0);
			} else {
				RGBy = ((LUVPoints[i].getCoord_X() + 16) / 116 - 16.0 / 116) / 7.787;
			}

			double RGBx = 0;

			double RGBz = 0;

			if (v != 0) {
				RGBx = 9 * u / (4 * v) * RGBy;
				RGBz = (3 - 0.75 * u - 5 * v) / v * RGBy;
			}

			R = (int) ((tansformMatrixLUVToRGB[0][0] * RGBx
					+ tansformMatrixLUVToRGB[0][1] * RGBy + tansformMatrixLUVToRGB[0][2]
					* RGBz) * 255);

			G = (int) ((tansformMatrixLUVToRGB[1][0] * RGBx
					+ tansformMatrixLUVToRGB[1][1] * RGBy + tansformMatrixLUVToRGB[1][2]
					* RGBz) * 255);

			B = (int) ((tansformMatrixLUVToRGB[2][0] * RGBx
					+ tansformMatrixLUVToRGB[2][1] * RGBy + tansformMatrixLUVToRGB[2][2]
					* RGBz) * 255);
			// 将新得到的(R,G,B)值存到原LUVPoints[]像素点对象
			LUVPoints[i].setCoord_X(R);
			LUVPoints[i].setCoord_Y(G);
			LUVPoints[i].setCoord_Z(B);

		}

		return LUVPoints;
	}

	
	/**
	 * 范围域聚类逻辑处理函数process 
	 * 1、初始密度估计； 
	 * 2、自适应性密度估计；
	 * 3、范围域聚类 ；
	 * 4、将范围域聚类结果映射到图像空间；
	 */
	public void process() {		
		/**
		 * 第5步 范围域聚类
		 * 1、定义噪声水平noiseLevel;
		 * 2、确定初始聚类中心小网格。确定某个基准小网格，比较与周围六个小网格梯度的差值，并将梯度差值最大的下一个网格作为基准网格，继续判断
		 * 直到相邻网格的梯度差值最大的那个<0.01，那么这个中心小网格就是某个初始聚类中心initialClusterCenter；
		 * 3、依据初始聚类中心，搜索各中心完整的网格域，并按聚类结果设置范围域内的所有数据点的类编号；
		 * 4、判断初始聚类中心小网格的自适应性密度函数adaptive[x][y][z] < noiseLevel的话，则该聚类结果归结为噪声；
		 */

		double constant1 = 0.1;
		noiseLevel = constant1 * geoMean2;
		System.out.println("噪声水平 noiseLevel = " + noiseLevel);
		
		/**
		 * 过程：
		 * 1、首先确定某个小网格作为基准小网格； 问题：如果用三个for循环遍历每一个小网格的话，怎么
		 * 2、由基准小网格寻找下一个next小网格；调用searchNextCell(Cells c)；
		 * 3、把next小网格作为新的基准小网格继续寻找下一个next小网格；
		 * 4、直到next小网格的Max < 0.01，寻找结束，当前小网格作为聚类中心center小网格；
		 * 5、找到以center为中心的所有的小网格path，把找到的所有小网格作为一个聚类结果；
		 * 6、为这个聚类结果编写类编号。
		 */
		
		// 计算各小网格的x偏导，y偏导，z偏导
		partialDerivative();

		//System.out.println("计算偏导结束！");
		
		// 存储初始聚类中心小网格集合(每一个聚类对应一个中心小网格)
		ArrayList<Cells> centers = new ArrayList<Cells>();
		
		//初始聚类集合
		ArrayList<Cluster> initialClusters = new ArrayList<Cluster>();
		
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					
					// alreadysearched == 0 表示没有被搜索过。
					if(cells[x][y][z].getAlreadysearched() == 0){						
						
						//新建一个聚类对象
						Cluster cluster = new Cluster(-1, -1, -1);						
						
						// 聚类中心
						Cells center = new Cells(-1, -1, -1);

						// 存储路径上的所有的小网格对象
						LinkedList<Cells> path = new LinkedList<Cells>();
						//LinkedList<Cells> contentList = new LinkedList<Cells>();
						//System.out.println("打印路径长度(添加自身之前): path.size() = " + path.size());
						
						// 加入自身
						path.add(new Cells(x, y, z));

						//System.out.println("打印路径长度（添加自身之后）: path.size() = " + path.size());
						System.out.println("cells[" + path.get(0).getX() + "][" + path.get(0).getY() + "][" + path.get(0).getZ() +"]");
						
						// 旧的基准点初始化(带参数的构造函数)
						Cells oldBase = new Cells(-1, -1, -1);
						//oldBase = null;
						
						// 新的基准点初始化
						Cells newBase =new Cells(-1, -1, -1);
						//newBase = null;
						
						// 计数器
						int counter = 0;
						
						//因子
						double peak = 1.0;
						
						//初始化数组
						double[] v = { -1, -1, -1, -1, -1, -1 };	
						do {							
							oldBase = path.get(counter);

							// 基准小网格对象的坐标值
							int x1 = oldBase.getX();
							int y1 = oldBase.getY();
							int z1 = oldBase.getZ();

							System.out.println("基准点为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							//基准点进来就先标记为1
							cells[x1][y1][z1].setAlreadysearched(1);
							
							v[0] = gapLeftCell(x1, y1, z1);
							v[1] = gapRightCell(x1, y1, z1);
							v[2] = gapUpCell(x1, y1, z1);
							v[3] = gapDownCell(x1, y1, z1);
							v[4] = gapFrontCell(x1, y1, z1);
							v[5] = gapBackCell(x1, y1, z1);
							
							//依次打印出来数组值
							for(int i = 0; i < v.length; i++){
								System.out.println("v[" + i +"] = "+ v[i]);
							}
							
							// v[0]最大
							if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
									&& v[0] >= v[4] && v[0] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[0] != -1){
									if(v[0] < peak){										
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
									} else {										
										newBase = new Cells(x1 - 1, y1, z1);
										System.out.println("新的基准点为:cells["+ (x1 - 1) +"][" + y1 + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1 - 1][y1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}// if结束

							
							// v[1]最大
							if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
									&& v[1] >= v[4] && v[1] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[1] != -1){
									if(v[1] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//如何把一个对象赋给另一个对象
										//newBase = cells[x1 + 1][y1][z1];
										newBase = new Cells(x1 + 1, y1, z1);
										System.out.println("新的基准点为:cells["+ (x1 + 1) +"][" + y1 + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1 + 1][y1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[2]最大
							if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
									&& v[2] >= v[4] && v[2] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[2] != -1){
									if(v[2] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1 - 1][z1];
										newBase = new Cells(x1, y1 - 1, z1);
										System.out.println("新的基准点为:cells["+ x1 +"][" + (y1 - 1) + "][" + z1 + "]");
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1][y1 - 1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}							
							

							// v[3]最大
							if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
									&& v[3] >= v[4] && v[3] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[3] != -1){
									if(v[3] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1 + 1][z1];
										newBase = new Cells(x1, y1 + 1, z1);
										System.out.println("新的基准点为:cells["+ x1 +"][" + (y1 + 1) + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1][y1 + 1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[4]最大
							if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
									&& v[4] >= v[3] && v[4] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[4] != -1){
									if(v[4] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1][z1 - 1];
										newBase = new Cells(x1, y1, z1 - 1);
										System.out.println("新的基准点为:cells["+ x1 +"][" + y1 + "][" + (z1 - 1) + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1][y1][z1 - 1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[5]最大
							if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
									&& v[5] >= v[3] && v[5] >= v[4]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[5] != -1){
									if(v[5] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1][z1 + 1];
										newBase = new Cells(x1, y1, z1 + 1);
										System.out.println("新的基准点为:cells["+ x1 +"][" + y1 + "][" + (z1 + 1) + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1][y1][z1 + 1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}	
							}
							
							//如果基准网格周围六个小网格都被搜索过了，单独作为一类
							if(v[0] == -1 && v[1] == -1 && v[2] == -1 && 
									v[3] == -1 && v[4] == -1 && v[4] == -1){
								center = new Cells(x1, y1, z1);
								System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
								newBase = null;
							}

						} while (counter <= path.size() && newBase != null && oldBase !=newBase);

						//聚类中心统计过程
						if (!centers.contains(center)){					
							cluster = new Cluster(center.getX(), center.getY(), center.getZ());
							cluster.setContentList(path);
							System.out.println("该聚类范围域内小网格的个数：" + cluster.getContentList().size());
							centers.add(center);
							initialClusters.add(cluster);
							System.out.println("中心的数量为：" + centers.size());
							System.out.println("初始聚类集合内聚类的数量为：" + initialClusters.size());							
						}//内层if结束
					}//外层if结束
					
				}
			}
		}// for循环结束
		
		System.out.println("有多少个初始聚类中心网格：" + centers.size());
		
		// 初始聚类集合中聚类的个数
		initalClusterNum = 1;

		// 存储聚类后的所有图像数据点
		DataPoint[] toWrite;
		
		//初始聚类集合中聚类的个数
		initalClusterNum = initialClusters.size();

		System.out.println("初始聚类集合中聚类的个数：" + initalClusterNum);
		
		//为每一个初始聚类设定一个编号，编号从1开始
		for(int x = 0, clusterNumber = 1; x < initalClusterNum; x++){
			initialClusters.get(x).setClusterLabel(clusterNumber++);
		}

		// 符合的条件的最终聚类集合
		LinkedList<Cluster> finnalInitList = new LinkedList<Cluster>();
		
		// 最终聚类集合中聚类的个数
		finalClusterNum = 1;
		
		for (int i = 0; i < initalClusterNum; i++) {
			// 如果初始聚类中心小网格的自适应性密度函数adaptive[x][y][z] < noiseLevel的话，则该聚类结果归结为噪声

			// 获取某个聚类的聚类中心小网格的x、y、z坐标
			int x = initialClusters.get(i).initialClusterCenter.getX();
			int y = initialClusters.get(i).initialClusterCenter.getY();
			int z = initialClusters.get(i).initialClusterCenter.getZ();
			double a = adaptive[x][y][z];
			System.out.println("adaptive[x][y][z] = " + a);
			
			if (adaptive[x][y][z] >= noiseLevel) {
				finnalInitList.add(initialClusters.get(i));
			}

		}// for循环结束

		// 最终聚类的个数
		finalClusterNum = finnalInitList.size();
		System.out.println("最终聚类集合中聚类的个数：" + finalClusterNum);
		
		// allPointsNum是原始图像数据点的个数
		toWrite = new DataPoint[allPointsNum];
		System.out.println("原始图像数据点个数：(toWrite.length)  :" + toWrite.length);		
	
		// 初始化
		for (int i = 0; i < allPointsNum; i++) {
			toWrite[i] = LUVPoints[i];
		}
		
		// 统计聚类算法处理的数据点的总量
		int toWriteS = 0;

		// 统计聚类整理过程中处理的数据点的总量
		int j = 0;

		for (int i = 0; i < finalClusterNum; i++) {
			//System.out.println("第"+ (i + 1) + "个聚类情况如下: ");
			Cluster cluster = finnalInitList.get(i);// 获取第(i+1)个聚类
						
			// 聚类范围域(上下文)网格坐标集,也就是这个聚类结果范围内所有的小网格对象的集合
			LinkedList<Cells> list = cluster.getContentList();
			
			System.out.println("第"+ (i + 1) + "个聚类范围域内包含 " + list.size() +"个小网格");
			
			// 为这个聚类范围域内所有小网格(坐标)设置类号,类号从1开始
			// coordNum为小网格的数量
			for (int coordNum = 0; coordNum < list.size(); coordNum++) {
				cells[list.get(coordNum).getX()][list.get(coordNum).getY()][list
						.get(coordNum).getZ()].setClusterNum(i + 1);//设置小网格的类号：clusterNum
			}

			int clusterSize = 0;// 聚类内数据点个数

			// 按聚类结果整理所有聚类算法处理到的数据点的相关属性(L,U,V)，类编号ClusterLabel
			for (Cells c : list) {    // 遍历聚类范围域内所有的小网格				
				Cells grid = cells[c.getX()][c.getY()][c.getZ()];					
				ArrayList<DataPoint> plist = grid.list;// 存储每个小网格内数据点的集合				
				clusterSize += plist.size();// 每个小网格内数据点个数累加和 = 聚类内数据点的个数
				
				// 更改点类编号
				for (DataPoint p : plist) {
					// 问题：??
					p.setClusterLabel(i);// 聚类内数据点的聚类标签号与它所在的聚类的编号相同，与初始聚类编号不同
					toWrite[(p.getLineNum() - 1)] = p;
					toWriteS++;
				}// 内层foreach循环结束
				
			}// 外层foreach循环结束

			j += clusterSize;
			System.out.println(" 第  " + (int) (i + 1) + " 类" + "共有 "
					+ clusterSize + " 个数据点");
			cluster.setPointsAll(clusterSize);
		}
		
		
		System.out.println("------共有 "
				+ j
				+ " 个点记录,"
				+ " 平滑前损失点 "
				+ (LUVPoints.length - j)
				+ " 个,"
				+ " 占总数 "
				+ new DecimalFormat("0.0000")
						.format((double) (LUVPoints.length - j)
								/ LUVPoints.length) + " ------");

		System.out.println(" -----共有  " + toWriteS + " 个点写出------");

		/**
		 * 第6步 将范围域聚类结果映射到图像空间 0、调用函数transformToImagePoints，
		 * 找回聚类过程中的所有损失点并将边界网格点按其与已知各类欧式距离最短的原则重新分配类编号；
		 * 1、调用函数smoothpix对分类结果进行平滑，去噪，修改所有数据点的对应类编号；
		 * 2、调用函数initSegments，依据数据点的类编号,初始化所有图片分块,确定每个块初始块号,上下文点集合和相邻块号的集合；
		 * 3、调用函数segmentsArrange整理统计分块结果,为图像数据点设定组号并更按分组改点的特征三维坐标；
		 * 4、调用函数LUVToRGB特征坐标回归为(R,G,B)；
		 * 5、初始化绘图API，调用函数show_DataPoint_on_Image(BufferedImage image,DataPoint
		 * point,Color color)绘图显示聚类结果
		 */


		// 去噪函数:1.损失点找回 2.边界网格点重新分类 3.将数据点初始化输出图像的像素点
		//toWrite = transformToimagePoints(LUVPoints, finnalInitList);

		// 当前聚类结果拥有的所有数据点(包括损失点)
		DataPoint[] imagePoints = new DataPoint[toWrite.length];
		
		// 数据点的类号
		int clusterLabel = -1;
		
		for(int i = 0; i < toWrite.length; i++){
			// 以下是未损失点
			// 如果数据点不是损失点，那么获得该数据点的类号
			// 注明：损失点的类号是-1
			if (toWrite[i].getClusterLabel() != -1) {
				clusterLabel = toWrite[i].getClusterLabel();
			} 
			//将数据点初始化输出图像的像素点
			imagePoints[i] = new DataPoint(toWrite[i].getLine(),
					toWrite[i].getRow(), toWrite[i].getCoord_X(),
					toWrite[i].getCoord_Y(), toWrite[i].getCoord_Z(), 1, 1,
					toWrite[i].getColor(), toWrite[i].getLineNum(), clusterLabel,
					toWrite[i].getSegmentLabel(), toWrite[i].getSmoothed(),
					toWrite[i].getCombined(), toWrite[i].getSegmentChecked());
		}
		System.out.println("imagePoints.length = " + imagePoints.length);
				
		// 对聚类后的数据点进行平滑
		toWrite = smoothpix(toWrite, finnalInitList);

		System.out.println("平滑过程完成！！");
		// 初始化所有图片分块
		segments = initSegments(toWrite, segments);
		
		DataPoint[] classifiedPoints = new DataPoint[toWrite.length];
		
		// 依据分组结果修改数据点的特征空间值，组内数据点的坐标值与该组的几何中心坐标值相同
		for(Segment s : segments){
			// 指定组的几何中心坐标值
			double centerL = s.getXloc();
			double centerU = s.getYloc();
			double centerV = s.getZloc();
			// s内的所有数据点的LUV特征空间值都相同，等于s的几何中心坐标值
			for (DataPoint point : s.getContentList()) {
				point.setCoord_X(centerL);
				point.setCoord_Y(centerU);
				point.setCoord_Z(centerV);
				classifiedPoints[point.getLineNum() - 1] = point;
			}
		}
		System.out.println("最终分块的个数是 = " + segments.size());
		// 整理分块信息,为图像数据点设定组号
		// 返回分组后的数据点集合
		//toWrite = segmentsArrange(segments, toWrite);

		// 特征坐标回归为(R,G,B)，返回GRB坐标值的数据点集合
		toWrite = LUVToRGB(toWrite);
		System.out.println("toWrite= " + toWrite.length);
		// 激活绘图API, 准备绘图
		BufferedImage img = ImageGenerator.drawImage(null, DataBase.sLen, DataBase.sHei);
		Graphics g = img.getGraphics();
		end = System.currentTimeMillis();
		String time = Long.toString(end - start) + "ms";
		for (DataPoint point : toWrite) {
			ImageGenerator.show_DataPoint_on_Image(img, point, new Color(
					(int) point.getCoord_X(), (int) point.getCoord_Y(),
					(int) point.getCoord_Z()));
		}

		g.setColor(Color.black);

		g.drawString("共" + toWrite.length + "个点，聚成" + finnalInitList.size()
				+ "个类", 0, DataBase.sHei - 45);

		g.drawString(time, 0, DataBase.sHei - 5);
		g.dispose();

		g.dispose();
		DataField.updateImagePanel(img);

	}// pocess函数结束

	
	/**
	 * 函数segmentsArrange 
	 * 将块中的数据点按其组编号归并 
	 * Parameters:
	 * allSegmentsGrouped:所有块分组后的结果集 
	 * allPoints:所有点集
	 * Return:分块结束,修改特征三维属性值后的所有数据点集合
	 */

	public DataPoint[] segmentsArrange(LinkedList<Segment> allSegmentsGrouped,
			DataPoint[] allPoints) {
		// 分组后的数据点集合
		DataPoint[] classifiedPoints = new DataPoint[allPoints.length];
		//所有组集合
		LinkedList<Group> allGroups = new LinkedList<Group>();

		// allSegmentsGrouped.size()是组内块的个数
		for (int i = 0; i < allSegmentsGrouped.size(); i++) {
			int included = -1;// 标志位，如果该块已经在某一分组中，则把该分组号码赋给它

			// 判断该块是否已在已建立的组中（通过组号判断）
			for (int j = 0; j < allGroups.size(); j++) {
				// 如果该块已经在某个分组中，则对应的组号相同
				if (allGroups.get(j).getGroupLabel() == allSegmentsGrouped.get(
						i).getGroupLabel()) {
					included = j;
				}
			}// for结束

			/*
			 * 下面的操作是把块加入到相应的分组中
			 */
			// 如果该块不在某个分组中，则建立新的分组
			if (included == -1) {
				// 疑问：它已经有组号了，怎么还没有相应的分组呢？？答：这是一个不断建立新分组，并且要避免分组重复的过程
				Group newGroup = new Group(allSegmentsGrouped.get(i)
						.getGroupLabel());
				// 为新的组建立数据点集合,for{...}在此处防止指针传递(浅复制)
				LinkedList<DataPoint> segentsContentMirrow = new LinkedList<DataPoint>();// 存放新的分组中所有的数据点
				for (int all = 0; all < allSegmentsGrouped.get(i)
						.getContentList().size(); all++) {
					segentsContentMirrow.add(allSegmentsGrouped.get(i)
							.getContentList().get(all));
				}
				newGroup.setContentList(segentsContentMirrow);
				newGroup.setNumAll();
				allGroups.add(newGroup);// 把新的分组加入到分组集合中
				// 否则加入现有分组
			} else {
				Group alreadyExistedGroup = allGroups.get(included);// 获取该块所在的组
				LinkedList<DataPoint> newContentList = alreadyExistedGroup
						.getContentList();// 获取组内所有数据点

				/*
				 * public boolean addAll(Collection<? extends E> c) 添加指定
				 * collection 中的所有元素到此列表的结尾，顺序是指定 collection 的迭代器返回这些元素的顺序。
				 * 如果指定的 collection 在操作过程中被修改，则此操作的行为是不确定的。
				 */
				newContentList.addAll(allSegmentsGrouped.get(i)
						.getContentList());
				alreadyExistedGroup.setContentList(newContentList);// 设定该块所在组的数据点集合
				alreadyExistedGroup.setNumAll();
				/*
				 * public E set(int index, E element) 将此列表中指定位置的元素替换为指定的元素。
				 * index - 要替换的元素的索引 element - 要在指定位置存储的元素 return 以前在指定位置的元素
				 */
				allGroups.set(included, alreadyExistedGroup);// 把指定组修改之后，放回原来的位置
			}
		}// 外层for循环结束

		// 依据分组结果修改数据点的特征空间值，组内数据点的坐标值与该组的几何中心坐标值相同
		for (Group group : allGroups) {
			// 指定组的几何中心坐标值
			double centerL = group.getXloc();
			double centerU = group.getYloc();
			double centerV = group.getZloc();

			// group内的所有数据点的LUV特征空间值都相同，等于group的几何中心坐标值
			for (DataPoint point : group.getContentList()) {
				point.setCoord_X(centerL);
				point.setCoord_Y(centerU);
				point.setCoord_Z(centerV);
				classifiedPoints[point.getLineNum() - 1] = point;
			}
		}

		// 最终分块的数目
		finalSegmentNum = allGroups.size();
		System.out.print("After arrange allSegmentsGrouped groupLabel: [");
		for (Group g : allGroups) {
			System.out.print(g.getGroupLabel() + " ");
		}
		System.out.println("]");
		System.out.println("After arrange allGroups size()= "
				+ allGroups.size());

		return classifiedPoints;
	}

	
	/**
	 * 函数outputGraphStructuretoFile 将无向图状结构输入文本,此处顶点为所有块，若两个块相邻则它们之间存在一条边
	 */
	public static void outputGraphStructuretoFile(LinkedList<Segment> segments) {
		String singleLine = new String();
		try {

			File f = new File("D:\\Graduation Project\\src\\0331\\"
					+ "ClusterForTD\\GraphStructure.txt");

			if (f.exists()) {
				System.out.println("GraphStructure存在，删除");
				f.delete();
			} else {
				System.out.println("GraphStructure不存在，正在创建...");
				if (f.createNewFile()) {
					System.out.println("GraphStructure创建成功！");
				} else {
					System.out.println("GraphStructure创建失败！");
				}
			}

			BufferedWriter output = new BufferedWriter(new FileWriter(f));
			output.write("*Vertices " + segments.size() + "\r\n");
			for (int num = 0; num < segments.size(); num++) {
				singleLine = "";
				singleLine = (num + 1) + "\r\n";
				output.write(singleLine);
			}
			output.write("*Edges" + "\r\n");
			for (int i = 0; i < segments.size(); i++) {
				singleLine = "";
				for (int j = 0; j < segments.get(i).getNeighborSegLabelList()
						.size(); j++) {
					Integer neighboor = segments.get(i)
							.getNeighborSegLabelList().get(j);
					if ((int) neighboor > i) {
						singleLine += (i + 1) + "   " + ((int) neighboor + 1);
						singleLine += "\r\n";
					}
				}
				output.write(singleLine);
			}
			output.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * 函数initSegments 
	 * 依据数据点的聚类结果初始化所有分块,包括各块的编号,数据点集合和邻块编号集合 
	 * Parameters:
	 * points:聚类平滑后的所有数据 
	 * segments:初始块列表 
	 * Return:初始分块结果集
	 */
	public LinkedList<Segment> initSegments(DataPoint[] points,
			LinkedList<Segment> segments) {
		// 块编号（自增长）
		int segmentIndex = 0;

		for (int i = 0; i < points.length; i++) {
			// 单块大小
			int sizeofOnePiece = 0;
			int k = 0;
			
			if (points[i].getSegmentChecked() == 0) {//检查是否已经被分块
				boolean tag;
				// 图像中单个块
				Segment oneSeg = new Segment(segmentIndex);

				// 加入该点自身
				oneSeg.getContentList().add(points[i]);// 块内数据点的集合
				points[i].setSegmentChecked(1);// 标记第i个数据点已经分配块了
				points[i].setSegmentLabel(segmentIndex);// 第i个数据点分配块号

				// 上下文及边界的初始镜像
				// 上下文数据点集合
				LinkedList<DataPoint> mirrorContentList = oneSeg.getContentList();
				System.out.println("初始上下文数据点集合内数据点个数：" + mirrorContentList.size());
				// 边界数据点集合
				LinkedList<DataPoint> mirrorMarginList = oneSeg.getMarginList();
				do {
					tag = true;
					// 用与八邻域判定，temp记录该数据点的全部相邻数据点,及完整的块
					//LinkedList<DataPoint> temp = allFourNeighborFeilds(mirrorContentList.get(k));//返回第k个数据点
					LinkedList<DataPoint> temp = allEightNeighborFeilds(mirrorContentList.get(k));//返回第k个数据点
					System.out.println("temp.size() =" + temp.size());
					
					for (DataPoint p : temp) {
						//System.out.println("p.getClusterLabel() =" + p.getClusterLabel());
						//System.out.println("points[i].getClusterLabel() =" + points[i].getClusterLabel());
						// 把与被考察点属于同一聚类(即同一块)的数据点加入邻域
						if (p.getClusterLabel() == points[i].getClusterLabel() ) {
							// 只有那些没有被分块的数据点
							if (p.getSegmentChecked() == 0) {
								mirrorContentList.add(p);
								p.setSegmentChecked(1);
								p.setSegmentLabel(segmentIndex);
							}
						}
						// 若某个点的四邻域(或八邻域)中有点属于其他类或者块,则该点为所属块的边界点
						else {
							if (!mirrorMarginList.contains(mirrorContentList.get(k)))
								// 把那些标准点八邻域范围内与标准点不同聚类的点，加入到边界点集中
								mirrorMarginList.add(mirrorContentList.get(k));
						}
					}
					k++;
					sizeofOnePiece = mirrorContentList.size();
					
					if (k == sizeofOnePiece) {
						// 块自身已经搜索完毕不再扩大，跳出循环
						tag = false;
					}
				} while (tag == true);
				
				// 引用传递,已修改
				oneSeg.setContentList(mirrorContentList);
				//System.out.println("分块之后上下文数据点集合内数据点个数：" + mirrorContentList.size());
				oneSeg.setMarginList(mirrorMarginList);
				//System.out.println("分块之后边界数据点集合内数据点个数：" + mirrorMarginList.size());
				oneSeg.setNumAll();
				segments.add(oneSeg);
				//System.out.println("分块个数：" + segments.size());
				segmentIndex++;//块编号不断增加
			}//if结束

		}//for循环结束

		System.out.println("分块个数：" + segments.size());
		
		// 第二次扫描，初始化每个块的相邻块信息
		for (int j = 0; j < segments.size(); j++) {
			// 边界及邻域的初始镜像
			// 相邻块的编号集合
			LinkedList<Integer> mirrorNeighborSegLabelList = segments.get(j)
					.getNeighborSegLabelList();
			
			// 被考察数据块的邻域数据点集合
			LinkedList<DataPoint> mirrorMarginList = segments.get(j)
					.getMarginList();

			//遍历邻域所有数据点
			for (DataPoint marginPoint : mirrorMarginList) {
				//统计邻域数据点四邻域数据点集合temp2
				//LinkedList<DataPoint> temp2 = allFourNeighborFeilds(marginPoint);
				LinkedList<DataPoint> temp2 = allEightNeighborFeilds(marginPoint);
				//遍历四邻域数据点集合
				for (DataPoint neighboorPoint : temp2) {
					// 把与被考察点不同块的数据块号记录下来
					if (neighboorPoint.getSegmentLabel() != marginPoint
							.getSegmentLabel()) {
						// 获取与被考察点块号不同的那些块号
						Integer neighboorSegLab = new Integer(
								neighboorPoint.getSegmentLabel());
						if (!mirrorNeighborSegLabelList
								.contains(neighboorSegLab)
								|| mirrorNeighborSegLabelList.size() == 0) {
							// 与被考察点不同块的数据块号的集合
							mirrorNeighborSegLabelList.add(neighboorSegLab);							
						}
					}
				}
			}
			segments.get(j).setNeighborSegLabelList(mirrorNeighborSegLabelList);
			System.out.println("与第" + j + "块相邻块的编号有："+ 
					mirrorNeighborSegLabelList.size() + "个");
		}
		return segments;
	}

	
	/**
	 * 函数transformToImagePoints 
	 * 将数据点初始化输出图像的像素点 
	 * Parameters:
	 * points:当前图像LUV坐标系所有的数据点 
	 * clusters：聚类集合 
	 * cells：小网格对象 
	 * return 输出图像的像素点集合
	 */
	public static DataPoint[] transformToimagePoints(DataPoint[] points,
			LinkedList<Cluster> clusters) {

		// 当前聚类结果拥有的所有数据点(包括损失点)
		DataPoint[] imagePoints = new DataPoint[points.length];

		// 数据点的类号
		int clusterLabel;

		for (int i = 0; i < points.length; i++) {
			
			// 以下是未损失点
			// 如果数据点不是损失点，那么获得该数据点的类号(损失点没有类号：clusterLabel )
			// 注明：损失点的类号是-1
			if (points[i].getClusterLabel() != -1) {
				clusterLabel = points[i].getClusterLabel();
			} 
			
			// 以下是损失点
			// 如果数据点本身就是损失点，那么计算该数据点到每个聚类几何中心坐标欧式几何距离，找回损失点
			else { 

				// 每个聚类中心坐标值(clusterCenterX, clusterCenterY, clusterCenterZ)
				double clusterCenterX = 0;
				double clusterCenterY = 0;
				double clusterCenterZ = 0;

				// 损失点与聚类中心间的最短欧式距离
				double minDistance = Double.MAX_VALUE;
				//System.out.println("minDistance = " +minDistance);
				
				// 归并类号,及满足损失点与聚类中心间的欧式距离最短的聚类的类号
				int minDistanceCluster = 0;

				//遍历所有的聚类
				for (int j = 0; j < clusters.size(); j++) {
					double dis = 0;
					
					// 该聚类几何中心小网格（不是山峰封顶的小网格，也就是说不是聚类中心）
					Cluster cluster = clusters.get(j);
					System.out.println("第" + (j + 1) + "个聚类范围域内有： " +cluster.contentList.size()+ "个小网格，" +
							"有： " +cluster.getPointsAll()+ " 个数据点" );
					//确定了该聚类的几何中心小网格
					cluster.adjustClusterCenter();
					//获取几何中心小网格
					Cells geoCenter = cluster.getClusterCenter();

					// 获取聚类几何中心坐标值，把这个坐标值作为聚类的坐标值
					/*
					clusterCenterX = cells[geoCenter.getX()][geoCenter.getY()][geoCenter
							.getZ()].getXloc();
					clusterCenterY = cells[geoCenter.getX()][geoCenter.getY()][geoCenter
							.getZ()].getYloc();
					clusterCenterZ = cells[geoCenter.getX()][geoCenter.getY()][geoCenter
							.getZ()].getZloc();
					*/
					clusterCenterX = geoCenter.getXloc();
					clusterCenterY = geoCenter.getYloc();
					clusterCenterZ = geoCenter.getZloc();

					// 计算距离(欧式距离)
					dis = StrictMath.pow(
							(clusterCenterX - points[i].getCoord_X()), 2.0)
							+ StrictMath.pow(
									(clusterCenterY - points[i].getCoord_Y()),
									2.0)
							+ StrictMath.pow(
									(clusterCenterZ - points[i].getCoord_Z()),
									2.0);

					if (j == 0 || dis < minDistance) {
						minDistance = dis;
						minDistanceCluster = j;// j是聚类类号
					}

				}// for循环结束
				//以上完成了寻找距离最短的聚类的类号
				clusterLabel = minDistanceCluster;
			}// else结束
			//以上完成了为图像中没有一个数据点分配一个类号
			System.out.println("第" + i + "个数据点的类号是：" +clusterLabel + "， 为第"+ i +"个数据点分配类号完成！！");
			
			//将数据点初始化输出图像的像素点
			imagePoints[i] = new DataPoint(points[i].getLine(),
					points[i].getRow(), points[i].getCoord_X(),
					points[i].getCoord_Y(), points[i].getCoord_Z(), 1, 1,
					points[i].getColor(), points[i].getLineNum(), clusterLabel,
					points[i].getSegmentLabel(), points[i].getSmoothed(),
					points[i].getCombined(), points[i].getSegmentChecked());
		}// for循环结束，完成了对每个数据点的操作

		return imagePoints;
	}

	/**
	 * 函数smoothpix 
	 * 完成数据聚类的平滑去噪 
	 * Parameters: 
	 * points: 聚类搜索后待平滑的数据点集 
	 * clusters: 搜索后的聚类列表 
	 * Return 平滑后的数据集合
	 */
	public DataPoint[] smoothpix(DataPoint[] points,
			LinkedList<Cluster> clusters) {
		smoothValve = 40;
		// 初始化allPix,原始图像的坐标值
		for (int i = 0; i < points.length; i++) {
			allPix[(int) points[i].getLine()][(int) points[i].getRow()] = points[i];
		}

		for (int i = 0; i < points.length; i++) {

			// 记录被考察点邻域外围(边界)的所有数据点的聚类信息
			int[] neighboorCluster = new int[finalClusterNum + 1];

			// 被考察点的邻域所有数据点的集合
			LinkedList<DataPoint> neighboor = new LinkedList<DataPoint>();

			// 被考察点的邻域的大小(被考察点邻域聚类内数据点的个数)
			int NBsize = 0;

			// 包含被考察点邻域外围(边界)数据点数量最多的那一类的编号，需要平滑的数据点将被平滑至该类(基准类)
			int NBMax = 0;

			int k = 0;// 邻域计数

			// 初始化被考察点邻域外围的所有数据点的聚类信息（初始化数组对象）,neighboorCluster.length表示邻域个数
			for (int j = 0; j < neighboorCluster.length; j++) {
				neighboorCluster[j] = 0;// 第j个聚类对象初始化为0
			}

			if (points[i].getSmoothed() == 0) {// 如果该数据点没有被平滑
				boolean tag;
				// 加入该点自身
				neighboor.add(points[i]);
				do {
					tag = true;
					// 先找到邻域的编号为k的数据点，调用allFourNeighborFeilds函数，返回指定数据点四邻域数据点集合
					//LinkedList<DataPoint> temp = allFourNeighborFeilds(neighboor.get(k));// 用与 四邻域判定，temp记录该数据点的全部相邻数据点
					LinkedList<DataPoint> temp = allEightNeighborFeilds(neighboor.get(k));// 用与 四邻域判定，temp记录该数据点的全部相邻数据点

					//依次遍历四邻域内每一个数据点 
					for (DataPoint p : temp) {
						// 把与被考察点属于同一聚类的数据点加入邻域(根据聚类标签判断)
						if (p.getClusterLabel() == points[i].getClusterLabel()) {
							if (!neighboor.contains(p)) {// 如果指定数据点的邻域不包含四邻域数据点集合中的点，那么将这个点加入到指定点的邻域
								neighboor.add(p);
							}
						}
						// 如果被考察数据点与它的八邻域数据点不属于同一个聚类的话，将其计入邻域边界，并定位基准类
						else {
							int clusterLabel = p.getClusterLabel();
							if (clusterLabel != -1) {
								clusterLabel++;
							} else {
								clusterLabel = 0;
							}
							neighboorCluster[clusterLabel]++;
							// 包含被考察点邻域外围(边界)数据点数量最多的那一类的编号，需要平滑的数据点将被平滑至该类(基准类)
							if (NBMax == -1) {
								NBMax++;
							}
							if (neighboorCluster[clusterLabel] > neighboorCluster[NBMax]) {
								clusterLabel--;
								NBMax = clusterLabel;
							}
						}// else结束
					}// foreach循环结束
					k++;
					// 被考察点邻域聚类内数据点的个数
					NBsize = neighboor.size();
					System.out.println("被考察点邻域聚类内数据点的个数： " + NBsize);
					if (k == NBsize) {
						// 邻域已经搜索完毕不再扩大，跳出循环
						tag = false;
					}
					// 邻域大小超过平滑界限40或邻域停止增长，循环终止
				} while (NBsize < smoothValve && tag == true);

				/*
				 * 除噪方法:向最大临近域归并 将需要平滑和除噪的数据点及其邻域所有点平滑至标准类 需要平滑和除噪的点 
				 * 1、NBsize < 40 表示邻域低于平滑界限的噪声点； 
				 * 2、neighboor.get(0).getClusterLabel() == -1表示该邻域中均为聚类算法损失的图像数据点；
				 */
				if (NBsize < smoothValve
						|| neighboor.get(0).getClusterLabel() == -1) {
					// 如果基准类的类号不等于-1
					if (NBMax != -1) {
						for (DataPoint p : neighboor) {
							// 将被考察点的邻域所有数据点的类号设置为被考察点邻域外围(边界)数据点数量最多的那一类的编号
							p.setClusterLabel(NBMax);
						}
					}
					// 聚类过程中损失的噪声点设为缺省值
					else {
						for (DataPoint p : neighboor) {
							p.setClusterLabel(NBMax);
						}
					}
				}// 外层if结束

				// 对于不需要平滑的点，由于邻域的互通性，不需要再次参与下次判定
				else {
					for (DataPoint p : neighboor) {
						p.setSmoothed(1);// 不需要平滑的数据点，设置为无需平滑
					}
				}
			}
		}// 最外层for循环结束

		// 原图像数据矩阵更新
		for (int i = 0; i < points.length; i++) {
			allPix[(int) points[i].getLine()][(int) points[i].getRow()] = points[i];
		}

		return points;
	}

	/**
	 * 函数allEightNeighborFeilds 搜索指定数据点在原图像中的八邻域 Parameters: point:指定的数据点
	 * Return:八邻域点集
	 */
	public LinkedList<DataPoint> allEightNeighborFeilds(DataPoint point) {
		LinkedList<DataPoint> allEightNeighboor = new LinkedList<DataPoint>();
		// 加入自身
		allEightNeighboor.add(point);
		// 添加上方行的点
		if (point.getLine() > 1) {
			allEightNeighboor.add(allPix[(int) point.getLine() - 1][(int) point
					.getRow()]);
			if (point.getRow() > 1)
				allEightNeighboor
						.add(allPix[(int) point.getLine() - 1][(int) point
								.getRow() - 1]);
			if (point.getRow() < row_num)
				allEightNeighboor
						.add(allPix[(int) point.getLine() - 1][(int) point
								.getRow() + 1]);
		}
		// 添加本行的点
		if (point.getRow() > 1)
			allEightNeighboor.add(allPix[(int) point.getLine()][(int) point
					.getRow() - 1]);

		if (point.getRow() < row_num)
			allEightNeighboor.add(allPix[(int) point.getLine()][(int) point
					.getRow() + 1]);
		// 添加下方行的点
		if (point.getLine() < line_num) {
			allEightNeighboor.add(allPix[(int) point.getLine() + 1][(int) point
					.getRow()]);
			if (point.getRow() > 1)
				allEightNeighboor
						.add(allPix[(int) point.getLine() + 1][(int) point
								.getRow() - 1]);
			if (point.getRow() < row_num)
				allEightNeighboor
						.add(allPix[(int) point.getLine() + 1][(int) point
								.getRow() + 1]);
		}
		return allEightNeighboor;
	}

	/**
	 * 函数allFourNeighborFeilds 
	 * 搜索指定数据点在原图像中的四邻域 
	 * Parameters: point:指定的数据点
	 * Return:四邻域点集
	 */
	public LinkedList<DataPoint> allFourNeighborFeilds(DataPoint point) {
		LinkedList<DataPoint> allFourNeighboor = new LinkedList<DataPoint>();
		// 加入自身
		allFourNeighboor.add(point);
		// 添加上方的点
		if (point.getLine() > 1)
			allFourNeighboor.add(allPix[(int) point.getLine() - 1][(int) point
					.getRow()]);
		// 添加左边的点
		if (point.getRow() > 1)
			allFourNeighboor.add(allPix[(int) point.getLine()][(int) point
					.getRow() - 1]);
		// 添加右边的点
		if (point.getRow() < row_num)
			allFourNeighboor.add(allPix[(int) point.getLine()][(int) point
					.getRow() + 1]);
		// 添加下方的点
		if (point.getLine() < line_num)
			allFourNeighboor.add(allPix[(int) point.getLine() + 1][(int) point
					.getRow()]);

		return allFourNeighboor;
	}

	

	/**
	 * 函数initialClusterCenter 
	 * 搜索初始聚类中心 
	 * Parameters: 
	 * c 指定的某个小网格对象 
	 * return 初始聚类中心小网格center对象
	 * author weil
	 */
	public Cells initialClusterCenter(Cells c) {
		// 计算各小网格的x偏导，y偏导，z偏导
		partialDerivative();

		// 聚类中心
		Cells center = new Cells();

		// 存储路径上的所有的小网格对象
		LinkedList<Cells> path = new LinkedList<Cells>();

		// 加入自身
		path.add(c);

		// 旧的基准点
		Cells oldBase = new Cells();

		// 新的基准点
		Cells newBase = new Cells();

		// 计数器
		int counter = 0;

		do {
			oldBase = path.get(counter);

			// 基准小网格对象的坐标值
			int x1 = oldBase.getX();
			int y1 = oldBase.getY();
			int z1 = oldBase.getZ();

			double[] v = { -1, -1, -1, -1, -1, -1 };

			v[0] = gapLeftCell(x1, y1, z1);
			v[1] = gapRightCell(x1, y1, z1);
			v[2] = gapUpCell(x1, y1, z1);
			v[3] = gapDownCell(x1, y1, z1);
			v[4] = gapFrontCell(x1, y1, z1);
			v[5] = gapBackCell(x1, y1, z1);
			
			// v[0]最大
			if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
					&& v[0] >= v[4] && v[0] >= v[5]) {
				/*
				 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
				 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
				 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
				 * == null)时不成立
				 */
				if (v[0] < 0.5) {
					//center = oldBase;
					center = new Cells(x1, y1, z1);
					newBase = null;
				} else {
					//newBase = cells[x1 - 1][y1][z1];
					newBase = new Cells(x1 - 1, y1, z1);
					if (!path.contains(newBase)) {
						path.add(newBase);
						// 爬山过程中已经搜索过的小网格标记为“已搜索”
						cells[x1 -1][y1][z1].setAlreadysearched(1);
						counter++;
					}
				}

			}// if结束

			// v[1]最大
			if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
					&& v[1] >= v[4] && v[1] >= v[5]) {
				/*
				 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
				 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
				 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
				 * == null)时不成立
				 */
				if (v[1] < 0.5) {
					//center = oldBase;
					center = new Cells(x1, y1, z1);
					newBase = null;
				} else {
					//如何把一个对象赋给另一个对象
					//newBase = cells[x1 + 1][y1][z1];
					newBase = new Cells(x1 + 1, y1, z1);
					if (!path.contains(newBase)) {
						path.add(newBase);
						// 爬山过程中已经搜索过的小网格标记为“已搜索”
						cells[x1 + 1][y1][z1].setAlreadysearched(1);
						counter++;
					}
				}

			}

			// v[2]最大
			if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
					&& v[2] >= v[4] && v[2] >= v[5]) {
				/*
				 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
				 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
				 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
				 * == null)时不成立
				 */
				if (v[2] < 0.5) {
					//center = oldBase;
					center = new Cells(x1, y1, z1);
					newBase = null;
				} else {
					//newBase = cells[x1][y1 - 1][z1];
					newBase = new Cells(x1, y1 - 1, z1);
					if (!path.contains(newBase)) {
						path.add(newBase);
						// 爬山过程中已经搜索过的小网格标记为“已搜索”
						cells[x1][y1 - 1][z1].setAlreadysearched(1);
						counter++;
					}
				}

			}

			// v[3]最大
			if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
					&& v[3] >= v[4] && v[3] >= v[5]) {
				/*
				 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
				 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
				 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
				 * == null)时不成立
				 */
				if (v[3] < 0.5) {
					//center = oldBase;
					center = new Cells(x1, y1, z1);
					newBase = null;
				} else {
					//newBase = cells[x1][y1 + 1][z1];
					newBase = new Cells(x1, y1 + 1, z1);
					if (!path.contains(newBase)) {
						path.add(newBase);
						// 爬山过程中已经搜索过的小网格标记为“已搜索”
						cells[x1][y1 + 1][z1].setAlreadysearched(1);
						counter++;
					}
				}

			}

			// v[4]最大
			if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
					&& v[4] >= v[3] && v[4] >= v[5]) {
				/*
				 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
				 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
				 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
				 * == null)时不成立
				 */
				if (v[4] < 0.5) {
					//center = oldBase;
					center = new Cells(x1, y1, z1);
					newBase = null;
				} else {
					//newBase = cells[x1][y1][z1 - 1];
					newBase = new Cells(x1, y1, z1 - 1);
					if (!path.contains(newBase)) {
						path.add(newBase);
						// 爬山过程中已经搜索过的小网格标记为“已搜索”
						cells[x1][y1][z1 - 1].setAlreadysearched(1);
						counter++;
					}
				}

			}

			// v[5]最大
			if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
					&& v[5] >= v[3] && v[5] >= v[4]) {
				/*
				 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
				 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
				 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
				 * == null)时不成立
				 */
				if (v[5] < 0.5) {
					//center = oldBase;
					center = new Cells(x1, y1, z1);
					newBase = null;
				} else {
					//newBase = cells[x1][y1][z1 + 1];
					newBase = new Cells(x1, y1, z1 + 1);
					if (!path.contains(newBase)) {
						path.add(newBase);
						// 爬山过程中已经搜索过的小网格标记为“已搜索”
						cells[x1][y1][z1 + 1].setAlreadysearched(1);
						counter++;
					}
				}

			}

		} while (counter <= path.size() && newBase != null);
		
		return center;
	}
	
	

	/**
	 * 函数initialCluster 
	 * 搜索初始聚类中心 
	 * Parameters: null 
	 * return 所有初始聚类中心小网格center对象集合
	 * author weil
	 */
	public ArrayList<Cells> initialCluster() {
		// 计算各小网格的x偏导，y偏导，z偏导
		partialDerivative();

		System.out.println("计算偏导结束！");
		
		// 存储初始聚类中心小网格集合
		ArrayList<Cells> centers = new ArrayList<Cells>();
		
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					// 聚类中心
					Cells center = new Cells();

					// 存储路径上的所有的小网格对象
					LinkedList<Cells> path = new LinkedList<Cells>();
					
					System.out.println("打印路径长度(添加自身之前): path.size() = " + path.size());
					
					// 加入自身
					path.add(new Cells(x, y, z));

					System.out.println("打印路径长度（添加自身之后）: path.size() = " + path.size());
					System.out.println("cells[" + path.get(0).getX() + "][" + path.get(0).getY() + "][" + path.get(0).getZ() +"]");
					
					// 旧的基准点初始化
					Cells oldBase = new Cells();
					oldBase = null;
					
					// 新的基准点初始化
					Cells newBase =new Cells();
					newBase = null;
					
					// 计数器
					int counter = 0;

					//alreadysearched == 0 表示没有被搜索过。
					if(cells[x][y][z].getAlreadysearched() == 0){
						double[] v = { -1, -1, -1, -1, -1, -1 };
						
						//LinkedList<Cells> contentList = new LinkedList<Cells>();
						
						do {
							
							oldBase = path.get(counter);

							// 基准小网格对象的坐标值
							int x1 = oldBase.getX();
							int y1 = oldBase.getY();
							int z1 = oldBase.getZ();

							System.out.println("基准点为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							cells[x1][y1][z1].setAlreadysearched(1);
							
							v[0] = gapLeftCell(x1, y1, z1);
							v[1] = gapRightCell(x1, y1, z1);
							v[2] = gapUpCell(x1, y1, z1);
							v[3] = gapDownCell(x1, y1, z1);
							v[4] = gapFrontCell(x1, y1, z1);
							v[5] = gapBackCell(x1, y1, z1);
							//依次打印出来数组值
							for(int i = 0; i < v.length; i++){
								System.out.println("v[" + i +"] = "+ v[i]);
							}
							
							// v[0]最大
							if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
									&& v[0] >= v[4] && v[0] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[0] != -1){
									if(v[0] < 1.0){
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
									} else {
										//newBase = cells[x1 - 1][y1][z1];
										newBase = new Cells(x1 - 1, y1, z1);
										System.out.println("新的基准点为:cells["+ (x1 - 1) +"][" + y1 + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1 -1][y1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}// if结束

							
							// v[1]最大
							if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
									&& v[1] >= v[4] && v[1] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[1] != -1){
									if(v[1] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//如何把一个对象赋给另一个对象
										//newBase = cells[x1 + 1][y1][z1];
										newBase = new Cells(x1 + 1, y1, z1);
										System.out.println("新的基准点为:cells["+ (x1 + 1) +"][" + y1 + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1 + 1][y1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[2]最大
							if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
									&& v[2] >= v[4] && v[2] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[2] != -1){
									if(v[2] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1 - 1][z1];
										newBase = new Cells(x1, y1 - 1, z1);
										System.out.println("新的基准点为:cells["+ x1 +"][" + (y1 - 1) + "][" + z1 + "]");
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1][y1 - 1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}							
							

							// v[3]最大
							if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
									&& v[3] >= v[4] && v[3] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[3] != -1){
									if(v[3] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1 + 1][z1];
										newBase = new Cells(x1, y1 + 1, z1);
										System.out.println("新的基准点为:cells["+ x1 +"][" + (y1 + 1) + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1][y1 + 1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[4]最大
							if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
									&& v[4] >= v[3] && v[4] >= v[5]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[4] != -1){
									if(v[4] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1][z1 - 1];
										newBase = new Cells(x1, y1, z1 - 1);
										System.out.println("新的基准点为:cells["+ x1 +"][" + y1 + "][" + (z1 - 1) + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1][y1][z1 - 1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[5]最大
							if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
									&& v[5] >= v[3] && v[5] >= v[4]) {
								/*
								 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
								 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
								 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
								 * == null)时不成立
								 */
								if(v[5] != -1){
									if(v[5] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1][z1 + 1];
										newBase = new Cells(x1, y1, z1 + 1);
										System.out.println("新的基准点为:cells["+ x1 +"][" + y1 + "][" + (z1 + 1) + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// 爬山过程中已经搜索过的小网格标记为“已搜索”
											cells[x1][y1][z1 + 1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
										}
									}
								}	
							}

						} while (counter <= path.size() && newBase != null && oldBase !=newBase);

						//聚类中心统计过程
						if (!centers.contains(center)) {					

							centers.add(center);
							System.out.println("中心的数量为：" + centers.size());
						}//内层if结束
					}//外层if结束
					
				}
			}
		}// for循环结束

		return centers;
	}

	
	/**
	 * 函数searchPathCells 
	 * 每一条搜索路径上所有小网格对象 
	 * Parameters: c 指定路径上的某个小网格对象
	 * return 指定小网格对象所在路径上的所有的小网格对象 
	 * author weil
	 */
	public LinkedList<Cells> searchPathCells(Cells c) {

		// 存储路径上的所有的小网格对象
		LinkedList<Cells> path = new LinkedList<Cells>();

		// 加入自身
		path.add(c);

		// 路径上的聚类中心小网格对象
		Cells center = new Cells();
		
		// 旧的基准点初始化
		Cells oldBase = new Cells();

		// 新的基准点初始化
		Cells newBase =new Cells();

		// 计数器
		int counter = 0;

		//alreadysearched == 0 表示没有被搜索过。
		if(cells[c.getX()][c.getY()][c.getZ()].getAlreadysearched() == 0){
			double[] v = { -1, -1, -1, -1, -1, -1 };
			
			do {
				oldBase = path.get(counter);

				// 基准小网格对象的坐标值
				int x1 = oldBase.getX();
				int y1 = oldBase.getY();
				int z1 = oldBase.getZ();

				System.out.println("基准点为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
				cells[x1][y1][z1].setAlreadysearched(1);
				
				v[0] = gapLeftCell(x1, y1, z1);
				v[1] = gapRightCell(x1, y1, z1);
				v[2] = gapUpCell(x1, y1, z1);
				v[3] = gapDownCell(x1, y1, z1);
				v[4] = gapFrontCell(x1, y1, z1);
				v[5] = gapBackCell(x1, y1, z1);
				//依次打印出来数组值
				for(int i = 0; i < v.length; i++){
					System.out.println("v[" + i +"] = "+ v[i]);
				}
				
				// v[0]最大
				if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
						&& v[0] >= v[4] && v[0] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[0] != -1){
						if(v[0] < 1.0){
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1 - 1][y1][z1];
							newBase = new Cells(x1 - 1, y1, z1);
							System.out.println("新的基准点为:cells["+ (x1 - 1) +"][" + y1 + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1 -1][y1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}// if结束

				
				// v[1]最大
				if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
						&& v[1] >= v[4] && v[1] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[1] != -1){
						if(v[1] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//如何把一个对象赋给另一个对象
							//newBase = cells[x1 + 1][y1][z1];
							newBase = new Cells(x1 + 1, y1, z1);
							System.out.println("新的基准点为:cells["+ (x1 + 1) +"][" + y1 + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1 + 1][y1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[2]最大
				if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
						&& v[2] >= v[4] && v[2] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[2] != -1){
						if(v[2] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1 - 1][z1];
							newBase = new Cells(x1, y1 - 1, z1);
							System.out.println("新的基准点为:cells["+ x1 +"][" + (y1 - 1) + "][" + z1 + "]");
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1][y1 - 1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}							
				

				// v[3]最大
				if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
						&& v[3] >= v[4] && v[3] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[3] != -1){
						if(v[3] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1 + 1][z1];
							newBase = new Cells(x1, y1 + 1, z1);
							System.out.println("新的基准点为:cells["+ x1 +"][" + (y1 + 1) + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1][y1 + 1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[4]最大
				if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
						&& v[4] >= v[3] && v[4] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[4] != -1){
						if(v[4] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1][z1 - 1];
							newBase = new Cells(x1, y1, z1 - 1);
							System.out.println("新的基准点为:cells["+ x1 +"][" + y1 + "][" + (z1 - 1) + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1][y1][z1 - 1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[5]最大
				if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
						&& v[5] >= v[3] && v[5] >= v[4]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[5] != -1){
						if(v[5] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1][z1 + 1];
							newBase = new Cells(x1, y1, z1 + 1);
							System.out.println("新的基准点为:cells["+ x1 +"][" + y1 + "][" + (z1 + 1) + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1][y1][z1 + 1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}	
				}

			} while (counter <= path.size() && newBase != null && oldBase !=newBase);
		}

		return path;
	}
	
	
	/**
	 * 函数searchforCluster
	 * 搜索指定聚类范围域内的所有小网格
	 * Parameters: 
	 * cluster：指定聚类中心
	 * return : 聚类范围域内所有小网格
	 * 
	 */
	
	public LinkedList<Cells> searchforCluster(Cluster c){

		//获取该聚类的聚类中心
		Cells center = c.initialClusterCenter;
		System.out.println("x = " + center.getX() +", y = " + center.getY() + ", z = " + center.getZ());
		
		// 存储路径上的所有的小网格对象
		LinkedList<Cells> path = new LinkedList<Cells>();
		
		System.out.println("打印路径长度(添加自身之前): path.size() = " + path.size());
		
		// 加入自身
		path.add(new Cells(center.getX(), center.getY(), center.getZ()));

		System.out.println("打印路径长度（添加自身之后）: path.size() = " + path.size());
		System.out.println("cells[" + path.get(0).getX() + "][" + path.get(0).getY() + "][" + path.get(0).getZ() +"]");
		
		// 旧的基准点初始化
		Cells oldBase = new Cells();
		oldBase = null;
		
		// 新的基准点初始化
		Cells newBase =new Cells();
		newBase = null;
		
		// 计数器
		int counter = 0;

		//alreadysearched == 0 表示没有被搜索过。
		if(cells[center.getX()][center.getY()][center.getZ()].getAlreadysearched() == 0){
			double[] v = { -1, -1, -1, -1, -1, -1 };
				
			do {
				oldBase = path.get(counter);

				// 基准小网格对象的坐标值
				int x1 = oldBase.getX();
				int y1 = oldBase.getY();
				int z1 = oldBase.getZ();

				System.out.println("基准点为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
				cells[x1][y1][z1].setAlreadysearched(1);
				
				v[0] = gapLeftCell(x1, y1, z1);
				v[1] = gapRightCell(x1, y1, z1);
				v[2] = gapUpCell(x1, y1, z1);
				v[3] = gapDownCell(x1, y1, z1);
				v[4] = gapFrontCell(x1, y1, z1);
				v[5] = gapBackCell(x1, y1, z1);
				//依次打印出来数组值
				for(int i = 0; i < v.length; i++){
					System.out.println("v[" + i +"] = "+ v[i]);
				}
				
				// v[0]最大
				if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
						&& v[0] >= v[4] && v[0] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[0] != -1){
						if(v[0] < 1.0){
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
						} else {
							//newBase = cells[x1 - 1][y1][z1];
							newBase = new Cells(x1 - 1, y1, z1);
							System.out.println("新的基准点为:cells["+ (x1 - 1) +"][" + y1 + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1 -1][y1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}// if结束

				
				// v[1]最大
				if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
						&& v[1] >= v[4] && v[1] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[1] != -1){
						if(v[1] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//如何把一个对象赋给另一个对象
							//newBase = cells[x1 + 1][y1][z1];
							newBase = new Cells(x1 + 1, y1, z1);
							System.out.println("新的基准点为:cells["+ (x1 + 1) +"][" + y1 + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1 + 1][y1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[2]最大
				if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
						&& v[2] >= v[4] && v[2] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[2] != -1){
						if(v[2] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1 - 1][z1];
							newBase = new Cells(x1, y1 - 1, z1);
							System.out.println("新的基准点为:cells["+ x1 +"][" + (y1 - 1) + "][" + z1 + "]");
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1][y1 - 1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}							
				

				// v[3]最大
				if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
						&& v[3] >= v[4] && v[3] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[3] != -1){
						if(v[3] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1 + 1][z1];
							newBase = new Cells(x1, y1 + 1, z1);
							System.out.println("新的基准点为:cells["+ x1 +"][" + (y1 + 1) + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1][y1 + 1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[4]最大
				if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
						&& v[4] >= v[3] && v[4] >= v[5]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[4] != -1){
						if(v[4] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1][z1 - 1];
							newBase = new Cells(x1, y1, z1 - 1);
							System.out.println("新的基准点为:cells["+ x1 +"][" + y1 + "][" + (z1 - 1) + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1][y1][z1 - 1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[5]最大
				if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
						&& v[5] >= v[3] && v[5] >= v[4]) {
					/*
					 * 这里注明：当寻找到两个相邻小网格的梯度差值 < 0.01
					 * 时，认为该oldBase是路径的顶峰（即聚类中心小网格）
					 * 那么爬山过程不要继续下去了，也就是说newBase = null;这样在检查if(oldBase
					 * == null)时不成立
					 */
					if(v[5] != -1){
						if(v[5] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("该条爬山路径的顶峰为:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1][z1 + 1];
							newBase = new Cells(x1, y1, z1 + 1);
							System.out.println("新的基准点为:cells["+ x1 +"][" + y1 + "][" + (z1 + 1) + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// 爬山过程中已经搜索过的小网格标记为“已搜索”
								cells[x1][y1][z1 + 1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("该条路径上的小网格个数:path.size() =  " + path.size());
							}
						}
					}	
				}

			} while (counter <= path.size() && newBase != null && oldBase !=newBase);
		}
					
		return path;		
	}
	
	// 下面的六个函数是计算六邻域梯度差值的
	/**
	 * 函数gapLeftCell 
	 * 计算当前小网格与其左边小网格的差值 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return 当前小网格与其左边小网格的差值 
	 * author weil
	 */
	public double gapLeftCell(int x, int y, int z) {

		if (x == 0)// 最左边的网格，没有左边的网格了！
			return -1;

		// 左边的小网格已经标记为“已搜索”，则返回-1
		if (cells[x - 1][y][z].getAlreadysearched() != 0) {
			return -1;
		}

		// 计算两个小网格的梯度差值
		double xd = 0;// x方向偏导数差值的平方

		double yd = 0;

		double zd = 0;

		double gap = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXdp() - cells[x - 1][y][z].getXdp(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYdp() - cells[x - 1][y][z].getYdp(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZdp() - cells[x - 1][y][z].getZdp(), 2.0);

		gap = StrictMath.pow(xd + yd + zd, 1.0 / 2);

		return gap;
	}
	

	/**
	 * 函数gapRightCell 
	 * 计算当前小网格与其右边小网格的差值 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标
	 * z 当前小网格的z坐标 
	 * return 当前小网格与其右边小网格的差值 
	 * author weil
	 */
	public double gapRightCell(int x, int y, int z) {

		if (x == K - 1)// 最右边的网格，没有右边的网格了！
			return -1;

		// 右边的小网格已经标记为“已搜索”，则返回-1
		if (cells[x + 1][y][z].getAlreadysearched() != 0) {
			return -1;
		}

		// 计算两个小网格的梯度差值
		double xd = 0;// x方向偏导数差值的平方

		double yd = 0;

		double zd = 0;

		double gap = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXdp() - cells[x + 1][y][z].getXdp(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYdp() - cells[x + 1][y][z].getYdp(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZdp() - cells[x + 1][y][z].getZdp(), 2.0);

		gap = StrictMath.pow(xd + yd + zd, 1.0 / 2);

		return gap;
	}
	

	/**
	 * 函数gapUpCell 
	 * 计算当前小网格与其上边小网格的差值 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标
	 * return 当前小网格与其上边小网格的差值 
	 * author weil
	 */
	public double gapUpCell(int x, int y, int z) {

		if (y == 0)// 最上边的网格，没有上边的网格了！
			return -1;

		// 上边的小网格已经标记为“已搜索”，则返回-1
		if (cells[x][y - 1][z].getAlreadysearched() != 0) {
			return -1;
		}

		// 计算两个小网格的梯度差值
		double xd = 0;// x方向偏导数差值的平方

		double yd = 0;

		double zd = 0;

		double gap = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXdp() - cells[x][y - 1][z].getXdp(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYdp() - cells[x][y - 1][z].getYdp(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZdp() - cells[x][y - 1][z].getZdp(), 2.0);

		gap = StrictMath.pow(xd + yd + zd, 1.0 / 2);

		return gap;
	}
	

	/**
	 * 函数gapDownCell 
	 * 计算当前小网格与其下边小网格的差值 
	 * Parameters: x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return 当前小网格与其下边小网格的差值 
	 * author weil
	 */
	public double gapDownCell(int x, int y, int z) {

		if (y == K - 1)// 最前边的网格，没有前边的网格了！
			return -1;

		// 前边的小网格已经标记为“已搜索”，则返回-1
		if (cells[x][y + 1][z].getAlreadysearched() != 0) {
			return -1;
		}

		// 计算两个小网格的梯度差值
		double xd = 0;// x方向偏导数差值的平方

		double yd = 0;

		double zd = 0;

		double gap = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXdp() - cells[x][y + 1][z].getXdp(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYdp() - cells[x][y + 1][z].getYdp(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZdp() - cells[x][y + 1][z].getZdp(), 2.0);

		gap = StrictMath.pow(xd + yd + zd, 1.0 / 2);

		return gap;
	}
	

	/**
	 * 函数gapFrontCell 
	 * 计算当前小网格与其前边小网格的差值 
	 * Parameters: x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标
	 * z 当前小网格的z坐标 
	 * return 当前小网格与其前边小网格的差值 
	 * author weil
	 */
	public double gapFrontCell(int x, int y, int z) {

		if (z == 0)// 最前边的网格，没有前边的网格了！
			return -1;

		// 前边的小网格已经标记为“已搜索”，则返回-1
		if (cells[x][y][z - 1].getAlreadysearched() != 0) {
			return -1;
		}

		// 计算两个小网格的梯度差值
		double xd = 0;// x方向偏导数差值的平方

		double yd = 0;

		double zd = 0;

		double gap = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXdp() - cells[x][y][z - 1].getXdp(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYdp() - cells[x][y][z - 1].getYdp(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZdp() - cells[x][y][z - 1].getZdp(), 2.0);

		gap = StrictMath.pow(xd + yd + zd, 1.0 / 2);

		return gap;
	}
	

	/**
	 * 函数gapBackCell 
	 * 计算当前小网格与其后边小网格的差值 
	 * Parameters: x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return 当前小网格与其后边小网格的差值 
	 * author weil
	 */
	public double gapBackCell(int x, int y, int z) {

		if (z == K - 1)// 最后边的网格，没有后边的网格了！
			return -1;

		// 后边的小网格已经标记为“已搜索”，则返回-1
		if (cells[x][y][z + 1].getAlreadysearched() != 0) {
			return -1;
		}

		// 计算两个小网格的梯度差值
		double xd = 0;// x方向偏导数差值的平方

		double yd = 0;

		double zd = 0;

		double gap = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXdp() - cells[x][y][z + 1].getXdp(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYdp() - cells[x][y][z + 1].getYdp(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZdp() - cells[x][y][z + 1].getZdp(), 2.0);

		gap = StrictMath.pow(xd + yd + zd, 1.0 / 2);

		return gap;
	}

	
	/**
	 * 函数partialDerivative() 
	 * 计算每个小网格的x偏导，y偏导，z偏导 
	 * return NULL
	 * author weil
	 */
	public void partialDerivative() {

		// 非空小网格的坐标集
		ArrayList<Cells> notEmptyGridCoord = new ArrayList<Cells>();
		for (int y = 0; y < K; y++) {
			for (int z = 0; z < K; z++) {
				for (int x = 0; x < K; x++) {
					if (cells[x][y][z].pNum >= 1)
						notEmptyGridCoord.add(new Cells(x, y, z));
				}
			}
		}

		// 以下输出小网格自适应性密度函数至文本文件GridWeight.txt
		String s = new String();
		try {
			File fo = new File(
					"D:\\Graduation Project\\src\\0331\\GridWeight.txt");

			if (fo.exists()) {
				System.out.println("文件GridWeight存在,删除以前的...");
				fo.delete();
			}
			if (fo.createNewFile()) {
				System.out.println("文件GridWeight创建成功！");
			} else {
				System.out.println("文件GridWeight创建失败！");
			}
			BufferedWriter outputfo = new BufferedWriter(new FileWriter(fo));

			// 计算x偏导，y偏导，z偏导
			for (int x = 0; x < K; x++) {
				for (int y = 0; y < K; y++) {
					for (int z = 0; z < K; z++) {

						// 小网格的中心坐标
						double xloc = cells[x][y][z].getXloc();
						double yloc = cells[x][y][z].getYloc();
						double zloc = cells[x][y][z].getZloc();

						// 小网格对应的二维图像的中心坐标值
						double spacialX = cells[x][y][z].getSpacialX();
						double spacialY = cells[x][y][z].getSpacialY();

						// x、y、z方向的偏导数
						double dpx = 0;
						double dpy = 0;
						double dpz = 0;
						double potential = 0;
						double gridWeight = 0;

						// 按数据场势函数的定义，计算各小网格的势值和偏导势
						for (int i = 0; i < K; i++) {
							for (int j = 0; j < K; j++) {
								for (int k = 0; k < K; k++) {
									Cells c = cells[i][j][k];
									double dx = (c.getXloc() - xloc) / width;
									double dy = (c.getYloc() - yloc) / width;
									double dz = (c.getZloc() - zloc) / width;

									double sigSpacial1 = line_num / width;
									double sigSpacial2 = row_num / width;
									double spacialDx = (c.getSpacialX() - spacialX)
											/ sigSpacial1;
									double spacialDy = (c.getSpacialY() - spacialY)
											/ sigSpacial2;

									// 特征距离
									double dist2 = dx * dx + dy * dy + dz * dz;

									// 空间距离
									double spacialdist2 = spacialDx * spacialDx
											+ spacialDy * spacialDy;

									double temp = Math.pow(Math.E,
											(double) (-dist2))
											* Math.pow(Math.E,
													(double) (-spacialdist2));

									gridWeight += Math.pow(Math.E,
											(double) (-spacialdist2));

									potential += temp;
									dpx += dx * temp;
									dpy += dy * temp;
									dpz += dz * temp;
								}
							}
						}

						s = gridWeight + "\r\n";
						outputfo.write(s);

						// 设定偏导数值
						cells[x][y][z].setPotential(potential);
						cells[x][y][z].setXdp(dpx);
						cells[x][y][z].setYdp(dpy);
						cells[x][y][z].setZdp(dpz);
						//System.out.println("x方向偏导值 dpx = "+dpx);
						//System.out.println("y方向偏导值 dpy = "+dpy);
						//System.out.println("z方向偏导值 dpz = "+dpz);
					}
				}
			}

			outputfo.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	//以下是计算当前小网格与周围网格的距离
	
	/**
	 * 函数disLeftCell 
	 * 计算当前小网格与其左边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftCell(int x, int y, int z) {

		if (x == 0)// 最左边的网格，没有左边的网格了！
			return 0;

		double d = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;
		double Z = 0;//返回值

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y][z].getZloc(), 2.0);

		d = StrictMath.pow((xd + yd + zd), 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y][z] - 0.0) != 0){
			varFactor[x - 1][y][z] = ro	* StrictMath.pow(geoMean / initial[x - 1][y][z], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y][z]);
			
			if (d <= varFactor[x - 1][y][z] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y][z] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y][z],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y][z] = " + varBandwidth[x - 1][y][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y][z].getXloc()) 
										/ varBandwidth[x - 1][y][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y][z].getYloc())
										/ varBandwidth[x - 1][y][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y][z].getZloc())
										/ varBandwidth[x - 1][y][z]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y][z] = StrictMath.pow(varBandwidth[x - 1][y][z], 3.0);
				Z = (initial[x - 1][y][z] * kernel) / varHd[x - 1][y][z];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y][z]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disRightCell 
	 * 计算当前小网格与其右边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightCell(int x, int y, int z) {

		if (x == K - 1) // 最右边的网格，没有右边的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y][z].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y][z] - 0.0) != 0){
			varFactor[x + 1][y][z] = ro	* StrictMath.pow(geoMean / initial[x + 1][y][z], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y][z]);
			
			if (d <= varFactor[x + 1][y][z] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y][z] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y][z],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y][z] = " + varBandwidth[x + 1][y][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y][z].getXloc()) 
										/ varBandwidth[x + 1][y][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y][z].getYloc())
										/ varBandwidth[x + 1][y][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y][z].getZloc())
										/ varBandwidth[x + 1][y][z]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y][z] = StrictMath.pow(varBandwidth[x + 1][y][z], 3.0);
				Z = (initial[x + 1][y][z] * kernel) / varHd[x + 1][y][z];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y][z]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disUpCell 
	 * 计算当前小网格与其上边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disUpCell(int x, int y, int z) {

		if (y == 0)// 最上边的网格，没有上边的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x][y - 1][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x][y - 1][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x][y - 1][z].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x][y - 1][z] - 0.0) != 0){
			varFactor[x][y - 1][z] = ro	* StrictMath.pow(geoMean / initial[x][y - 1][z], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x][y - 1][z]);
			
			if (d <= varFactor[x][y - 1][z] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x][y - 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x][y - 1][z],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x][y - 1][z] = " + varBandwidth[x][y - 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y - 1][z].getXloc()) 
										/ varBandwidth[x][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y - 1][z].getYloc())
										/ varBandwidth[x][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y - 1][z].getZloc())
										/ varBandwidth[x][y - 1][z]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x][y - 1][z] = StrictMath.pow(varBandwidth[x][y - 1][z], 3.0);
				Z = (initial[x][y - 1][z] * kernel) / varHd[x][y - 1][z];
				
				System.out.println(" h(x)的d次方 : " + varHd[x][y - 1][z]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disDownCell 
	 * 计算当前小网格与其下边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disDownCell(int x, int y, int z) {

		if (y == K - 1)// 最下边的网格，没有下边的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x][y + 1][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x][y + 1][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x][y + 1][z].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x][y + 1][z] - 0.0) != 0){
			varFactor[x][y + 1][z] = ro	* StrictMath.pow(geoMean / initial[x][y + 1][z], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x][y + 1][z]);
			
			if (d <= varFactor[x][y + 1][z] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x][y + 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x][y + 1][z],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x][y + 1][z] = " + varBandwidth[x][y + 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y + 1][z].getXloc()) 
										/ varBandwidth[x][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y + 1][z].getYloc())
										/ varBandwidth[x][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y + 1][z].getZloc())
										/ varBandwidth[x][y + 1][z]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x][y + 1][z] = StrictMath.pow(varBandwidth[x][y + 1][z], 3.0);
				Z = (initial[x][y + 1][z] * kernel) / varHd[x][y + 1][z];
				
				System.out.println(" h(x)的d次方 : " + varHd[x][y + 1][z]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disFrontCell 
	 * 计算当前小网格与其前边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disFrontCell(int x, int y, int z) {

		if (z == 0)// 最前边的网格，没有前边的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x][y][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x][y][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x][y][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x][y][z - 1] - 0.0) != 0){
			varFactor[x][y][z - 1] = ro	* StrictMath.pow(geoMean / initial[x][y][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x][y][z - 1]);
			
			if (d <= varFactor[x][y][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x][y][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x][y][z - 1] = " + varBandwidth[x][y][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y][z - 1].getXloc()) 
										/ varBandwidth[x][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y][z - 1].getYloc())
										/ varBandwidth[x][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y][z - 1].getZloc())
										/ varBandwidth[x][y][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x][y][z - 1] = StrictMath.pow(varBandwidth[x][y][z - 1], 3.0);
				Z = (initial[x][y][z - 1] * kernel) / varHd[x][y][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x][y][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disBackCell 
	 * 计算当前小网格与其后边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disBackCell(int x, int y, int z) {

		if (z == K - 1)// 最前边的网格，没有前边的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x][y][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x][y][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x][y][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x][y][z + 1] - 0.0) != 0){
			varFactor[x][y][z + 1] = ro	* StrictMath.pow(geoMean / initial[x][y][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x][y][z + 1]);
			
			if (d <= varFactor[x][y][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x][y][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x][y][z + 1] = " + varBandwidth[x][y][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y][z + 1].getXloc()) 
										/ varBandwidth[x][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y][z + 1].getYloc())
										/ varBandwidth[x][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y][z + 1].getZloc())
										/ varBandwidth[x][y][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x][y][z + 1] = StrictMath.pow(varBandwidth[x][y][z + 1], 3.0);
				Z = (initial[x][y][z + 1] * kernel) / varHd[x][y][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x][y][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disLeftUpCell 
	 * 计算当前小网格与其左上边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftUpCell(int x, int y, int z) {

		if (x == 0 || y == 0)// 最左上的网格，没有左上的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y - 1][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y - 1][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y - 1][z].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y - 1][z] - 0.0) != 0){
			varFactor[x - 1][y - 1][z] = ro	* StrictMath.pow(geoMean / initial[x - 1][y - 1][z], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y - 1][z]);
			
			if (d <= varFactor[x - 1][y - 1][z] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y - 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y - 1][z],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y - 1][z] = " + varBandwidth[x - 1][y - 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y - 1][z].getXloc()) 
										/ varBandwidth[x - 1][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y - 1][z].getYloc())
										/ varBandwidth[x - 1][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y - 1][z].getZloc())
										/ varBandwidth[x - 1][y - 1][z]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y - 1][z] = StrictMath.pow(varBandwidth[x - 1][y - 1][z], 3.0);
				Z = (initial[x - 1][y - 1][z] * kernel) / varHd[x - 1][y - 1][z];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y - 1][z]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disLeftDownCell 
	 * 计算当前小网格与其左上边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftDownCell(int x, int y, int z) {

		if (x == 0 || y == K - 1)// 最左下的网格，没有左下的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y + 1][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y + 1][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y + 1][z].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y + 1][z] - 0.0) != 0){
			varFactor[x - 1][y + 1][z] = ro	* StrictMath.pow(geoMean / initial[x - 1][y + 1][z], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y + 1][z]);
			
			if (d <= varFactor[x - 1][y + 1][z] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y + 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y + 1][z],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y + 1][z] = " + varBandwidth[x - 1][y + 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y + 1][z].getXloc()) 
										/ varBandwidth[x - 1][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y + 1][z].getYloc())
										/ varBandwidth[x - 1][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y + 1][z].getZloc())
										/ varBandwidth[x - 1][y + 1][z]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y + 1][z] = StrictMath.pow(varBandwidth[x - 1][y + 1][z], 3.0);
				Z = (initial[x - 1][y + 1][z] * kernel) / varHd[x - 1][y + 1][z];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y + 1][z]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disLeftFrontCell 
	 * 计算当前小网格与其左前边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftFrontCell(int x, int y, int z) {

		if (x == 0 || z == 0)// 最左前的网格，没有左前的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y][z - 1] - 0.0) != 0){
			varFactor[x - 1][y][z - 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y][z - 1]);
			
			if (d <= varFactor[x - 1][y][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y][z - 1] = " + varBandwidth[x - 1][y][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y][z - 1].getXloc()) 
										/ varBandwidth[x - 1][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y][z - 1].getYloc())
										/ varBandwidth[x - 1][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y][z - 1].getZloc())
										/ varBandwidth[x - 1][y][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y][z - 1] = StrictMath.pow(varBandwidth[x - 1][y][z - 1], 3.0);
				Z = (initial[x - 1][y][z - 1] * kernel) / varHd[x - 1][y][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disLeftBackCell 
	 * 计算当前小网格与其左后边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftBackCell(int x, int y, int z) {

		if (x == 0 || z == K - 1)// 最左后的网格，没有左后的网格了！
			return 0;

		double d = 0;
		double Z = 0;	
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y][z + 1] - 0.0) != 0){
			varFactor[x - 1][y][z + 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y][z + 1]);
			
			if (d <= varFactor[x - 1][y][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y][z + 1] = " + varBandwidth[x - 1][y][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y][z + 1].getXloc()) 
										/ varBandwidth[x - 1][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y][z + 1].getYloc())
										/ varBandwidth[x - 1][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y][z + 1].getZloc())
										/ varBandwidth[x - 1][y][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y][z + 1] = StrictMath.pow(varBandwidth[x - 1][y][z + 1], 3.0);
				Z = (initial[x - 1][y][z + 1] * kernel) / varHd[x - 1][y][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disRightUpCell 
	 * 计算当前小网格与其右上边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightUpCell(int x, int y, int z) {

		if (x == K - 1 || y == 0)// 最右上的网格，没有右上的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y - 1][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y - 1][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y - 1][z].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y - 1][z] - 0.0) != 0){
			varFactor[x + 1][y - 1][z] = ro	* StrictMath.pow(geoMean / initial[x + 1][y - 1][z], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y - 1][z]);
			
			if (d <= varFactor[x + 1][y - 1][z] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y - 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y - 1][z],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y - 1][z] = " + varBandwidth[x + 1][y - 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y - 1][z].getXloc()) 
										/ varBandwidth[x + 1][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y - 1][z].getYloc())
										/ varBandwidth[x + 1][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y - 1][z].getZloc())
										/ varBandwidth[x + 1][y - 1][z]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y - 1][z] = StrictMath.pow(varBandwidth[x + 1][y - 1][z], 3.0);
				Z = (initial[x + 1][y - 1][z] * kernel) / varHd[x + 1][y - 1][z];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y - 1][z]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disRightDownCell 
	 * 计算当前小网格与其右下边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightDownCell(int x, int y, int z) {

		if (x == K - 1 || y == K - 1)// 最右下的网格，没有右下的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y + 1][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y + 1][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y + 1][z].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y + 1][z] - 0.0) != 0){
			varFactor[x + 1][y + 1][z] = ro	* StrictMath.pow(geoMean / initial[x + 1][y + 1][z], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y + 1][z]);
			
			if (d <= varFactor[x + 1][y + 1][z] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y + 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y + 1][z],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y + 1][z] = " + varBandwidth[x + 1][y + 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y + 1][z].getXloc()) 
										/ varBandwidth[x + 1][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y + 1][z].getYloc())
										/ varBandwidth[x + 1][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y + 1][z].getZloc())
										/ varBandwidth[x + 1][y + 1][z]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y + 1][z] = StrictMath.pow(varBandwidth[x + 1][y + 1][z], 3.0);
				Z = (initial[x + 1][y + 1][z] * kernel) / varHd[x + 1][y + 1][z];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y + 1][z]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disRightFrontCell 
	 * 计算当前小网格与其右前边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightFrontCell(int x, int y, int z) {

		if (x == K - 1 || z == 0)// 最右前的网格，没有右前的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y][z - 1] - 0.0) != 0){
			varFactor[x + 1][y][z - 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y][z - 1]);
			
			if (d <= varFactor[x + 1][y][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y][z - 1] = " + varBandwidth[x + 1][y][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y][z - 1].getXloc()) 
										/ varBandwidth[x + 1][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y][z - 1].getYloc())
										/ varBandwidth[x + 1][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y][z - 1].getZloc())
										/ varBandwidth[x + 1][y][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y][z - 1] = StrictMath.pow(varBandwidth[x + 1][y][z - 1], 3.0);
				Z = (initial[x + 1][y][z - 1] * kernel) / varHd[x + 1][y][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disRightBackCell 
	 * 计算当前小网格与其右后边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightBackCell(int x, int y, int z) {

		if (x == K - 1 || z == K- 1)// 最右后的网格，没有右后的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y][z + 1] - 0.0) != 0){
			varFactor[x + 1][y][z + 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y][z + 1]);
			
			if (d <= varFactor[x + 1][y][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y][z + 1] = " + varBandwidth[x + 1][y][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y][z + 1].getXloc()) 
										/ varBandwidth[x + 1][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y][z + 1].getYloc())
										/ varBandwidth[x + 1][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y][z + 1].getZloc())
										/ varBandwidth[x + 1][y][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y][z + 1] = StrictMath.pow(varBandwidth[x + 1][y][z + 1], 3.0);
				Z = (initial[x + 1][y][z + 1] * kernel) / varHd[x + 1][y][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	
	/**
	 * 函数disUpFrontCell 
	 * 计算当前小网格与其上前边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disUpFrontCell(int x, int y, int z) {

		if (y == 0 || z == 0)// 最上前的网格，没有上前的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x][y - 1][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x][y - 1][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x][y - 1][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x][y - 1][z - 1] - 0.0) != 0){
			varFactor[x][y - 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x][y - 1][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x][y - 1][z - 1]);
			
			if (d <= varFactor[x][y - 1][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x][y - 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y - 1][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x][y - 1][z - 1] = " + varBandwidth[x][y - 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y - 1][z - 1].getXloc()) 
										/ varBandwidth[x][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y - 1][z - 1].getYloc())
										/ varBandwidth[x][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y - 1][z - 1].getZloc())
										/ varBandwidth[x][y - 1][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x][y - 1][z - 1] = StrictMath.pow(varBandwidth[x][y - 1][z - 1], 3.0);
				Z = (initial[x][y - 1][z - 1] * kernel) / varHd[x][y - 1][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x][y - 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disUpBackCell 
	 * 计算当前小网格与其上后边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disUpBackCell(int x, int y, int z) {

		if (y == 0 || z == K - 1)// 最上后的网格，没有上后的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x][y - 1][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x][y - 1][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x][y - 1][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x][y - 1][z + 1] - 0.0) != 0){
			varFactor[x][y - 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x][y - 1][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x][y - 1][z + 1]);
			
			if (d <= varFactor[x][y - 1][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x][y - 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y - 1][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x][y - 1][z + 1] = " + varBandwidth[x][y - 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y - 1][z + 1].getXloc()) 
										/ varBandwidth[x][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y - 1][z + 1].getYloc())
										/ varBandwidth[x][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y - 1][z + 1].getZloc())
										/ varBandwidth[x][y - 1][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x][y - 1][z + 1] = StrictMath.pow(varBandwidth[x][y - 1][z + 1], 3.0);
				Z = (initial[x][y - 1][z + 1] * kernel) / varHd[x][y - 1][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x][y - 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disDownFrontCell 
	 * 计算当前小网格与其下前边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disDownFrontCell(int x, int y, int z) {

		if (y == K - 1 || z == 0)// 最下前的网格，没有下前的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x][y + 1][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x][y + 1][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x][y + 1][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x][y + 1][z - 1] - 0.0) != 0){
			varFactor[x][y + 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x][y + 1][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x][y + 1][z - 1]);
			
			if (d <= varFactor[x][y + 1][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x][y + 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y + 1][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x][y + 1][z - 1] = " + varBandwidth[x][y + 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y + 1][z - 1].getXloc()) 
										/ varBandwidth[x][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y + 1][z - 1].getYloc())
										/ varBandwidth[x][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y + 1][z - 1].getZloc())
										/ varBandwidth[x][y + 1][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x][y + 1][z - 1] = StrictMath.pow(varBandwidth[x][y + 1][z - 1], 3.0);
				Z = (initial[x][y + 1][z - 1] * kernel) / varHd[x][y + 1][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x][y + 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disDownBackCell 
	 * 计算当前小网格与其下后边小网格的自适应密度  
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disDownBackCell(int x, int y, int z) {

		if (y == K - 1 || z == K - 1)// 最下后的网格，没有下后的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x][y + 1][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x][y + 1][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x][y + 1][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x][y + 1][z + 1] - 0.0) != 0){
			varFactor[x][y + 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x][y + 1][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x][y + 1][z + 1]);
			
			if (d <= varFactor[x][y + 1][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x][y + 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y + 1][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x][y + 1][z + 1] = " + varBandwidth[x][y + 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y + 1][z + 1].getXloc()) 
										/ varBandwidth[x][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y + 1][z + 1].getYloc())
										/ varBandwidth[x][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y + 1][z + 1].getZloc())
										/ varBandwidth[x][y + 1][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x][y + 1][z + 1] = StrictMath.pow(varBandwidth[x][y + 1][z + 1], 3.0);
				Z = (initial[x][y + 1][z + 1] * kernel) / varHd[x][y + 1][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x][y + 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disLeftUpFrontCell 
	 * 计算当前小网格与其左上前边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftUpFrontCell(int x, int y, int z) {

		if (x == 0 || y == 0 || z == 0)// 最左上前的网格，没有左上前的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y - 1][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y - 1][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y - 1][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y - 1][z - 1] - 0.0) != 0){
			varFactor[x - 1][y - 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y - 1][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y - 1][z - 1]);
			
			if (d <= varFactor[x - 1][y - 1][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y - 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y - 1][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y - 1][z - 1] = " + varBandwidth[x - 1][y - 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y - 1][z - 1].getXloc()) 
										/ varBandwidth[x - 1][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y - 1][z - 1].getYloc())
										/ varBandwidth[x - 1][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y - 1][z - 1].getZloc())
										/ varBandwidth[x - 1][y - 1][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y - 1][z - 1] = StrictMath.pow(varBandwidth[x - 1][y - 1][z - 1], 3.0);
				Z = (initial[x - 1][y - 1][z - 1] * kernel) / varHd[x - 1][y - 1][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y - 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disLeftUpBackCell 
	 * 计算当前小网格与其左上后边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftUpBackCell(int x, int y, int z) {

		if (x == 0 || y == 0 || z == K - 1)// 最左上后的网格，没有左上后的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y - 1][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y - 1][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y - 1][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y - 1][z + 1] - 0.0) != 0){
			varFactor[x - 1][y - 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y - 1][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y - 1][z + 1]);
			
			if (d <= varFactor[x - 1][y - 1][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y - 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y - 1][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y - 1][z + 1] = " + varBandwidth[x - 1][y - 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y - 1][z + 1].getXloc()) 
										/ varBandwidth[x - 1][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y - 1][z + 1].getYloc())
										/ varBandwidth[x - 1][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y - 1][z + 1].getZloc())
										/ varBandwidth[x - 1][y - 1][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y - 1][z + 1] = StrictMath.pow(varBandwidth[x - 1][y - 1][z + 1], 3.0);
				Z = (initial[x - 1][y - 1][z + 1] * kernel) / varHd[x - 1][y - 1][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y - 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disLeftDownFrontCell 
	 * 计算当前小网格与其左下前边小网格的自适应密度  
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftDownFrontCell(int x, int y, int z) {

		if (x == 0 || y == K - 1 || z == 0)// 最左下前的网格，没有左下前的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y + 1][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y + 1][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y + 1][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y + 1][z - 1] - 0.0) != 0){
			varFactor[x - 1][y + 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y + 1][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y + 1][z - 1]);
			
			if (d <= varFactor[x - 1][y + 1][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y + 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y + 1][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y + 1][z - 1] = " + varBandwidth[x - 1][y + 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y + 1][z - 1].getXloc()) 
										/ varBandwidth[x - 1][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y + 1][z - 1].getYloc())
										/ varBandwidth[x - 1][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y + 1][z - 1].getZloc())
										/ varBandwidth[x - 1][y + 1][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y + 1][z - 1] = StrictMath.pow(varBandwidth[x - 1][y + 1][z - 1], 3.0);
				Z = (initial[x - 1][y + 1][z - 1] * kernel) / varHd[x - 1][y + 1][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y + 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	
	/**
	 * 函数disLeftDownBackCell 
	 * 计算当前小网格与其左下后边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disLeftDownBackCell(int x, int y, int z) {

		if (x == 0 || y == K - 1 || z == K - 1)// 最左下后的网格，没有左下后的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y + 1][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y + 1][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y + 1][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x - 1][y + 1][z + 1] - 0.0) != 0){
			varFactor[x - 1][y + 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y + 1][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x - 1][y + 1][z + 1]);
			
			if (d <= varFactor[x - 1][y + 1][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x - 1][y + 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y + 1][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x - 1][y + 1][z + 1] = " + varBandwidth[x - 1][y + 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y + 1][z + 1].getXloc()) 
										/ varBandwidth[x - 1][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y + 1][z + 1].getYloc())
										/ varBandwidth[x - 1][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y + 1][z + 1].getZloc())
										/ varBandwidth[x - 1][y + 1][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x - 1][y + 1][z + 1] = StrictMath.pow(varBandwidth[x - 1][y + 1][z + 1], 3.0);
				Z = (initial[x - 1][y + 1][z + 1] * kernel) / varHd[x - 1][y + 1][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x - 1][y + 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	
	
	/**
	 * 函数disRightUpFrontCell 
	 * 计算当前小网格与其右上前边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightUpFrontCell(int x, int y, int z) {

		if (x == K - 1 || y == 0 || z == 0)// 最右上前的网格，没有右上前的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y - 1][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y - 1][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y - 1][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y - 1][z - 1] - 0.0) != 0){
			varFactor[x + 1][y - 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y - 1][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y - 1][z - 1]);
			
			if (d <= varFactor[x + 1][y - 1][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y - 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y - 1][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y - 1][z - 1] = " + varBandwidth[x + 1][y - 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y - 1][z - 1].getXloc()) 
										/ varBandwidth[x + 1][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y - 1][z - 1].getYloc())
										/ varBandwidth[x + 1][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y - 1][z - 1].getZloc())
										/ varBandwidth[x + 1][y - 1][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y - 1][z - 1] = StrictMath.pow(varBandwidth[x + 1][y - 1][z - 1], 3.0);
				Z = (initial[x + 1][y - 1][z - 1] * kernel) / varHd[x + 1][y - 1][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y - 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disRightUpBackCell 
	 * 计算当前小网格与其右上后边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightUpBackCell(int x, int y, int z) {

		if (x == K - 1 || y == 0 || z == K - 1)// 最右上后的网格，没有右上后的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y - 1][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y - 1][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y - 1][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y - 1][z + 1] - 0.0) != 0){
			varFactor[x + 1][y - 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y - 1][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y - 1][z + 1]);
			
			if (d <= varFactor[x + 1][y - 1][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y - 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y - 1][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y - 1][z + 1] = " + varBandwidth[x + 1][y - 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y - 1][z + 1].getXloc()) 
										/ varBandwidth[x + 1][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y - 1][z + 1].getYloc())
										/ varBandwidth[x + 1][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y - 1][z + 1].getZloc())
										/ varBandwidth[x + 1][y - 1][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y - 1][z + 1] = StrictMath.pow(varBandwidth[x + 1][y - 1][z + 1], 3.0);
				Z = (initial[x + 1][y - 1][z + 1] * kernel) / varHd[x + 1][y - 1][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y - 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * 函数disRightDownFrontCell 
	 * 计算当前小网格与其右下前边小网格的自适应密度
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightDownFrontCell(int x, int y, int z) {

		if (x == K - 1 || y == K - 1 || z == 0)// 最右下前的网格，没有右下前的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y + 1][z - 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y + 1][z - 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y + 1][z - 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y + 1][z - 1] - 0.0) != 0){
			varFactor[x + 1][y + 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y + 1][z - 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y + 1][z - 1]);
			
			if (d <= varFactor[x + 1][y + 1][z - 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y + 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y + 1][z - 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y + 1][z - 1] = " + varBandwidth[x + 1][y + 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y + 1][z - 1].getXloc()) 
										/ varBandwidth[x + 1][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y + 1][z - 1].getYloc())
										/ varBandwidth[x + 1][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y + 1][z - 1].getZloc())
										/ varBandwidth[x + 1][y + 1][z - 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y + 1][z - 1] = StrictMath.pow(varBandwidth[x + 1][y + 1][z - 1], 3.0);
				Z = (initial[x + 1][y + 1][z - 1] * kernel) / varHd[x + 1][y + 1][z - 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y + 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	
	/**
	 * 函数disRightDownBackCell 
	 * 计算当前小网格与其右下后边小网格的自适应密度 
	 * Parameters: 
	 * x 当前小网格的x坐标(个数) 
	 * y 当前小网格的y坐标 
	 * z 当前小网格的z坐标 
	 * return Z:累加之前的结果
	 * author weil
	 */
	public double disRightDownBackCell(int x, int y, int z) {

		if (x == K - 1 || y == K - 1 || z == K - 1)// 最右下后的网格，没有右下后的网格了！
			return 0;

		double d = 0;
		double Z = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x + 1][y + 1][z + 1].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x + 1][y + 1][z + 1].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x + 1][y + 1][z + 1].getZloc(), 2.0);

		d = StrictMath.pow(xd + yd + zd, 1.0 / 2);
		
		// 计算varFactor参数
		//只计算非空小网格
		if((initial[x + 1][y + 1][z + 1] - 0.0) != 0){
			varFactor[x + 1][y + 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y + 1][z + 1], 1.0 / 2);
			System.out.println("比例因子:" + varFactor[x + 1][y + 1][z + 1]);
			
			if (d <= varFactor[x + 1][y + 1][z + 1] * width) {

				// 欧式距离平方和
				double dis = 0;
				// 可变带宽
				varBandwidth[x + 1][y + 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y + 1][z + 1],	1.0 / 2);
				System.out.println("可变带宽 varBandwidth[x + 1][y + 1][z + 1] = " + varBandwidth[x + 1][y + 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y + 1][z + 1].getXloc()) 
										/ varBandwidth[x + 1][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y + 1][z + 1].getYloc())
										/ varBandwidth[x + 1][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y + 1][z + 1].getZloc())
										/ varBandwidth[x + 1][y + 1][z + 1]), 2.0);
									
				System.out.println("欧式距离：" + dis);
				
				double kernel = 0;

				// 满足条件||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)的d次方
				varHd[x + 1][y + 1][z + 1] = StrictMath.pow(varBandwidth[x + 1][y + 1][z + 1], 3.0);
				Z = (initial[x + 1][y + 1][z + 1] * kernel) / varHd[x + 1][y + 1][z + 1];
				
				System.out.println(" h(x)的d次方 : " + varHd[x + 1][y + 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//内层if结束	
		}//外层if结束
		else {
			Z = 0;
		}
			
		return Z;
	}
	

}
