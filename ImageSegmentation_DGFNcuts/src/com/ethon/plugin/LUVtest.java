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
import com.ethon.io.DataWriter;
import com.ethon.model.DataPoint;
import com.ethon.tools.ImageGenerator;
import com.ethon.ui.DataField;

/**
 * LUVtest聚类份额算法，依据数据场理论,网格计算,矩阵运算和图片分隔对三维图像数据DataPoint[L,U,V]进行聚类分析,平滑处理,数据分块,
 * 组融合
 * 
 * @author Jinfei Yin
 * 
 */
public class LUVtest {

	// 图片的行宽和列宽
	int line_num, row_num;

	// 每个大网格的边所包含的小网格的数量
	private int IFP;

	// 小网格的边长、起点
	public double dx, dy, dz, xmin, ymin, zmin;

	// 网格划分的特征空间影响因子和坐标空间影响因子，推荐值通过matlab估计
	private double sig1, sig2, sig3, sigSpacial1, sigSpacial2;

	// 小网格类型的矩阵
	public Grid[][][] steps;

	// 小网格坐标对象的矩阵
	public Coord[][][] allCoords;

	// 程序的开始运行时间
	private long start;

	// 聚类核心算运行法结束时间
	private long end;

	// 聚类时过滤噪声的参数
	// private int thread;

	// 网格系统每边包含小网格的数量
	private int K;

	// 大网格类型的矩阵
	private BigGrid[][][] grids;

	// 所有原始图像数据点
	DataPoint[] rawPoints;

	// 所有LUV数据点
	DataPoint[] points;

	// 最终聚类数目
	private int finalClusterNum;

	// 原始图像数据点的数量
	private int allPointsNum = 0;

	// 测试用
	// 网格系统中所有非空小网格
	private ArrayList<Grid> notEmptyGrid = new ArrayList<Grid>();

	// 网格系统中所非空大网格
	private ArrayList<BigGrid> notEmptyBigGrid = new ArrayList<BigGrid>();

	// 网格系统中所有大网格
	private ArrayList<BigGrid> allBigGrid = new ArrayList<BigGrid>();

	// 该矩阵中，所有数据点按其在原图像中的坐标(x,y)记录于对应的位置，原图像为line_num*row_num，此处今后用可变参数代替
	DataPoint[][] allPix;

	// 存放平滑后图像中所有的块,也是图像按块划分的最终结果
	LinkedList<Segment> segments = new LinkedList<Segment>();

	// 块所属分组的编号
	int segmentGroupLabelIndex = 0;

	// 平滑阈值
	double smoothValve;

	// 矩阵块分隔阈值
	double nCutValve;

	// 块分割影响因子
	double matrixSig;

	// 每块生成的辅助节点个数
	int auxiliaryNum;

	// 全图片块矩阵W
	Matrix WPrior;

	// 记录所有数据的点
	private DataPoint[] points1 = null;

	/**
	 * 构造函数HastaAndSNC，接受用户输入参数K, IFP，将所有原始图像数据点映射入网格系统，并确保初始化全部小网格对象
	 * Parameters:K:网格系统每边包含小网格的数量 IFP：每个大网格的边所包含的小网格的数量 line_num:图片的行宽
	 * row_num:图片的列宽
	 * 
	 */
	public LUVtest(int K, int IFP, int line_num, int row_num, double sig1,
			double sig2, double sig3, double sigX, double sigY,
			double sigMatrix, double smoothValve, double nCutValve,
			int auxiliaryNum) {

		start = System.currentTimeMillis();
		this.IFP = IFP;
		this.K = K;
		this.line_num = line_num;
		this.row_num = row_num;

		this.sig1 = sig1;
		this.sig2 = sig2;
		this.sig3 = sig3;

		this.sigSpacial1 = sigX;
		this.sigSpacial2 = sigY;

		this.matrixSig = sigMatrix;

		this.smoothValve = smoothValve;
		this.nCutValve = nCutValve;
		this.auxiliaryNum = auxiliaryNum;

		allPix = new DataPoint[line_num + 1][row_num + 1];
		// 实例化数据数组
		rawPoints = DataBase.getInstance().getPoints();
		points1 = rawPoints;
		DataBase db = new DataBase();
		// 保存全局变量points1
		points = db.new_public_points(points1);

		// GRB值转换为LUV
		points = RGBToLUV(points);

		allPointsNum = points.length;

		xmin = Double.MAX_VALUE;
		ymin = Double.MAX_VALUE;
		zmin = Double.MAX_VALUE;
		double xmax = Double.MIN_VALUE;
		double ymax = Double.MIN_VALUE;
		double zmax = Double.MIN_VALUE;

		steps = new Grid[K][K][K];
		allCoords = new Coord[K][K][K];
		for (DataPoint p : points) {
			double x = p.getCoord_X();
			double y = p.getCoord_Y();
			double z = p.getCoord_Z();
			// 计算三维网格系统各维的起点
			if (xmin > x)
				xmin = x;
			if (ymin > y)
				ymin = y;
			if (zmin > z)
				zmin = z;
			// 计算三维网格系统各维的终点
			if (xmax < x)
				xmax = x;
			if (ymax < y)
				ymax = y;
			if (zmax < z)
				zmax = z;
		}
		// 计算三维网格系统各维的步长(窗宽)
		dx = (xmax - xmin) / K;
		dy = (ymax - ymin) / K;
		dz = (zmax - zmin) / K;

		// 影响因子的基本确定方法
		// double sigmaPre = Math.max(dx, dy);
		// sigma = Math.max(sigmaPre, dz);
		// sigma = sigma * IFP;

		// 将每个数据点分配到相应的小网格中，可能有一部分小网格为空
		for (DataPoint p : points) {
			double x = p.getCoord_X();
			double y = p.getCoord_Y();
			double z = p.getCoord_Z();
			int xl = (int) Math.floor((x - xmin) / dx);
			int yl = (int) Math.floor((y - ymin) / dy);
			int zl = (int) Math.floor((z - zmin) / dz);

			if (xl == K)
				xl = K - 1;
			if (yl == K)
				yl = K - 1;
			if (zl == K)
				zl = K - 1;

			if (steps[xl][yl][zl] == null)
				steps[xl][yl][zl] = new Grid();
			steps[xl][yl][zl].add(p);

			if (allCoords[xl][yl][zl] == null)
				allCoords[xl][yl][zl] = new Coord(xl, yl, zl);

			// 记录所有非空网格
			if (!notEmptyGrid.contains(steps[xl][yl][zl]))
				notEmptyGrid.add(steps[xl][yl][zl]);
		}

		// 确保所有小网格均被初始化
		// 统计小空网格的数量
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					double xloc = xmin + dx / 2 + dx * x;
					double yloc = ymin + dy / 2 + dy * y;
					double zloc = zmin + dz / 2 + dz * z;

					Grid grid = steps[x][y][z];
					if (grid == null) {
						grid = new Grid();
						// 小空网格的中心设为其几何中心
						grid.setCenter(xloc, yloc, zloc);
						steps[x][y][z] = grid;
					}
					// 创建小网坐标对象矩阵
					allCoords[x][y][z] = new Coord(x, y, z);
				}
			}
		}
	}

	/**
	 * 函数RGBToLUV，接受数组RGBPoints,将所有原始图像数据点的特征三维坐标值(R,G,B)按公式转换为(L,U,V),用于聚类和分割操作
	 * 
	 */
	public DataPoint[] RGBToLUV(DataPoint[] RGBPoints) {

		double L = 0;
		double U = 0;
		double V = 0;

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

			RGBPoints[i].setCoord_X(L);
			RGBPoints[i].setCoord_Y(U);
			RGBPoints[i].setCoord_Z(V);

		}

		// 以下输出数据点(L,U,V)值至文本文件
		String s = new String();
		try {
			File fo = new File("C:\\Users\\YinJF\\Work\\DataPointsLUV.txt");

			if (fo.exists()) {
				System.out.println("文件DataPointsLUV存在,删除以前的...");
				fo.delete();
			}
			if (fo.createNewFile()) {
				System.out.println("文件DataPointsLUV创建成功！");
			} else {
				System.out.println("文件DataPointsLUV创建失败！");
			}

			BufferedWriter outputfo = new BufferedWriter(new FileWriter(fo));
			for (DataPoint p : RGBPoints) {
				s = p.getCoord_X() + " " + p.getCoord_Y() + " "
						+ p.getCoord_Z() + "\r\n";
				outputfo.write(s);
			}
			outputfo.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		return RGBPoints;
	}

	/**
	 * 函数LUVToRGB，接受数组LUVPoints,将处理所得的图像数据点的特征三维坐标值(L,U,V)转换为(R,G,B),用于图像绘制
	 * 
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

		String s = new String();
		try {
			File fo = new File("C:\\Users\\YinJF\\Work\\DataPointsLUVRGBC.txt");

			if (fo.exists()) {
				System.out.println("文件DataPointsLUVRGBC存在,删除以前的...");
				fo.delete();
			}
			if (fo.createNewFile()) {
				System.out.println("文件DataPointsLUVRGBC创建成功！");
			} else {
				System.out.println("文件DataPointsLUVRGBC创建失败！");
			}

			BufferedWriter outputfo = new BufferedWriter(new FileWriter(fo));

			for (int i = 0; i < LUVPoints.length; i++) {

				s = LUVPoints[i].getCoord_X() + " " + LUVPoints[i].getCoord_Y()
						+ " " + LUVPoints[i].getCoord_Z() + " ";

				double u = LUVPoints[i].getCoord_Y()
						/ (13 * LUVPoints[i].getCoord_X()) + un;

				double v = LUVPoints[i].getCoord_Z()
						/ (13 * LUVPoints[i].getCoord_X()) + vn;

				double RGBy = 0;
				if ((LUVPoints[i].getCoord_X() + 16) / 116 > 0.206893) {
					RGBy = StrictMath.pow(
							(LUVPoints[i].getCoord_X() + 16) / 116, 3.0);
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

				s += R + " " + G + " " + B + " "
						+ LUVPoints[i].getClusterLabel() + "\r\n";
				outputfo.write(s);

				LUVPoints[i].setCoord_X(R);
				LUVPoints[i].setCoord_Y(G);
				LUVPoints[i].setCoord_Z(B);
			}
			outputfo.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		return LUVPoints;
	}

	/**
	 * 聚类逻辑处理与矩阵分块处理函数process，算法核心
	 * 
	 * 1.调用函数markCandidatePoints()，依据图像数据点在网格数据场中的偏导势大小关系寻找势值的极大点，确定所有初始聚类中心网格
	 * 2.对初始聚类中心网格进行排序，调用函数markCandidatePoints并按空间六邻域原则合并相邻的初始聚类中心，合并后各中心的初始上下文(
	 * 范围域)的重合将完全消除 3.调用函数searchforCluster(Cluster
	 * cluster)，依据合并后的聚类中心，搜索各中心完整的网格范围域，并按聚类结果设置范围域内的所有数据点的类编号
	 * 4.调用函数transformToImagePoints(toWrite, clusters, steps, allCoords,
	 * K),找回聚类过程中的所有损失点并将边界网格点按其与已知各类欧式距离最短的原则重新分配类编号
	 * 5.调用函数smoothpix(DataPoint[] points, List<Cluster>
	 * clusters)，对分类结果进行平滑，去噪，修改所有数据点的对应类编号 6.调用函数initSegments(toWrite,
	 * segments)依据数据点的类编号,初始化所有图片分块,确定每个块初始块号,上下文点集合和相邻块号的集合
	 * 7.调用函数generateAuxiliarySeg(segments,
	 * auxiliaryNum)为每个块生成指定的的额外辅助块,并改动对应的块号和邻块信息
	 * 8.调用函数matrixWPriorByAverageValue
	 * /matrixWPriorByCovariance/matrixWPriorBySegAverageValue(segments,
	 * auxiliaryNum)构造最外层分块特征矩阵W(初始矩阵)
	 * 9.调用函数matrixDeclareAndSovle(segments)依据图片块初始矩阵初始化当前的W和D
	 * ,求特征值和特征向量得到二分子块,将结果写入全局变量,并进行递归运算,以完成分块过程
	 * 10.调用函数segmentsArrange(segments,
	 * toWrite)整理统计分块结果,为图像数据点设定组号并更按分组改点的特征三维坐标
	 * 11.调用函数LUVToRGB(toWrite)特征坐标回归为(R,G,B)
	 * 12.初始化绘图API，调用函数show_DataPoint_on_Image(BufferedImage image,DataPoint
	 * point,Color color)绘图显示聚类结果
	 */
	public void process() {
		// 所有初始聚类中心网格的坐标
		ArrayList<Coord> centers = markCandidatePoints();
		// 所有初始聚类中心
		ArrayList<Cluster> initialClusters = new ArrayList<Cluster>();
		// 记录各初始聚类中是否被已合并过，-1表示尚未合并，否则表示已被合并
		ArrayList<Integer> initMerged = new ArrayList<Integer>();
		// 合并后的初始聚类中心
		LinkedList<Cluster> finnalInitList = new LinkedList<Cluster>();
		// 所有搜索后的聚类中心
		ArrayList<Cluster> clusters = new ArrayList<Cluster>();

		// 存储聚类后的所有图像数据点
		DataPoint[] toWrite;

		// 依据坐标值初始化所有初始聚类中心
		for (int x = 0, clusterNumber = 1; x < centers.size(); x++) {
			Cluster cluster = new Cluster(centers.get(x).getX(), centers.get(x)
					.getY(), centers.get(x).getZ(), clusterNumber++);
			initialClusters.add(cluster);
		}

		// 将聚类中心按照三维坐标值由小到大排序
		Cluster tempMin = new Cluster(-1, -1, -1, -1);
		for (int i = 0; i < initialClusters.size(); i++) {
			tempMin = initialClusters.get(i);
			int minIndex = i;
			for (int j = i + 1; j < initialClusters.size(); j++) {

				if (tempMin.initialClusterCenter.getX() > initialClusters
						.get(j).initialClusterCenter.getX()) {
					tempMin = initialClusters.get(j);
					minIndex = j;
				} else if (tempMin.initialClusterCenter.getX() == initialClusters
						.get(j).initialClusterCenter.getX()) {
					if (tempMin.initialClusterCenter.getY() > initialClusters
							.get(j).initialClusterCenter.getY()) {
						tempMin = initialClusters.get(j);
						minIndex = j;
					}

					if (tempMin.initialClusterCenter.getY() == initialClusters
							.get(j).initialClusterCenter.getY()) {
						if (tempMin.initialClusterCenter.getZ() > initialClusters
								.get(j).initialClusterCenter.getZ()) {
							tempMin = initialClusters.get(j);
							minIndex = j;
						}
					}
				}
			}
			initialClusters.set(minIndex, initialClusters.get(i));
			initialClusters.set(i, tempMin);
		}

		for (int i = 0; i < initialClusters.size(); i++)
			initMerged.add(-1);
		// 依据六领域合并初始聚类中心, 依据相邻网格顺序搜索先
		for (int i = 0; i < initialClusters.size(); i++) {
			// 每次合并完成后所得到的新聚类的范围域
			LinkedList<Coord> finnalInitCoordList = new LinkedList<Coord>();
			if (!initMerged.contains(i)) {
				int j = 0;
				// 合并后的聚类中心
				Cluster mergedCluster = initialClusters.get(i);
				finnalInitCoordList = initialClusters.get(i).contentList;
				while (j < initialClusters.size()) {
					Coord toVerify = initialClusters.get(j).initialClusterCenter;
					// 如果待判定的聚类中心坐标已经被另外一类的范围域所包含，则将这两类的范围域在去重的基础上进行合并
					if (!initMerged.contains(j) && i != j
							&& finnalInitCoordList.contains(toVerify)) {
						for (Coord coord : initialClusters.get(j).contentList) {
							if (!finnalInitCoordList.contains(coord)) {
								finnalInitCoordList.add(coord);
							}
						}
						// 标记成功合并的初始聚类中心，避免重复合并
						initMerged.set(j, j);
						// 若某一类并入了新的聚类中心，介于其上文(范围域)的扩大，将重头开始遍历原始聚类中心列表
						j = 0;
					} else {
						// 若遍历到该聚类中心自身，直接跳过范围与检测，简化操作
						if (i == j) {
							initMerged.set(j, j);
						}
						j++;
					}
				}
				int num = 0;
				// 确保每个网格最多被搜索一遍且最多归属一类，getAlreadysearched()返回值1表示已被搜索过，此类点将不参与步骤3
				while (num < finnalInitCoordList.size()) {
					Coord c = finnalInitCoordList.get(num);
					if (steps[c.getX()][c.getY()][c.getZ()]
							.getAlreadysearched() == 0) {
						steps[c.getX()][c.getY()][c.getZ()]
								.setAlreadysearched(1);
						num++;
					} else {
						finnalInitCoordList.remove(c);
					}
				}
				mergedCluster.setContentList(finnalInitCoordList);
				finnalInitList.add(mergedCluster);
			}
		}

		for (int i = 0; i < finnalInitList.size(); i++) {
			// 合并后的初始聚类中心
			Cluster inicluser = new Cluster(-1, -1, -1, -1);
			inicluser = finnalInitList.get(i);

			// 按合并后聚类中心的初始网格范围域搜索其完整的网格范围域
			inicluser = searchforCluster(inicluser);
			// 按搜索得到完整范围域调整各聚类的中心网格，即聚类重心
			inicluser.adjustClusterCenter();
			clusters.add(inicluser);

		}
		// 接受用户输入的密度阈值，去除噪声类
		// String input = JOptionPane.showInputDialog("请输入聚类的密度阀值",
		// Integer.toString(thread - 1 > 0 ? thread - 1 : 0));
		// thread = Integer.parseInt(input.trim());

		finalClusterNum = clusters.size();
		toWrite = new DataPoint[allPointsNum];

		// 初始化
		for (int i = 0; i < allPointsNum; i++) {
			toWrite[i] = points[i];
		}
		// 统计聚类算法处理的数据点的总量
		int toWriteS = 0;
		// 统计聚类整理过程中处理的数据点的总量
		int j = 0;

		for (int i = 0; i < clusters.size(); i++) {
			Cluster cluster = clusters.get(i);
			List<Coord> list = cluster.getCoordList();
			// 为所有网格(坐标)设置类号
			for (int coordNum = 0; coordNum < list.size(); coordNum++) {
				allCoords[list.get(coordNum).getX()][list.get(coordNum).getY()][list
						.get(coordNum).getZ()].setClusterNum(i + 1);
			}

			int clustersize = 0;
			// 按聚类结果整理所有聚类算法处理到的数据点的相关属性(L,U,V)，类编号ClusterLabel
			for (Coord c : list) {
				Grid grid = steps[c.getX()][c.getY()][c.getZ()];
				ArrayList<DataPoint> plist = grid.list;
				clustersize += plist.size();
				// 更改点类编号
				for (DataPoint p : plist) {
					p.setClusterLabel(i);
					toWrite[(p.getLineNum() - 1)] = p;
					toWriteS++;
				}
			}
			j += clustersize;
			System.out.println(" 第  " + (int) (i + 1) + " 类" + "共有 "
					+ clustersize + " 个点");
		}

		System.out.println("------共有 "
				+ j
				+ " 个点记录,"
				+ " 平滑前损失点 "
				+ (points.length - j)
				+ " 个,"
				+ " 占总数 "
				+ new DecimalFormat("0.0000")
						.format((double) (points.length - j) / points.length)
				+ " ------");

		System.out.println(" -----共有  " + toWriteS + " 个点写出------");

		// 去噪函数:1.损失点找回 2.边界网格点重新分类 3.将数据点初始化输出图像的像素点
		toWrite = DataWriter.transformToImagePoints(toWrite, clusters, steps,
				allCoords, K);

		// 对聚类后的数据点进行平滑
		toWrite = smoothpix(toWrite, clusters);

		// 初始化所有图片分块
		segments = initSegments(toWrite, segments);

		// 生成辅助节点
		generateAuxiliarySeg(segments, auxiliaryNum);

		// 构造块矩阵
		matrixWPriorByAverageValue(segments, auxiliaryNum);
		// matrixWPriorByCovariance(segments, auxiliaryNum);
		// matrixWPriorBySegAverageValue(segments, auxiliaryNum);

		// 依据图片块矩阵求特征值和特征向量得到二分子块,并进行迭代运算,以完成分块过程
		matrixDeclareAndSovle(segments);

		// 整理分块信息,为图像数据点设定组号
		toWrite = segmentsArrange(segments, toWrite);

		// 特征坐标回归为(R,G,B)
		toWrite = LUVToRGB(toWrite);

		// 激活绘图API, 准备绘图
		BufferedImage img = ImageGenerator.drawImage(null, DataBase.sLen,
				DataBase.sHei);
		Graphics g = img.getGraphics();
		end = System.currentTimeMillis();
		String time = Long.toString(end - start) + "ms";
		for (DataPoint point : toWrite) {
			ImageGenerator.show_DataPoint_on_Image(img, point, new Color(
					(int) point.getCoord_X(), (int) point.getCoord_Y(),
					(int) point.getCoord_Z()));
		}

		g.setColor(Color.black);
		g.drawString("K=" + K + ",IFP=" + IFP, 0, DataBase.sHei - 40);
		g.drawString("共" + toWrite.length + "个点，聚成" + clusters.size() + "个类",
				0, DataBase.sHei - 25);

		g.drawString(time, 0, DataBase.sHei - 10);
		g.dispose();

		g.dispose();
		DataField.updateImagePanel(img);

	}

	/**
	 * 函数smoothpix完成数据聚类的平滑去噪 Parameters:points: 聚类所搜后待平滑的数据集 clusters:
	 * 搜索后的聚类中心列表 Return 平滑后的数据集合
	 * 
	 */
	public DataPoint[] smoothpix(DataPoint[] points, List<Cluster> clusters) {

		// 初始化allPix
		for (int i = 0; i < points.length; i++) {
			allPix[(int) points[i].getLine()][(int) points[i].getRow()] = points[i];
		}

		for (int i = 0; i < points.length; i++) {
			// 记录被考察点邻域外围(边界)的所有数据点的聚类信息
			int[] neighboorCluster = new int[finalClusterNum + 1];
			// 被考察点的邻域
			LinkedList<DataPoint> neighboor = new LinkedList<DataPoint>();
			// 被考察点的邻域的大小
			int NBsize = 0;
			// 包含被考察点邻域外围(边界)数据点数量最多的那一类的编号，需要平滑的数据点将被平滑至该类(基准类)
			int NBMax = 0;

			int k = 0;
			for (int j = 0; j < neighboorCluster.length; j++) {
				neighboorCluster[j] = 0;
			}

			if (points[i].getSmoothed() == 0) {
				boolean tag;
				// 加入该点自身
				neighboor.add(points[i]);
				do {
					tag = true;
					// 用与八邻域判定，temp记录该数据点的全部相邻数据点
					// LinkedList<DataPoint> temp =
					// allEightNeighborFeilds(neighboor
					// .get(k));
					// 用与四邻域判定，temp记录该数据点的全部相邻数据点,将用于对比试验
					LinkedList<DataPoint> temp = allFourNeighborFeilds(neighboor
							.get(k));
					for (DataPoint p : temp) {
						// 把与被考察点属于同一聚类的数据点加入邻域
						if (p.getClusterLabel() == points[i].getClusterLabel()) {
							if (!neighboor.contains(p)) {
								neighboor.add(p);
							}
							// 否则将其计入邻域边界，并定位基准类
						} else {
							int clusterLabel = p.getClusterLabel();
							if (clusterLabel != -1) {
								clusterLabel++;
							} else {
								clusterLabel = 0;
							}
							neighboorCluster[clusterLabel]++;
							if (NBMax == -1) {
								NBMax++;
							}
							if (neighboorCluster[clusterLabel] > neighboorCluster[NBMax]) {
								clusterLabel--;
								NBMax = clusterLabel;
							}
						}
					}
					k++;
					NBsize = neighboor.size();
					if (k == NBsize) {
						// 邻域已经搜索完毕不再扩大，跳出循环
						tag = false;
					}
					// 邻域大小超过平滑界限40或邻域停止增长，循环终止
				} while (NBsize < smoothValve && tag == true);

				// 将需要平滑和除噪的数据点及其邻域所有点平滑至标准类，NBsize <
				// 40表示邻域低于平滑界限的噪声点，neighboor.get(0).getClusterLabel() ==
				// -1表示该邻域中均为聚类算法损失的图像数据点

				// 除噪方法:向最大临近域归并
				if (NBsize < smoothValve
						|| neighboor.get(0).getClusterLabel() == -1) {
					if (NBMax != -1) {
						for (DataPoint p : neighboor) {
							p.setClusterLabel(NBMax);
						}
						// 聚类过程中损失的噪声点设为缺省值
					} else {
						for (DataPoint p : neighboor) {
							p.setClusterLabel(NBMax);
						}
					}
				}
				// 对于不需要平滑的点，由于邻域的互通性，不需要再次参与下次判定
				else {
					for (DataPoint p : neighboor) {
						p.setSmoothed(1);
					}
				}
			}
		}
		// 数据矩阵更新
		for (int i = 0; i < points.length; i++) {
			allPix[(int) points[i].getLine()][(int) points[i].getRow()] = points[i];
		}
		// 聚类矩阵写入文件
		// DataWriter.outPutMatrixtoFile(allPix);
		// 聚类后数据点[R,G,B]写入文件
		// DataWriter.outPutPointsRGBtoFile(points);
		// System.out.println("smoothCount= " + smoothCount);
		return points;
	}

	/**
	 * 函数allEightNeighborFeilds搜索指定数据点在原图像中的八邻域 Parameters:point:指定的数据点
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
	 * 函数allFourNeighborFeilds搜索指定数据点在原图像中的四邻域 Parameters:point:指定的数据点
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
	 * 函数initSegments依据数据点的聚类结果初始化所有分块,包括各块的编号,数据点集合和邻块编号集合
	 * Parameters:points:聚类平滑后的所有数据 segments:初始块列表 Return:初始分块结果集
	 */
	public LinkedList<Segment> initSegments(DataPoint[] points,
			LinkedList<Segment> segments) {
		// 块编号
		int segmentIndex = 0;

		for (int i = 0; i < points.length; i++) {
			// 单块大小
			int sizeofOnePiece = 0;
			int k = 0;
			if (points[i].getSegmentChecked() == 0) {
				boolean tag;
				// 图像中单个块
				Segment oneSeg = new Segment(segmentIndex, -1);
				// 加入该点自身
				oneSeg.getContentList().add(points[i]);
				points[i].setSegmentChecked(1);
				points[i].setSegmentLabel(segmentIndex);
				// 上下文及边界的初始镜像
				LinkedList<DataPoint> mirrorContentList = oneSeg
						.getContentList();
				LinkedList<DataPoint> mirrorMarginList = oneSeg.getMarginList();
				do {
					tag = true;
					// 用与四邻域判定，temp记录该数据点的全部相邻数据点,及完整的块
					LinkedList<DataPoint> temp = allFourNeighborFeilds(mirrorContentList
							.get(k));
					for (DataPoint p : temp) {
						// 把与被考察点属于同一聚类(即同一块)的数据点加入邻域
						if (p.getClusterLabel() == points[i].getClusterLabel()) {
							// 执行效率
							if (p.getSegmentChecked() == 0) {
								mirrorContentList.add(p);
								p.setSegmentChecked(1);
								p.setSegmentLabel(segmentIndex);
							}
						}
						// 若某个点的四邻域(或八邻域)中有点属于其他类或者块,则该点为所属块的边界点
						else {
							if (!mirrorMarginList.contains(mirrorContentList
									.get(k)))
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
				oneSeg.setMarginList(mirrorMarginList);
				oneSeg.setNumAll();
				segments.add(oneSeg);
				segmentIndex++;
			}

		}

		// 第二次扫描，初始化每个块的相邻块信息
		for (int j = 0; j < segments.size(); j++) {
			// 边界及邻域的初始镜像
			LinkedList<Integer> mirrorNeighborSegLabelList = segments.get(j)
					.getNeighborSegLabelList();
			LinkedList<DataPoint> mirrorMarginList = segments.get(j)
					.getMarginList();

			for (DataPoint marginPoint : mirrorMarginList) {
				LinkedList<DataPoint> temp2 = allFourNeighborFeilds(marginPoint);
				for (DataPoint neighboorPoint : temp2) {
					// 把与被考察点不同块的数据点块号记录下来
					if (neighboorPoint.getSegmentLabel() != marginPoint
							.getSegmentLabel()) {
						Integer neighboorSegLab = new Integer(
								neighboorPoint.getSegmentLabel());
						if (!mirrorNeighborSegLabelList
								.contains(neighboorSegLab)
								|| mirrorNeighborSegLabelList.size() == 0) {
							mirrorNeighborSegLabelList.add(neighboorSegLab);
						}
					}
				}
			}
			segments.get(j).setNeighborSegLabelList(mirrorNeighborSegLabelList);
		}
		return segments;
	}

	/**
	 * 函数generateAuxiliarySeg为每块生成两个镜像辅助块并修改块邻域的相关信息 Parameters:oriSegments:原始分块
	 * auxNum:每块扩展后达到的新块数 Return:带辅助结点的分块结果集
	 */
	public void generateAuxiliarySeg(LinkedList<Segment> oriSegments, int auxNum) {
		// 拓展辅助结点后的块集
		LinkedList<Segment> extendedSegments = new LinkedList<Segment>();
		// 对于新的块集,块编号和邻块的编号会增长变化,其余块属性不变
		for (int i = 0; i < oriSegments.size(); i++) {
			LinkedList<Integer> extendedSegmentsNeighboors = new LinkedList<Integer>();
			extendedSegmentsNeighboors = oriSegments.get(i)
					.getNeighborSegLabelList();
			int beforeExtend = extendedSegmentsNeighboors.size();
			// 每块的新邻块由原来的邻块衍生得来
			for (int k = 0; k < beforeExtend; k++) {
				extendedSegmentsNeighboors.add(new Integer(
						extendedSegmentsNeighboors.get(k) * 3 + 1));
				extendedSegmentsNeighboors.add(new Integer(
						extendedSegmentsNeighboors.get(k) * 3 + 2));
				extendedSegmentsNeighboors.set(k, new Integer(
						extendedSegmentsNeighboors.get(k) * 3));
			}
			// 每块一分为三
			for (int j = 0; j < 3; j++) {
				Segment newSegment = new Segment(-1, -1);
				// newSegment =
				// segments.get(i);浅复制仅能得倒三个指向同一块：segments.get(i)的指针
				newSegment.setContentList(oriSegments.get(i).getContentList());
				newSegment.setGroupLabel(oriSegments.get(i).getGroupLabel());
				newSegment.setMarginList(oriSegments.get(i).getMarginList());
				newSegment.setNumAll();

				newSegment.setNeighborSegLabelList(extendedSegmentsNeighboors);
				newSegment.setSegmentLabel(i * 3 + j);
				extendedSegments.add(newSegment);
			}
		}
		segments = extendedSegments;
	}

	/**
	 * 函数matrixWPriorByAverageValue按照“数据点-数据块”均值公式计算全图片块特征矩阵W,
	 * 后续子矩阵的W的将参照该矩阵中对应位置的值进行构造 Parameters:fatherSegments:原始分块集(已进行辅助结点生成)
	 * auxiliaryNum:每块扩展后达到的新块数 Return:初始分块矩阵W(segment-point AverageValue based)
	 */
	public void matrixWPriorByAverageValue(LinkedList<Segment> fatherSegments,
			int auxiliaryNum) {
		// 矩阵的行数
		int sideLength = fatherSegments.size();
		// 奇异矩阵不可再分,将导致奇异的行(块)剔除
		if (sideLength > 1) {
			// 特征数组和特征矩阵
			double[][] arrayW = new double[sideLength][sideLength];
			// 初始全部置0,对角线上及由同一块拓展得到块间值置为1
			int index = sideLength / auxiliaryNum;
			for (int num = 0; num < index; num++) {
				for (int i = auxiliaryNum * num; i < auxiliaryNum * (num + 1); i++) {
					for (int j = auxiliaryNum * num; j < auxiliaryNum
							* (num + 1); j++) {
						arrayW[i][j] = 1;
					}
				}
			}
			// 依据邻域关系修改矩阵的值
			for (int row = 0; row < sideLength; row++) {
				// 矩阵W的值W[i,j]=W[j,i]
				double Wij = 0;
				// 与当前块相邻且属于同一组(即同一矩阵)的所有块的块号
				LinkedList<Integer> checkForNeighboors = new LinkedList<Integer>();
				// 与当前块相邻的所有块的块号
				LinkedList<Integer> allNeighboors = fatherSegments.get(row)
						.getNeighborSegLabelList();
				for (Integer segLabel : allNeighboors) {
					checkForNeighboors.add(segLabel);
				}
				// 矩阵的列索引
				int column = 0;
				// 块内的数据点集合
				LinkedList<DataPoint> checkForSeg = fatherSegments.get(row)
						.getContentList();
				// 当前块中心
				double segCenterX = fatherSegments.get(row).getXloc();
				double segCenterY = fatherSegments.get(row).getYloc();
				double segCenterZ = fatherSegments.get(row).getZloc();

				for (Integer neighboor : checkForNeighboors) {
					// 确定当前相邻块在矩阵中的列号
					for (int columnIndex = 0; columnIndex < fatherSegments
							.size(); columnIndex++) {
						if ((fatherSegments.get(columnIndex).getSegmentLabel() == (int) neighboor)) {
							column = columnIndex;
						}
					}
					// 相邻块的中心
					double segNeighboorCenterX = segments.get((int) neighboor)
							.getXloc();
					double segNeighboorCenterY = segments.get((int) neighboor)
							.getYloc();
					double segNeighboorCenterZ = segments.get((int) neighboor)
							.getZloc();
					LinkedList<DataPoint> oneNeighboorSeg = segments.get(
							(int) neighboor).getContentList();
					Wij = 0;
					// 当前块的中心与其某一相邻块中所有点的距离的平方和
					for (DataPoint neighboorPoint : oneNeighboorSeg) {

						double neighboorPointXDelta = (neighboorPoint
								.getCoord_X() - segCenterX) / matrixSig;
						double neighboorPointYDelta = (neighboorPoint
								.getCoord_Y() - segCenterY) / matrixSig;
						double neighboorPointZDelta = (neighboorPoint
								.getCoord_Z() - segCenterZ) / matrixSig;

						double exponent = neighboorPointXDelta
								* neighboorPointXDelta + neighboorPointYDelta
								* neighboorPointYDelta + neighboorPointZDelta
								* neighboorPointZDelta;
						double betweenTowSegA = Math.pow(Math.E,
								(double) (-exponent));
						Wij += betweenTowSegA;
					}
					double Wijb = Wij / oneNeighboorSeg.size();
					Wij = 0;
					// 当前块中的所有点与其某一相邻块中心的距离的平方和
					for (DataPoint segmentPoint : checkForSeg) {
						double segPointXDelta = (segmentPoint.getCoord_X() - segNeighboorCenterX)
								/ matrixSig;
						double segPointYDelta = (segmentPoint.getCoord_Y() - segNeighboorCenterY)
								/ matrixSig;
						double segPointZDelta = (segmentPoint.getCoord_Z() - segNeighboorCenterZ)
								/ matrixSig;

						double exponentInverse = segPointXDelta
								* segPointXDelta + segPointYDelta
								* segPointYDelta + segPointZDelta
								* segPointZDelta;

						double betweenTowSegB = Math.pow(Math.E,
								-exponentInverse);
						Wij += betweenTowSegB;
					}
					double Wija = Wij / checkForSeg.size();
					arrayW[row][column] = (Wija + Wijb) / 2;
					arrayW[column][row] = (Wija + Wijb) / 2;
				}
			}
			WPrior = new Matrix(arrayW);
		}
	}

	/**
	 * 函数matrixWPriorBySegAverageValue按照“数据块-数据块”均值公式计算全图片块特征矩阵W,
	 * 后续子矩阵的W的将参照该矩阵中对应位置的值进行构造 Parameters:fatherSegments:原始分块集(已进行辅助结点生成)
	 * auxiliaryNum:每块扩展后达到的新块数 Return:初始分块矩阵W(segment-segment AverageValue
	 * based)
	 */
	public void matrixWPriorBySegAverageValue(
			LinkedList<Segment> fatherSegments, int auxiliaryNum) {
		// 矩阵的行数
		int sideLength = fatherSegments.size();
		// 奇异矩阵不可再分,将导致奇异的行(块)剔除
		if (sideLength > 1) {
			// 特征数组和特征矩阵
			double[][] arrayW = new double[sideLength][sideLength];
			// 初始全部置0,对角线上及由同一块拓展得到块间值置为1
			int index = sideLength / auxiliaryNum;
			for (int num = 0; num < index; num++) {
				for (int i = auxiliaryNum * num; i < auxiliaryNum * (num + 1); i++) {
					for (int j = auxiliaryNum * num; j < auxiliaryNum
							* (num + 1); j++) {
						arrayW[i][j] = 1;
					}
				}
			}
			// 依据邻域关系修改矩阵的值
			for (int row = 0; row < sideLength; row++) {
				// 矩阵W的值W[i,j]=W[j,i]
				double Wij = 0;
				// 与当前块相邻且属于同一组(即同一矩阵)的所有块的块号
				LinkedList<Integer> checkForNeighboors = new LinkedList<Integer>();
				// 与当前块相邻的所有块的块号
				LinkedList<Integer> allNeighboors = fatherSegments.get(row)
						.getNeighborSegLabelList();
				for (Integer segLabel : allNeighboors) {
					checkForNeighboors.add(segLabel);
				}
				// 矩阵的列索引
				int column = 0;

				// 块中心
				double segCenterX = fatherSegments.get(row).getXloc();
				double segCenterY = fatherSegments.get(row).getYloc();
				double segCenterZ = fatherSegments.get(row).getZloc();

				for (Integer neighboor : checkForNeighboors) {
					// 确定当前相邻块在矩阵中的列号
					for (int columnIndex = 0; columnIndex < fatherSegments
							.size(); columnIndex++) {
						if ((fatherSegments.get(columnIndex).getSegmentLabel() == (int) neighboor)) {
							column = columnIndex;
						}
					}
					// 相邻块的中心
					double segNeighboorCenterX = segments.get((int) neighboor)
							.getXloc();
					double segNeighboorCenterY = segments.get((int) neighboor)
							.getYloc();
					double segNeighboorCenterZ = segments.get((int) neighboor)
							.getZloc();

					double DletaX = (segCenterX - segNeighboorCenterX)
							/ matrixSig;
					double DletaY = (segCenterY - segNeighboorCenterY)
							/ matrixSig;
					double DletaZ = (segCenterZ - segNeighboorCenterZ)
							/ matrixSig;

					double exponent = DletaX * DletaX + DletaY * DletaY
							+ DletaZ * DletaZ;

					Wij = Math.pow(Math.E, -exponent);
					arrayW[row][column] = Wij;
					arrayW[column][row] = Wij;
				}
			}
			WPrior = new Matrix(arrayW);
		}
	}

	/**
	 * 函数matrixWPriorBySegAverageValue按照“数据
	 * 块-数据点”协方差公式计算全图片块特征矩阵W,后续子矩阵的W的将参照该矩阵中对应位置的值进行构造
	 * Parameters:fatherSegments:原始分块集(已进行辅助结点生成) auxiliaryNum:每块扩展后达到的新块数
	 * Return:初始分块矩阵W(segment-point Covariance based)
	 */
	public void matrixWPriorByCovariance(LinkedList<Segment> fatherSegments,
			int auxiliaryNum) {
		// 矩阵的行数
		int sideLength = fatherSegments.size();
		// 三维单位矩阵
		Matrix I = new Matrix(new double[][] { { 1, 0, 0 }, { 0, 1, 0 },
				{ 0, 0, 1 } });

		if (sideLength > 1) {
			// 特征数组和特征矩阵
			double[][] arrayW = new double[sideLength][sideLength];
			// 初始全部置0,对角线上及由同一块拓展得到块间值置为1
			int index = sideLength / auxiliaryNum;
			for (int num = 0; num < index; num++) {
				for (int i = auxiliaryNum * num; i < auxiliaryNum * (num + 1); i++) {
					for (int j = auxiliaryNum * num; j < auxiliaryNum
							* (num + 1); j++) {
						arrayW[i][j] = 1;
					}
				}
			}
			// 依据邻域关系修改矩阵的值
			for (int row = 0; row < sideLength; row++) {
				// 矩阵W的值W[i,j]=W[j,i]
				double Wij = 0;
				// 协方差矩阵
				Matrix onePointCovR = new Matrix(new double[][] { { 0, 0, 0 },
						{ 0, 0, 0 }, { 0, 0, 0 } });

				Matrix CovarianceSegRow = new Matrix(new double[][] {
						{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } });

				// 与当前块相邻且属于同一组(即同一矩阵)的所有块的块号
				LinkedList<Integer> checkForNeighboors = new LinkedList<Integer>();
				// 与当前块相邻的所有块的块号
				LinkedList<Integer> allNeighboors = fatherSegments.get(row)
						.getNeighborSegLabelList();
				for (Integer segLabel : allNeighboors) {
					checkForNeighboors.add(segLabel);
				}
				// 矩阵的列索引
				int column = 0;
				// 块内的数据点集合
				LinkedList<DataPoint> checkForSeg = fatherSegments.get(row)
						.getContentList();
				// 当前块(L,U,V)均值
				double segCenterX = fatherSegments.get(row).getXloc();
				double segCenterY = fatherSegments.get(row).getYloc();
				double segCenterZ = fatherSegments.get(row).getZloc();

				// 构造当前块的协方差矩阵Matrix CovarianceSegRow
				for (DataPoint segmentPoint : checkForSeg) {
					double segPointXDelta = (segmentPoint.getCoord_X() - segCenterX)
							/ matrixSig;
					double segPointYDelta = (segmentPoint.getCoord_Y() - segCenterY)
							/ matrixSig;
					double segPointZDelta = (segmentPoint.getCoord_Z() - segCenterZ)
							/ matrixSig;

					Matrix covRow = new Matrix(new double[] { segPointXDelta,
							segPointYDelta, segPointZDelta }, 3);
					Matrix covRowT = covRow.transpose();
					onePointCovR = covRow.times(covRowT);
					CovarianceSegRow = CovarianceSegRow.plus(onePointCovR);
				}
				CovarianceSegRow = CovarianceSegRow.times((double) 1
						/ checkForSeg.size());

				// 当前块均值向量
				Matrix RowMeanVector = new Matrix(new double[] { segCenterX,
						segCenterY, segCenterZ }, 3);

				for (Integer neighboor : checkForNeighboors) {
					// 确定当前相邻块在矩阵中的列号
					for (int columnIndex = 0; columnIndex < fatherSegments
							.size(); columnIndex++) {
						if ((fatherSegments.get(columnIndex).getSegmentLabel() == (int) neighboor)) {
							column = columnIndex;
						}
					}

					Matrix onePointCovC = new Matrix(new double[][] {
							{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } });

					Matrix CovarianceSegColumn = new Matrix(new double[][] {
							{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } });

					// 构造相邻块的协方差矩阵Matrix CovarianceSegColumn
					// 相邻块的点集
					LinkedList<DataPoint> oneNeighboorSeg = segments.get(
							(int) neighboor).getContentList();
					// 相邻块的中心(L,U,V)均值
					double segNeighboorCenterX = segments.get((int) neighboor)
							.getXloc();
					double segNeighboorCenterY = segments.get((int) neighboor)
							.getYloc();
					double segNeighboorCenterZ = segments.get((int) neighboor)
							.getZloc();

					for (DataPoint neighboorPoint : oneNeighboorSeg) {
						double neighboorPointXDelta = (neighboorPoint
								.getCoord_X() - segNeighboorCenterX)
								/ matrixSig;
						double neighboorPointYDelta = (neighboorPoint
								.getCoord_Y() - segNeighboorCenterY)
								/ matrixSig;
						double neighboorPointZDelta = (neighboorPoint
								.getCoord_Z() - segNeighboorCenterZ)
								/ matrixSig;

						Matrix covCol = new Matrix(new double[] {
								neighboorPointXDelta, neighboorPointYDelta,
								neighboorPointZDelta }, 3);
						Matrix covColT = covCol.transpose();
						onePointCovC = covCol.times(covColT);
						CovarianceSegColumn.plusEquals(onePointCovC);
					}
					CovarianceSegColumn = CovarianceSegColumn.times((double) 1
							/ oneNeighboorSeg.size());

					// 相邻块均值向量
					Matrix ColMeanVector = new Matrix(new double[] {
							segNeighboorCenterX, segNeighboorCenterX,
							segNeighboorCenterX }, 3);
					// 计算Wij
					Matrix traceMatrix;
					Matrix RI_C;
					Matrix CI_R;
					Matrix differentialVector;
					Matrix differentialVectorTranspose;
					Matrix differentialRI_CI;

					RI_C = CovarianceSegRow.inverse()
							.times(CovarianceSegColumn);
					CI_R = CovarianceSegColumn.inverse()
							.times(CovarianceSegRow);
					traceMatrix = RI_C.plus(CI_R).minus(I.times(2));

					differentialVector = RowMeanVector.minus(ColMeanVector);
					differentialVectorTranspose = differentialVector
							.transpose();
					differentialRI_CI = (CovarianceSegRow.inverse())
							.plus(CovarianceSegColumn.inverse());

					Matrix x = differentialVectorTranspose.times(
							differentialRI_CI).times(differentialVector);

					double dij = 0.5 * (traceMatrix.trace() + x.get(0, 0));
					Wij = 1 - Math.pow(Math.E, (double) (-dij / 8));

					arrayW[row][column] = Wij;
					arrayW[column][row] = Wij;
				}
			}
			WPrior = new Matrix(arrayW);
		}
	}

	/**
	 * 函数matrixDeclareAndSovle按照当前父节点生成特征矩阵D和W,求特征向量和特征值,并得倒二分子节点,子节点进行递归,
	 * 当超过分块阈值时递归返回 Parameters:fatherSegments:当前父分块集 Return:递归标识
	 */
	public int matrixDeclareAndSovle(LinkedList<Segment> fatherSegments) {
		// 本次二分的左右两个子块集
		LinkedList<Segment> leftChildSegments = new LinkedList<Segment>();
		LinkedList<Segment> rightChildSegments = new LinkedList<Segment>();
		// 矩阵D,W,D-W,pow(D,-1/2)(D-W)pow(D-1/2)
		Matrix D;
		Matrix W;
		Matrix dInverseSqar;
		Matrix d_Minus_w;
		Matrix equator;

		// 矩阵的行数
		int sideLength = 0;
		sideLength = fatherSegments.size();

		// 奇异矩阵不可再分,将导致奇异的行(块)剔除
		int checkForSingular = 0;
		while (checkForSingular < fatherSegments.size()) {
			boolean noNeighboorInMatrix = true;
			// 与当前块相邻的所有块的块号
			LinkedList<Integer> allNeighboors = fatherSegments.get(
					checkForSingular).getNeighborSegLabelList();
			for (Integer segLabel : allNeighboors) {
				if (segments.get((int) segLabel).getGroupLabel() == fatherSegments
						.get(checkForSingular).getGroupLabel()) {
					noNeighboorInMatrix = false;
				}
			}
			if (noNeighboorInMatrix == true) {
				fatherSegments.remove(checkForSingular);
			} else {
				checkForSingular++;
			}
		}
		sideLength = fatherSegments.size();

		if (sideLength > 1) {
			// 特征数组和特征矩阵
			double[][] arrayD = new double[sideLength][sideLength];
			double[][] arrayW = new double[sideLength][sideLength];
			// 初始全部置0
			for (int i = 0; i < sideLength; i++)
				for (int j = 0; j < sideLength; j++) {
					arrayD[i][j] = 0;
					arrayW[i][j] = 0;
				}

			// 依据当前矩阵块的排列顺序对应查找初始W的值,构造当前矩阵的W和D
			// 矩阵的行索引
			int row;
			int column;
			for (row = 0; row < sideLength; row++) {
				// 矩阵D的对角线上的值D[k,k]
				double Dk = 0;
				// 矩阵的列索引
				for (column = 0; column < sideLength; column++) {
					arrayW[row][column] = WPrior.get(fatherSegments.get(row)
							.getSegmentLabel(), fatherSegments.get(column)
							.getSegmentLabel());

					Dk += arrayW[row][column];
				}
				arrayD[row][row] = Dk;
			}

			D = new Matrix(arrayD);
			W = new Matrix(arrayW);
			d_Minus_w = D.minus(W);
			dInverseSqar = D.inverse();

			for (int k = 0; k < sideLength; k++) {
				double value = Math.pow(dInverseSqar.get(k, k), (double) 1 / 2);
				dInverseSqar.set(k, k, value);
			}

			equator = dInverseSqar.times(d_Minus_w);
			equator = equator.times(dInverseSqar);

			Matrix eigenVectors = equator.eig().getV();

			// 特征向量=pow(D,-1/2)×eigenVectors
			eigenVectors = dInverseSqar.times(eigenVectors);

			double[] eigenValues = equator.eig().getRealEigenvalues();

			LinkedList<LinkedList<Segment>> checkForSubDivsion = computeNcutValue(
					eigenValues, eigenVectors, fatherSegments, D, W);

			// 判定当前块集是否满足继续划分的阈值条件,满足则获取子块集
			if ((fatherSegments.size() >= 2)
					&& (checkForSubDivsion.size() == 2)) {
				leftChildSegments = checkForSubDivsion.get(0);
				rightChildSegments = checkForSubDivsion.get(1);
				// 左右子块集更新组编号
				for (Segment seg : leftChildSegments) {
					seg.setGroupLabel(segmentGroupLabelIndex);
				}
				for (Segment seg : rightChildSegments) {
					int nextSegmentGroupLabelIndex = segmentGroupLabelIndex + 1;
					seg.setGroupLabel(nextSegmentGroupLabelIndex);
				}
				// 将二分结果写入全局变量segments
				for (Segment seg : segments) {
					if (leftChildSegments.contains(seg)) {
						seg.setGroupLabel(segmentGroupLabelIndex);
					} else if (rightChildSegments.contains(seg)) {
						int nextSegmentGroupLabelIndex = segmentGroupLabelIndex + 1;
						seg.setGroupLabel(nextSegmentGroupLabelIndex);
					}

				}
				// 组编号递增
				segmentGroupLabelIndex = segmentGroupLabelIndex + 2;
				// 左右子块集递归二分
				matrixDeclareAndSovle(leftChildSegments);
				matrixDeclareAndSovle(rightChildSegments);
				return 0;
			} else {
				return -1;
			}
		} else
			return -1;
	}

	/**
	 * 函数computeNcutValue每次在所有的二分策略中寻找最优划分，找不到则返回空集
	 * Parameters:eigenValues:父块矩阵特征值数组 eigenVectors:父块矩阵特征向量矩阵
	 * fatherSegments:父块集 D_of_FatherSeg:父块D矩阵 W_of_FatherSeg:父块W矩阵 Return:二分子块
	 */
	public LinkedList<LinkedList<Segment>> computeNcutValue(
			double[] eigenValues, Matrix eigenVectors,
			LinkedList<Segment> fatherSegments, Matrix D_of_FatherSeg,
			Matrix W_of_FatherSeg) {

		LinkedList<LinkedList<Segment>> result = new LinkedList<LinkedList<Segment>>();
		// 期望的最优划分
		double minNcutValue = Double.MAX_VALUE;
		// 特征值从小到大的排序结果记录对应的块号
		double[] segLabelSortByEigenValues = new double[eigenValues.length];
		// 所有可能的二分策略
		int possibleDivision = 0;
		// 第二小的特征值对应的特征向量对应的块号排列
		Matrix secondSmallestEigenVectors = new Matrix(eigenValues.length, 1);
		// 第二小的特征值对应的特征向量
		Matrix secondSmallestEigenVectorsMirror = new Matrix(
				eigenValues.length, 1);

		// 特征值排序
		for (int i = 0; i < eigenValues.length; i++) {
			segLabelSortByEigenValues[i] = eigenValues[i];
		}

		double switchEV = 0;
		for (int i = 0; i < segLabelSortByEigenValues.length; i++) {
			int minEigenValueSeg = i;
			for (int j = i + 1; j < segLabelSortByEigenValues.length; j++) {
				if (segLabelSortByEigenValues[minEigenValueSeg] > segLabelSortByEigenValues[j]) {
					minEigenValueSeg = j;
				}
			}
			switchEV = segLabelSortByEigenValues[i];
			segLabelSortByEigenValues[i] = segLabelSortByEigenValues[minEigenValueSeg];
			segLabelSortByEigenValues[minEigenValueSeg] = switchEV;
		}

		// 按其值从小到大的排序结果记录对应的块号,找到第二小的特征值对应的特征向量
		for (int i = 0; i < segLabelSortByEigenValues.length; i++) {
			// 若有多个特征值相等的情况，服从首次满足原则
			boolean tag = true;
			int eigenValueSegLabel = i;
			for (int j = 0; j < eigenValues.length; j++) {
				if (segLabelSortByEigenValues[eigenValueSegLabel] == eigenValues[j]
						&& tag) {
					tag = false;
					eigenValueSegLabel = j;
				}
			}
			segLabelSortByEigenValues[i] = eigenValueSegLabel;
		}

		secondSmallestEigenVectors = eigenVectors.getMatrix(0,
				eigenValues.length - 1, (int) segLabelSortByEigenValues[1],
				(int) segLabelSortByEigenValues[1]);

		for (int i = 0; i < secondSmallestEigenVectors.getRowDimension(); i++) {
			secondSmallestEigenVectorsMirror.set(i, 0,
					secondSmallestEigenVectors.get(i, 0));
		}

		// 按特征向量中的值从小到大的顺序排序
		double switchSMVE = 0;
		for (int i = 0; i < eigenValues.length; i++) {
			int secMinEigenVectorSeg = i;
			for (int j = i + 1; j < eigenValues.length; j++) {
				if (secondSmallestEigenVectors.get(secMinEigenVectorSeg, 0) > secondSmallestEigenVectors
						.get(j, 0)) {
					secMinEigenVectorSeg = j;
				}
			}
			switchSMVE = secondSmallestEigenVectors
					.get(secMinEigenVectorSeg, 0);
			secondSmallestEigenVectors.set(secMinEigenVectorSeg, 0,
					secondSmallestEigenVectors.get(i, 0));
			secondSmallestEigenVectors.set(i, 0, switchSMVE);
		}
		// 记录块在矩阵中的号
		for (int i = 0; i < secondSmallestEigenVectors.getRowDimension(); i++) {
			boolean tag = true;
			int eigenVectorSegLabel = i;
			for (int j = 0; j < secondSmallestEigenVectorsMirror
					.getRowDimension(); j++) {
				if ((secondSmallestEigenVectors.get(i, 0) == secondSmallestEigenVectorsMirror
						.get(j, 0)) && tag) {
					tag = false;
					eigenVectorSegLabel = j;
					secondSmallestEigenVectorsMirror
							.set(j, 0, Double.MAX_VALUE);
				}
			}
			secondSmallestEigenVectors.set(i, 0, eigenVectorSegLabel);
		}

		while (possibleDivision < eigenValues.length - 1) {
			// 某一种二分法构造的子块集,候选解
			LinkedList<Segment> leftChildSegmentsCandidate = new LinkedList<Segment>();
			LinkedList<Segment> rightChildSegmentsCandidate = new LinkedList<Segment>();
			// 单次划分的Ncut值和相关的计算中间值
			double certainNcutValue = 0;

			double cutLR = 0;
			double assocLF = 0;
			double cutRL = 0;
			// double cutRL1 = 0;
			double assocRF = 0;
			// 构造单次划分的子块集
			for (int i = 0; i < secondSmallestEigenVectors.getRowDimension(); i++) {
				if (i <= possibleDivision)
					leftChildSegmentsCandidate.add(fatherSegments
							.get((int) secondSmallestEigenVectors.get(i, 0)));
				else
					rightChildSegmentsCandidate.add(fatherSegments
							.get((int) secondSmallestEigenVectors.get(i, 0)));
			}

			for (int i = 0; i < leftChildSegmentsCandidate.size(); i++) {
				for (int j = 0; j < rightChildSegmentsCandidate.size(); j++) {
					cutLR += WPrior.get(leftChildSegmentsCandidate.get(i)
							.getSegmentLabel(), rightChildSegmentsCandidate
							.get(j).getSegmentLabel());
				}
				for (int k = 0; k < fatherSegments.size(); k++) {
					assocLF += WPrior.get(leftChildSegmentsCandidate.get(i)
							.getSegmentLabel(), fatherSegments.get(k)
							.getSegmentLabel());
				}
			}
			certainNcutValue += cutLR / assocLF;

			for (int m = 0; m < rightChildSegmentsCandidate.size(); m++) {
				for (int n = 0; n < leftChildSegmentsCandidate.size(); n++) {
					cutRL += WPrior
							.get(rightChildSegmentsCandidate.get(m)
									.getSegmentLabel(),
									leftChildSegmentsCandidate.get(n)
											.getSegmentLabel());
				}
				for (int x = 0; x < fatherSegments.size(); x++) {
					assocRF += WPrior.get(rightChildSegmentsCandidate.get(m)
							.getSegmentLabel(), fatherSegments.get(x)
							.getSegmentLabel());
				}
			}
			certainNcutValue += cutRL / assocRF;

			if (certainNcutValue < minNcutValue) {
				minNcutValue = certainNcutValue;
				if (result.size() != 0) {
					result.remove(1);
					result.remove(0);
				}
				result.add(leftChildSegmentsCandidate);
				result.add(rightChildSegmentsCandidate);
			}
			possibleDivision++;

		}

		// 二分法临界判定
		if (minNcutValue >= nCutValve) {
			result.remove(1);
			result.remove(0);
		}

		return result;
	}

	/**
	 * 函数segmentsArrange将块中的数据点按其组编号归并 Parameters:allSegmentsGrouped:所有块分组后的结果集
	 * allPoints:所有点集
	 * 
	 * Return:分块结束,修改特征三维属性值后的所有数据点集合
	 */

	public DataPoint[] segmentsArrange(LinkedList<Segment> allSegmentsGrouped,
			DataPoint[] allPoints) {
		// 分组后的数据点
		DataPoint[] classifiedPoints = new DataPoint[allPoints.length];
		LinkedList<Group> allGroups = new LinkedList<Group>();

		for (int i = 0; i < allSegmentsGrouped.size(); i++) {
			int included = -1;
			// 该块是否已在已建立的组中
			for (int j = 0; j < allGroups.size(); j++) {
				if (allGroups.get(j).getGroupLabel() == allSegmentsGrouped.get(
						i).getGroupLabel()) {
					included = j;

				}
			}
			// 不在则建立新的分组
			if (included == -1) {
				Group newGroup = new Group(allSegmentsGrouped.get(i)
						.getGroupLabel());
				// 为新的组建立数据点集合,for{...}在此处防止指针传递(浅复制)
				LinkedList<DataPoint> segentsContentMirrow = new LinkedList<DataPoint>();
				for (int all = 0; all < allSegmentsGrouped.get(i)
						.getContentList().size(); all++) {
					segentsContentMirrow.add(allSegmentsGrouped.get(i)
							.getContentList().get(all));
				}
				newGroup.setContentList(segentsContentMirrow);
				newGroup.setNumAll();
				allGroups.add(newGroup);
				// 否则加入现有分组
			} else {
				Group alreadyExistedGroup = allGroups.get(included);
				LinkedList<DataPoint> newContentList = alreadyExistedGroup
						.getContentList();
				newContentList.addAll(allSegmentsGrouped.get(i)
						.getContentList());
				alreadyExistedGroup.setContentList(newContentList);
				alreadyExistedGroup.setNumAll();
				allGroups.set(included, alreadyExistedGroup);
			}
		}

		// 依据分组结果修改数据点的特征空间值
		for (Group group : allGroups) {
			double centerL = group.getXloc();
			double centerU = group.getYloc();
			double centerV = group.getZloc();

			for (DataPoint point : group.getContentList()) {
				point.setCoord_X(centerL);
				point.setCoord_Y(centerU);
				point.setCoord_Z(centerV);
				classifiedPoints[point.getLineNum() - 1] = point;
			}
		}

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
	 * 搜索指定聚类中心所代表的聚类簇包含的所有网格(范围域) cluster：指定聚类中心 return 完整的聚类簇
	 */
	public Cluster searchforCluster(Cluster cluster) {
		// 待搜索聚类的范围域
		LinkedList<Coord> list = new LinkedList<Coord>();
		list = cluster.contentList;

		while (true) {
			int loopS = list.size();
			// 搜索空间网格的上下方向
			List<Coord> upDown = checkUpAndDown(list);
			if (upDown.size() != 0) {
				for (Coord c : upDown) {
					if (!list.contains(c))
						list.add(c);
				}
				upDown.clear();
			}
			// 搜索空间网格的左右方向
			List<Coord> leftRight = checkLeftAndRight(list);
			if (leftRight.size() != 0) {
				for (Coord c : leftRight) {
					if (!list.contains(c))
						list.add(c);
				}
				leftRight.clear();
			}
			// 搜索空间网格的前后方向
			List<Coord> frontBack = checkFrontAndBack(list);
			for (Coord c : frontBack) {
				if (!list.contains(c))
					list.add(c);
			}
			frontBack.clear();
			int loopE = list.size();
			// 当该类的范围域不再增大时，搜索结束
			if (loopE == loopS)
				break;
		}
		cluster.contentList = list;
		return cluster;
	}

	/**
	 * 在空间网格的左右方向搜索聚类所包含的所有网格 list：聚类的初始范围域 return result：左右方向完整的网格范围域
	 * 
	 */
	List<Coord> checkLeftAndRight(List<Coord> list) {
		List<Coord> result = new LinkedList<Coord>();
		for (Coord gnumber : list) {
			int x = gnumber.getX();
			int y = gnumber.getY();
			int z = gnumber.getZ();
			// 边界条件检测变量
			boolean flag = true;
			while (flag && x >= 0) {
				flag = checkLeftGrid(x, y, z, list, result);
				x--;
			}
			x = gnumber.getX();
			flag = true;
			while (flag && x < K - 1) {
				flag = checkRightGrid(x, y, z, list, result);
				x++;
			}
		}
		return result;
	}

	/**
	 * 在空间网格的上下方向搜索聚类所包含的所有网格 list：聚类的初始范围域 return result：上下方向完整的网格范围域
	 * 
	 */
	private List<Coord> checkUpAndDown(List<Coord> list) {
		List<Coord> result = new LinkedList<Coord>();
		for (Coord gnumber : list) {
			int x = gnumber.getX();
			int y = gnumber.getY();
			int z = gnumber.getZ();

			boolean flag = true;
			while (flag && y >= 0) {
				flag = checkUpGrid(x, y, z, list, result);
				y--;
			}
			y = gnumber.getZ();
			flag = true;
			while (flag && y < K - 1) {
				flag = checkDownGrid(x, y, z, list, result);
				y++;
			}
		}
		return result;
	}

	/**
	 * 在空间网格的前后方向搜索聚类所包含的所有网格 list：聚类的初始范围域 return result：前后方向完整的网格范围域
	 * 
	 */
	private List<Coord> checkFrontAndBack(List<Coord> list) {
		List<Coord> result = new LinkedList<Coord>();
		for (Coord gnumber : list) {
			int x = gnumber.getX();
			int y = gnumber.getY();
			int z = gnumber.getZ();

			boolean flag = true;
			while (flag && z >= 0) {
				flag = checkFrontGrid(x, y, z, list, result);
				z--;
			}
			z = gnumber.getZ();
			flag = true;
			while (flag && z < K - 1) {
				flag = checkBackGrid(x, y, z, list, result);
				z++;
			}
		}
		return result;
	}

	/**
	 * 检查左方的网格，如果左方的grid的偏导值梯度大于当前网格的偏导值梯度，就返回true，否则返回false x 当前网格的X坐标 y
	 * 当前网格的Y坐标 z 当前网格的Z坐标 list 当前聚类所包含的网格列表(网格范围域) return 左方的网格是否符合条件
	 */
	private boolean checkLeftGrid(int x, int y, int z, List<Coord> list,
			List<Coord> result) {
		if (x == 0)
			return false;
		if (steps[x - 1][y][z].getAlreadysearched() != 0)
			return false;
		// 计算两点的梯度值
		double d1 = steps[x][y][z].getXdp() * steps[x][y][z].getXdp()
				+ steps[x][y][z].getYdp() * steps[x][y][z].getYdp()
				+ steps[x][y][z].getZdp() * steps[x][y][z].getZdp();
		double d2 = steps[x - 1][y][z].getXdp() * steps[x - 1][y][z].getXdp()
				+ steps[x - 1][y][z].getYdp() * steps[x - 1][y][z].getYdp()
				+ steps[x - 1][y][z].getZdp() * steps[x - 1][y][z].getZdp();
		if (d1 > d2)
			return false;

		Coord coord = new Coord(x - 1, y, z);
		if (!list.contains(coord) && !result.contains(coord)) {
			result.add(coord);
			steps[x - 1][y][z].setAlreadysearched(1);
		}
		return true;
	}

	/**
	 * 检查右方的网格，如果右方的grid的偏导值梯度大于当前网格的偏导值梯度，就返回true，否则返回false x 当前网格的X坐标 y
	 * 当前网格的Y坐标 z 当前网格的Z坐标 list 当前聚类所包含的网格列表(网格范围域) return 右方的网格是否符合条件
	 */
	private boolean checkRightGrid(int x, int y, int z, List<Coord> list,
			List<Coord> result) {
		if (x == K - 1)
			return false;
		if (steps[x + 1][y][z].getAlreadysearched() != 0)
			return false;
		double d1 = steps[x][y][z].getXdp() * steps[x][y][z].getXdp()
				+ steps[x][y][z].getYdp() * steps[x][y][z].getYdp()
				+ steps[x][y][z].getZdp() * steps[x][y][z].getZdp();
		double d2 = steps[x + 1][y][z].getXdp() * steps[x + 1][y][z].getXdp()
				+ steps[x + 1][y][z].getYdp() * steps[x + 1][y][z].getYdp()
				+ steps[x + 1][y][z].getZdp() * steps[x + 1][y][z].getZdp();
		if (d1 > d2)
			return false;

		Coord coord = new Coord(x + 1, y, z);
		if (!list.contains(coord) && !result.contains(coord)) {
			result.add(coord);
			steps[x + 1][y][z].setAlreadysearched(1);
		}
		return true;
	}

	/**
	 * 检查上方的网格，如果上方的grid的偏导值梯度大于当前网格的偏导值梯度，就返回true，否则返回false x 当前网格的X坐标 y
	 * 当前网格的Y坐标 z 当前网格的Z坐标 list 当前聚类所包含的网格列表(网格范围域) return 上方的网格是否符合条件
	 */
	private boolean checkUpGrid(int x, int y, int z, List<Coord> list,
			List<Coord> result) {
		if (y == 0)
			return false;
		if (steps[x][y - 1][z].getAlreadysearched() != 0)
			return false;
		double d1 = steps[x][y][z].getXdp() * steps[x][y][z].getXdp()
				+ steps[x][y][z].getYdp() * steps[x][y][z].getYdp()
				+ steps[x][y][z].getZdp() * steps[x][y][z].getZdp();
		double d2 = steps[x][y - 1][z].getXdp() * steps[x][y - 1][z].getXdp()
				+ steps[x][y - 1][z].getYdp() * steps[x][y - 1][z].getYdp()
				+ steps[x][y - 1][z].getZdp() * steps[x][y - 1][z].getZdp();
		if (d1 > d2)
			return false;

		Coord coord = new Coord(x, y - 1, z);
		if (!list.contains(coord) && !result.contains(coord)) {
			result.add(coord);
			steps[x][y - 1][z].setAlreadysearched(1);
		}
		return true;
	}

	/**
	 * 检查下方的网格，如果下方的grid的偏导值梯度大于当前网格的偏导值梯度，就返回true，否则返回false x 当前网格的X坐标 y
	 * 当前网格的Y坐标 z 当前网格的Z坐标 list 当前聚类所包含的网格列表(网格范围域) return 下方的网格是否符合条件
	 */
	private boolean checkDownGrid(int x, int y, int z, List<Coord> list,
			List<Coord> result) {
		if (y == K - 1)
			return false;
		if (steps[x][y + 1][z].getAlreadysearched() != 0)
			return false;
		double d1 = steps[x][y][z].getXdp() * steps[x][y][z].getXdp()
				+ steps[x][y][z].getYdp() * steps[x][y][z].getYdp()
				+ steps[x][y][z].getZdp() * steps[x][y][z].getZdp();
		double d2 = steps[x][y + 1][z].getXdp() * steps[x][y + 1][z].getXdp()
				+ steps[x][y + 1][z].getYdp() * steps[x][y + 1][z].getYdp()
				+ steps[x][y + 1][z].getZdp() * steps[x][y + 1][z].getZdp();
		if (d1 > d2)
			return false;

		Coord coord = new Coord(x, y + 1, z);
		if (!list.contains(coord) && !result.contains(coord)) {
			result.add(coord);
			steps[x][y + 1][z].setAlreadysearched(1);
		}
		return true;
	}

	/**
	 * 检查前方的网格，如果前方的grid的偏导值梯度大于当前网格的偏导值梯度，就返回true，否则返回false x 当前网格的X坐标 y
	 * 当前网格的Y坐标 z 当前网格的Z坐标 list 当前聚类所包含的网格列表(网格范围域) return 前方的网格是否符合条件
	 */
	private boolean checkFrontGrid(int x, int y, int z, List<Coord> list,
			List<Coord> result) {
		if (z == 0)
			return false;
		if (steps[x][y][z - 1].getAlreadysearched() != 0)
			return false;
		double d1 = steps[x][y][z].getXdp() * steps[x][y][z].getXdp()
				+ steps[x][y][z].getYdp() * steps[x][y][z].getYdp()
				+ steps[x][y][z].getZdp() * steps[x][y][z].getZdp();
		double d2 = steps[x][y][z - 1].getXdp() * steps[x][y][z - 1].getXdp()
				+ steps[x][y][z - 1].getYdp() * steps[x][y][z - 1].getYdp()
				+ steps[x][y][z - 1].getZdp() * steps[x][y][z - 1].getZdp();
		if (d1 > d2)
			return false;

		Coord coord = new Coord(x, y, z - 1);
		if (!list.contains(coord) && !result.contains(coord)) {
			result.add(coord);
			steps[x][y][z - 1].setAlreadysearched(1);
		}
		return true;
	}

	/**
	 * 检查后方的网格，如果后方的grid的偏导值梯度大于当前网格的偏导值梯度，就返回true，否则返回false x 当前网格的X坐标 y
	 * 当前网格的Y坐标 z 当前网格的Z坐标 list 当前聚类所包含的网格列表(网格范围域) return 后方的网格是否符合条件
	 */
	private boolean checkBackGrid(int x, int y, int z, List<Coord> list,
			List<Coord> result) {
		if (z == K - 1)
			return false;
		if (steps[x][y][z + 1].getAlreadysearched() != 0)
			return false;
		double d1 = steps[x][y][z].getXdp() * steps[x][y][z].getXdp()
				+ steps[x][y][z].getYdp() * steps[x][y][z].getYdp()
				+ steps[x][y][z].getZdp() * steps[x][y][z].getZdp();
		double d2 = steps[x][y][z + 1].getXdp() * steps[x][y][z + 1].getXdp()
				+ steps[x][y][z + 1].getYdp() * steps[x][y][z + 1].getYdp()
				+ steps[x][y][z + 1].getZdp() * steps[x][y][z + 1].getZdp();
		if (d1 > d2)
			return false;

		Coord coord = new Coord(x, y, z + 1);
		if (!list.contains(coord) && !result.contains(coord)) {
			result.add(coord);
			steps[x][y][z + 1].setAlreadysearched(1);
		}
		return true;
	}

	/**
	 * 依据网格的偏导势大小关系，搜索初始聚类中心 return 初始聚类中心坐标集
	 */
	ArrayList<Coord> markCandidatePoints() {
		// 计算各小网格的x偏导，y偏导，z偏导
		countPotential();
		// 沿各坐标轴找出初始聚类中心点集合
		ArrayList<Coord> xlist = getCanPointsOnX();
		ArrayList<Coord> ylist = getCanPointsOnY();
		ArrayList<Coord> zlist = getCanPointsOnZ();

		ArrayList<Coord> centers = new ArrayList<Coord>();
		// 求集合的交集得到初始聚类中心
		for (Coord c : xlist) {
			if (ylist.contains(c) && zlist.contains(c)) {
				if (steps[c.getX()][c.getY()][c.getZ()].pnumber > 0)
					centers.add(c);
			}
		}
		return centers;
	}

	/**
	 * 沿X轴方向搜索候选聚类中心 return X轴方向初始聚类中心坐标集
	 */
	private ArrayList<Coord> getCanPointsOnX() {
		ArrayList<Coord> points = new ArrayList<Coord>();
		// 比较空间相邻网格的X偏导势来选取候选聚类中心
		for (int y = 0; y < K; y++) {
			for (int z = 0; z < K; z++) {
				for (int x = 1; x < K; x++) {
					double dpx1 = steps[x - 1][y][z].getXdp();
					double dpx2 = steps[x][y][z].getXdp();
					if (dpx1 == 0 && dpx2 < 0)
						// dpx1是候选聚类中心
						points.add(new Coord(x - 1, y, z));
					else if (dpx1 > 0 && dpx2 == 0)
						// dpx2是候选聚类中心
						points.add(new Coord(x, y, z));
					// dpx1,dpx2均是候选聚类中心
					else if (dpx1 > 0 && dpx2 < 0) {
						points.add(new Coord(x, y, z));
						points.add(new Coord(x - 1, y, z));
					}
				}
			}
		}
		return points;
	}

	/**
	 * 沿Y轴方向搜索候选聚类中心 return X轴方向初始聚类中心坐标集
	 */
	private ArrayList<Coord> getCanPointsOnY() {
		ArrayList<Coord> points = new ArrayList<Coord>();
		// 比较空间相邻网格的Y偏导势来选取候选聚类中心
		for (int z = 0; z < K; z++) {
			for (int x = 0; x < K; x++) {
				for (int y = 1; y < K; y++) {
					double dpy1 = steps[x][y - 1][z].getYdp();
					double dpy2 = steps[x][y][z].getYdp();
					if (dpy1 == 0 && dpy2 < 0)
						// dpy1是候选聚类中心
						points.add(new Coord(x, y - 1, z));
					else if (dpy1 > 0 && dpy2 == 0)
						// dpy2是候选聚类中心
						points.add(new Coord(x, y, z));
					else if (dpy1 > 0 && dpy2 < 0) {
						points.add(new Coord(x, y, z));
						points.add(new Coord(x, y - 1, z));
					}
				}
			}
		}
		return points;
	}

	/**
	 * 沿Z轴方向搜索候选聚类中心 return X轴方向初始聚类中心坐标集
	 */
	private ArrayList<Coord> getCanPointsOnZ() {
		ArrayList<Coord> points = new ArrayList<Coord>();
		// 比较空间相邻网格的Z偏导势来选取候选聚类中心
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 1; z < K; z++) {
					double dpy1 = steps[x][y][z - 1].getZdp();
					double dpy2 = steps[x][y][z].getZdp();
					if (dpy1 == 0 && dpy2 < 0)
						// dpy1是候选聚类中心
						points.add(new Coord(x, y, z - 1));
					else if (dpy1 > 0 && dpy2 == 0)
						// dpy2是候选聚类中心
						points.add(new Coord(x, y, z));
					else if (dpy1 > 0 && dpy2 < 0) {
						points.add(new Coord(x, y, z));
						points.add(new Coord(x, y, z - 1));
					}
				}
			}
		}
		return points;
	}

	/**
	 * 计算每个小网格的x偏导，y偏导，z偏导
	 */
	private void countPotential() {
		// 根据小网格初始化大网格
		int BK = K / IFP;
		grids = new BigGrid[BK][BK][BK];
		for (int kz = 0; kz < BK; kz++) {
			for (int ky = 0; ky < BK; ky++) {
				for (int kx = 0; kx < BK; kx++) {
					double xall = 0;
					double yall = 0;
					double zall = 0;
					double spacialXall = 0;
					double spacialYall = 0;
					int pall = 0;
					// 小网格属性的累加
					for (int z = 0; z < IFP; z++) {
						for (int y = 0; y < IFP; y++) {
							for (int x = 0; x < IFP; x++) {
								int xnum = kx * IFP + x;
								int ynum = ky * IFP + y;
								int znum = kz * IFP + z;

								xall += steps[xnum][ynum][znum].getXall();
								yall += steps[xnum][ynum][znum].getYall();
								zall += steps[xnum][ynum][znum].getZall();
								spacialXall += steps[xnum][ynum][znum]
										.getSpacialXall();
								spacialYall += steps[xnum][ynum][znum]
										.getSpacialYall();
								pall += steps[xnum][ynum][znum].getWeight();
							}
						}
					}
					if (pall > 0) {
						BigGrid grid = new BigGrid(xall / pall, yall / pall,
								zall / pall, spacialXall / pall, spacialYall
										/ pall, pall);
						grids[kx][ky][kz] = grid;
						notEmptyBigGrid.add(grids[kx][ky][kz]);
						allBigGrid.add(grids[kx][ky][kz]);

					} else {
						// 空大网格的特征中心设为几何中心
						BigGrid grid = new BigGrid(xmin + kx * IFP * dx + IFP
								* dx / 2, ymin + ky * IFP * dy + IFP * dy / 2,
								zmin + kz * IFP * dz + IFP * dz / 2, 0, 0, pall);
						grids[kx][ky][kz] = grid;
						allBigGrid.add(grids[kx][ky][kz]);
					}

				}
			}
		}

		// 非空小网格的坐标集
		ArrayList<Coord> notEmptyGridCoord = new ArrayList<Coord>();
		for (int y = 0; y < K; y++) {
			for (int z = 0; z < K; z++) {
				for (int x = 0; x < K; x++) {
					if (steps[x][y][z].pnumber >= 1)
						notEmptyGridCoord.add(new Coord(x, y, z));
				}
			}
		}
		// 非空大网格的坐标集
		ArrayList<Coord> notEmptyBigGridCoord = new ArrayList<Coord>();
		for (int y = 0; y < BK; y++) {
			for (int z = 0; z < BK; z++) {
				for (int x = 0; x < BK; x++) {
					if (grids[x][y][z].pnumber >= 1)
						notEmptyBigGridCoord.add(new Coord(x, y, z));
				}
			}
		}

		// 以下输出大,小网格的信息至文本文件
		String s = new String();
		try {
			File fo = new File("C:\\Users\\YinJF\\Work\\CenterOfGravity.txt");
			File fb = new File(
					"C:\\Users\\YinJF\\Work\\CenterOfGravityBigGrid.txt");
			if (fo.exists()) {
				System.out.println("文件og存在,删除以前的...");
				fo.delete();
			}
			if (fb.exists()) {
				System.out.println("文件bg存在,删除以前的...");
				fb.delete();
			}
			if (fo.createNewFile()) {
				System.out.println("文件og创建成功！");
			} else {
				System.out.println("文件og创建失败！");
			}
			if (fb.createNewFile()) {
				System.out.println("文件bg创建成功！");
			} else {
				System.out.println("文件bg创建失败！");
			}

			BufferedWriter outputfb = new BufferedWriter(new FileWriter(fb));
			// 大网格重心输出
			for (BigGrid bigGrid : allBigGrid) {
				System.out.println("------所有大网格: ["
						+ (int) Math.floor((bigGrid.getXloc() - xmin)
								/ (dx * IFP))
						+ "]["
						+ (int) Math.floor((bigGrid.getYloc() - ymin)
								/ (dy * IFP))
						+ "]["
						+ (int) Math.floor((bigGrid.getZloc() - zmin)
								/ (dz * IFP)) + "]" + " 有点"
						+ bigGrid.getWeight() + " 个------");
				s = Math.round(bigGrid.getXloc() * 100) / 100.0 + " "
						+ Math.round(bigGrid.getYloc() * 100) / 100.0 + " "
						+ Math.round(bigGrid.getZloc() * 100) / 100.0 + " "
						+ +bigGrid.pnumber + "\r\n";
				// 网格重心输出
				outputfb.write(s);
			}
			outputfb.close();
			// 小网格重心及所含数据点数输出
			BufferedWriter outputfo = new BufferedWriter(new FileWriter(fo));
			for (Grid steps : notEmptyGrid) {
				System.out.println("------非空小网格 : ["
						+ (int) Math.floor((steps.getXloc() - xmin) / dx)
						+ "]["
						+ (int) Math.floor((steps.getYloc() - ymin) / dy)
						+ "]["
						+ (int) Math.floor((steps.getZloc() - zmin) / dz) + "]"
						+ " 有点" + steps.getWeight() + " 个------");
				s = Math.round(steps.getXloc() * 100) / 100.0 + " "
						+ Math.round(steps.getYloc() * 100) / 100.0 + " "
						+ Math.round(steps.getZloc() * 100) / 100.0 + " "
						+ steps.getWeight() + "\r\n";
				// 网格重心输出
				outputfo.write(s);
			}
			outputfo.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		// 以下输出小网格权重至文本文件
		String sw = new String();
		try {
			File fo = new File("C:\\Users\\YinJF\\Work\\GridWeight.txt");

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
						double xloc = steps[x][y][z].getXloc();
						double yloc = steps[x][y][z].getYloc();
						double zloc = steps[x][y][z].getZloc();

						double spacialX = steps[x][y][z].getSpacialX();
						double spacialY = steps[x][y][z].getSpacialY();

						double dpx = 0;
						double dpy = 0;
						double dpz = 0;
						double potential = 0;
						double gridWeight = 0;

						// 按数据场势函数的定义，计算各小网格的势值和偏导势
						for (int bx = 0; bx < BK; bx++) {
							for (int by = 0; by < BK; by++) {
								for (int bz = 0; bz < BK; bz++) {
									BigGrid bg = grids[bx][by][bz];
									double dx = (bg.getXloc() - xloc) / sig1;
									double dy = (bg.getYloc() - yloc) / sig2;
									double dz = (bg.getZloc() - zloc) / sig3;

									double spacialDx = (bg.getSpacialX() - spacialX)
											/ sigSpacial1;
									double spacialDy = (bg.getSpacialY() - spacialY)
											/ sigSpacial2;
									// 特征距离
									double dist2 = dx * dx + dy * dy + dz * dz;
									// 空间距离
									double spacialdist2 = spacialDx * spacialDx
											+ spacialDy * spacialDy;

									double temp = bg.getWeight()
											* Math.pow(Math.E,
													(double) (-dist2))
											* Math.pow(Math.E,
													(double) (-spacialdist2));

									gridWeight += bg.getWeight()
											* Math.pow(Math.E,
													(double) (-spacialdist2));

									potential += temp;
									dpx += dx * temp;
									dpy += dy * temp;
									dpz += dz * temp;
								}
							}
						}

						sw = gridWeight + "\r\n";
						outputfo.write(sw);
						steps[x][y][z].setPotential(potential);
						steps[x][y][z].setXdp(dpx);
						steps[x][y][z].setYdp(dpy);
						steps[x][y][z].setZdp(dpz);
					}
				}
			}

			outputfo.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * 类Coord定义空间网格的坐标对象
	 */
	public class Coord {
		// X，Y，Z轴坐标值
		private int x;
		private int y;
		private int z;

		// 单个网格所属的类号
		private int clusterNum;

		public int getClusterNum() {
			return clusterNum;
		}

		public void setClusterNum(int clusterNum) {
			this.clusterNum = clusterNum;
		}

		public Coord(int x, int y, int z) {
			this.setX(x);
			this.setY(y);
			this.setZ(z);
		}

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
			Coord other = (Coord) obj;
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
			info = "Coord [" + this.getX() + "]" + "[" + this.getY() + "]"
					+ "[" + this.getZ() + "]";
			return info;
		}

		private LUVtest getOuterType() {
			return LUVtest.this;
		}

		public void setX(int x) {
			this.x = x;
		}

		public int getX() {
			return x;
		}

		public void setY(int y) {
			this.y = y;
		}

		public int getY() {
			return y;
		}

		public void setZ(int z) {
			this.z = z;
		}

		public int getZ() {
			return z;
		}
	}

	/**
	 * 类BigGrid定义空间大网格对象
	 */
	private class BigGrid {
		// 大网格的X，Y，Z轴坐标值
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		private double spacialX = -1;

		private double spacialY = -1;

		// 大网格内包含的数据点量
		private int pnumber = 0;

		public BigGrid(double xloc, double yloc, double zloc, double spacialX,
				double spacialY, int pnumber) {
			this.xloc = xloc;
			this.yloc = yloc;
			this.zloc = zloc;

			this.spacialX = spacialX;
			this.spacialY = spacialY;
			this.pnumber = pnumber;
		}

		public double getSpacialX() {
			return spacialX;
		}

		public double getSpacialY() {
			return spacialY;
		}

		public double getXloc() {
			return xloc;
		}

		public double getYloc() {
			return yloc;
		}

		public double getZloc() {
			return zloc;
		}

		public int getWeight() {
			return pnumber;
		}
	}

	/**
	 * 类Cluster定义网格聚类对象
	 */
	public class Cluster {
		// 聚类编号
		int clusterLabel;

		// 三轴坐标
		double xloc = -1;

		double yloc = -1;

		double zloc = -1;

		// 各坐标轴的数据值统计
		double xAll = 0;

		double yAll = 0;

		double zAll = 0;

		double distanceSQ = 0;

		// 聚类中包含的数据点总数
		double numAll = 0;

		// 聚类范围域(上下文)网格坐标集
		LinkedList<Coord> contentList;

		// 聚类中心网格
		Coord clusterCenter;

		// 初始聚类中心网格
		Coord initialClusterCenter;

		// 找聚类初始中心的六邻域
		public Cluster(int x, int y, int z, int clusterNumber) {
			initialClusterCenter = new Coord(x, y, z);
			this.clusterLabel = clusterNumber;
			contentList = new LinkedList<Coord>();
			contentList.add(initialClusterCenter);
			if (x > 0)
				contentList.add(new Coord(x - 1, y, z));
			if (x < K - 1)
				contentList.add(new Coord(x + 1, y, z));
			if (y > 0)
				contentList.add(new Coord(x, y - 1, z));
			if (y < K - 1)
				contentList.add(new Coord(x, y + 1, z));
			if (z > 0)
				contentList.add(new Coord(x, y, z - 1));
			if (z < K - 1)
				contentList.add(new Coord(x, y, z + 1));
		}

		public void setContentList(LinkedList<Coord> list) {
			contentList = list;
		}

		public Coord getClusterCenter() {
			return clusterCenter;
		}

		public void setClusterCenter(Coord clusterCenter) {
			this.clusterCenter = clusterCenter;
		}

		// 合并类后计算各自重心
		public void adjustClusterCenter() {
			this.setNumAll();

			int x = (int) Math.floor((this.getxAll() / this.getNumAll() - xmin)
					/ dx);
			if (x == K)
				x = K - 1;
			int y = (int) Math.floor((this.getyAll() / this.getNumAll() - ymin)
					/ dy);
			if (y == K)
				y = K - 1;
			int z = (int) Math.floor((this.getzAll() / this.getNumAll() - zmin)
					/ dz);
			if (z == K)
				z = K - 1;

			this.setClusterCenter(new Coord(x, y, z));
		}

		public double getxAll() {
			xAll = 0;
			for (Coord c : contentList) {
				xAll += steps[c.getX()][c.getY()][c.getZ()].getXall();
			}
			return xAll;
		}

		public double getyAll() {
			yAll = 0;
			for (Coord c : contentList) {
				yAll += steps[c.getX()][c.getY()][c.getZ()].getYall();
			}
			return yAll;
		}

		public double getzAll() {
			zAll = 0;
			for (Coord c : contentList) {
				zAll += steps[c.getX()][c.getY()][c.getZ()].getZall();
			}
			return zAll;
		}

		public void setNumAll() {
			for (Coord c : contentList)
				numAll += steps[c.getX()][c.getY()][c.getZ()].getWeight();
		}

		public double getNumAll() {
			return numAll;
		}

		public double computeDistance(Cluster c) {

			distanceSQ = (this.getxAll() / this.getNumAll() - c.getxAll()
					/ c.getNumAll())
					* (this.getxAll() / this.getNumAll() - c.getxAll()
							/ c.getNumAll())
					+ (this.getyAll() / this.getNumAll() - c.getyAll()
							/ c.getNumAll())
					* (this.getyAll() / this.getNumAll() - c.getyAll()
							/ c.getNumAll())
					+ (this.getzAll() / this.getNumAll() - c.getzAll()
							/ c.getNumAll())
					* (this.getzAll() / this.getNumAll() - c.getzAll()
							/ c.getNumAll());
			return distanceSQ;
		}

		public Cluster(int clusterLabel) {
			this.clusterLabel = clusterLabel;
			contentList = new LinkedList<Coord>();
		}

		void addCoord(Coord c) {
			if (xloc == -1 && yloc == -1 && zloc == -1) {
				xloc = steps[c.getX()][c.getY()][c.getZ()].getXloc();
				yloc = steps[c.getX()][c.getY()][c.getZ()].getYloc();
				zloc = steps[c.getX()][c.getY()][c.getZ()].getZloc();
			}
			contentList.add(c);
		}

		public List<Coord> getCoordList() {
			return contentList;
		}

		public void addCoords(List<Coord> coords) {
			for (Coord c : coords) {
				contentList.add(c);
			}
		}

		public Coord getInitialClusterCenter() {
			return initialClusterCenter;
		}
	}

	/**
	 * 类Segment定义图像分块对象
	 */
	private class Segment {
		// 块编号
		int segmentLabel;

		// 组编号
		int groupLabel;

		// 块中心特征空间三轴坐标
		double xloc = -1;

		double yloc = -1;

		double zloc = -1;

		// 特征空间各坐标值统计
		double xAll = 0;

		double yAll = 0;

		double zAll = 0;;

		// 块中包含的数据点总数
		double numAll = 0;

		// 块内数据点集合
		LinkedList<DataPoint> contentList;

		// 块内边界点集合
		LinkedList<DataPoint> marginList;

		// 相邻块的编号集合
		LinkedList<Integer> neighborSegLabelList;

		public void setSegmentLabel(int segmentLabel) {
			this.segmentLabel = segmentLabel;
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

		public Segment(int segmentLabel, int groupLabel) {
			this.segmentLabel = segmentLabel;
			this.groupLabel = groupLabel;
			contentList = new LinkedList<DataPoint>();
			marginList = new LinkedList<DataPoint>();
			neighborSegLabelList = new LinkedList<Integer>();
		}

	}

	/**
	 * 类Group定义图像猪对象,及块的的集合
	 */
	private class Group {
		// 组编号
		int groupLabel;

		// 住中心特征空间三轴坐标
		double xloc = -1;

		double yloc = -1;

		double zloc = -1;

		// 特征空间各坐标值统计
		double xAll = 0;

		double yAll = 0;

		double zAll = 0;

		// 块中包含的数据点总数
		double numAll = 0;

		// 块内数据点集合
		LinkedList<DataPoint> contentList;

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

		public Group(int groupLabel) {
			this.groupLabel = groupLabel;
			contentList = new LinkedList<DataPoint>();
		}

	}

	/**
	 * 类Grid定义空间小网格对象
	 */
	public class Grid {
		// 各坐标轴的数据值统计
		private double xall = 0;

		private double yall = 0;

		private double zall = 0;

		private double spacialXall = 0;

		private double spacialYall = 0;

		// 小网格的X，Y，Z轴坐标值,及特征坐标平均值
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// 小网格空间坐标平均值
		private double spacialX = -1;

		private double spacialY = -1;

		// 三个偏导势
		private double xdp = -1;

		private double ydp = -1;

		private double zdp = -1;

		// 小网格包含的数据点总数
		private int pnumber = 0;

		// 小网格包含的所有数据点集
		private ArrayList<DataPoint> list = new ArrayList<DataPoint>();

		// 小网格是否被搜索过
		private int alreadysearched = 0;

		// 小网格的势值
		private double potential = -1;

		public double getPotential() {
			return potential;
		}

		public void setPotential(double potential) {
			this.potential = potential;
		}

		public void add(DataPoint point) {
			pnumber++;
			xall += point.getCoord_X();
			yall += point.getCoord_Y();
			zall += point.getCoord_Z();

			spacialXall += point.getLine();
			spacialYall += point.getRow();

			list.add(point);
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

		public int getAlreadysearched() {
			return alreadysearched;
		}

		public void setAlreadysearched(int alreadysearched) {
			this.alreadysearched = alreadysearched;
		}

		public void setCenter(double xloc, double yloc, double zloc) {
			this.xloc = xloc;
			this.yloc = yloc;
			this.zloc = zloc;
		}

		public double getXall() {
			return xall;
		}

		public double getYall() {
			return yall;
		}

		public double getZall() {
			return zall;
		}

		public double getXloc() {
			if (xloc == -1) {
				if (pnumber != 0)
					xloc = xall / pnumber;
			}
			return xloc;
		}

		public double getYloc() {
			if (yloc == -1) {
				if (pnumber != 0)
					yloc = yall / pnumber;
			}
			return yloc;
		}

		public double getZloc() {
			if (zloc == -1) {
				if (pnumber != 0)
					zloc = zall / pnumber;
			}
			return zloc;
		}

		public double getSpacialX() {
			if (spacialX == -1) {
				if (pnumber != 0)
					spacialX = spacialXall / pnumber;
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
				if (pnumber != 0)
					spacialY = spacialYall / pnumber;
				else
					spacialY = 0;
			}
			return spacialY;
		}

		public void setSpacialY(double spacilaY) {
			this.spacialY = spacilaY;
		}

		public int getWeight() {
			return pnumber;
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

		public ArrayList<DataPoint> getList() {
			return list;
		}

		public void setList(ArrayList<DataPoint> list) {
			this.list = list;
		}
	}
}
