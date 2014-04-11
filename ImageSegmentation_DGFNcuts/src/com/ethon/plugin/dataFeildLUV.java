package com.ethon.plugin;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
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
public class dataFeildLUV {

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
	public dataFeildLUV(int K, int IFP, int line_num, int row_num, double sig1,
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

			LUVPoints[i].setCoord_X(R);
			LUVPoints[i].setCoord_Y(G);
			LUVPoints[i].setCoord_Z(B);

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
		toWrite = DataWriter.transformToImagePoints2(toWrite, clusters, steps,
				allCoords, K);

		// 特征坐标回归为(R,G,B)
		toWrite = LUVToRGB(toWrite);

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
		g.drawString("K=" + K + ",IFP=" + IFP, 0, DataBase.sHei - 40);
		g.drawString("共" + toWrite.length + "个点，聚成" + clusters.size() + "个类",
				0, DataBase.sHei - 25);

		g.drawString(time, 0, DataBase.sHei - 10);
		g.dispose();

		g.dispose();
		DataField.updateImagePanel(img);

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
										* Math.pow(Math.E, (double) (-dist2))
										* Math.pow(Math.E,
												(double) (-spacialdist2));
								potential += temp;
								dpx += dx * temp;
								dpy += dy * temp;
								dpz += dz * temp;
							}
						}
					}
					steps[x][y][z].setPotential(potential);
					steps[x][y][z].setXdp(dpx);
					steps[x][y][z].setYdp(dpy);
					steps[x][y][z].setZdp(dpz);
				}
			}
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

		private dataFeildLUV getOuterType() {
			return dataFeildLUV.this;
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
