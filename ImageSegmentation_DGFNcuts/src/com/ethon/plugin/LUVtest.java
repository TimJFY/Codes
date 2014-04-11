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
 * LUVtest����ݶ��㷨���������ݳ�����,�������,���������ͼƬ�ָ�����άͼ������DataPoint[L,U,V]���о������,ƽ������,���ݷֿ�,
 * ���ں�
 * 
 * @author Jinfei Yin
 * 
 */
public class LUVtest {

	// ͼƬ���п���п�
	int line_num, row_num;

	// ÿ��������ı���������С���������
	private int IFP;

	// С����ı߳������
	public double dx, dy, dz, xmin, ymin, zmin;

	// ���񻮷ֵ������ռ�Ӱ�����Ӻ�����ռ�Ӱ�����ӣ��Ƽ�ֵͨ��matlab����
	private double sig1, sig2, sig3, sigSpacial1, sigSpacial2;

	// С�������͵ľ���
	public Grid[][][] steps;

	// С�����������ľ���
	public Coord[][][] allCoords;

	// ����Ŀ�ʼ����ʱ��
	private long start;

	// ������������з�����ʱ��
	private long end;

	// ����ʱ���������Ĳ���
	// private int thread;

	// ����ϵͳÿ�߰���С���������
	private int K;

	// ���������͵ľ���
	private BigGrid[][][] grids;

	// ����ԭʼͼ�����ݵ�
	DataPoint[] rawPoints;

	// ����LUV���ݵ�
	DataPoint[] points;

	// ���վ�����Ŀ
	private int finalClusterNum;

	// ԭʼͼ�����ݵ������
	private int allPointsNum = 0;

	// ������
	// ����ϵͳ�����зǿ�С����
	private ArrayList<Grid> notEmptyGrid = new ArrayList<Grid>();

	// ����ϵͳ�����ǿմ�����
	private ArrayList<BigGrid> notEmptyBigGrid = new ArrayList<BigGrid>();

	// ����ϵͳ�����д�����
	private ArrayList<BigGrid> allBigGrid = new ArrayList<BigGrid>();

	// �þ����У��������ݵ㰴����ԭͼ���е�����(x,y)��¼�ڶ�Ӧ��λ�ã�ԭͼ��Ϊline_num*row_num���˴�����ÿɱ��������
	DataPoint[][] allPix;

	// ���ƽ����ͼ�������еĿ�,Ҳ��ͼ�񰴿黮�ֵ����ս��
	LinkedList<Segment> segments = new LinkedList<Segment>();

	// ����������ı��
	int segmentGroupLabelIndex = 0;

	// ƽ����ֵ
	double smoothValve;

	// �����ָ���ֵ
	double nCutValve;

	// ��ָ�Ӱ������
	double matrixSig;

	// ÿ�����ɵĸ����ڵ����
	int auxiliaryNum;

	// ȫͼƬ�����W
	Matrix WPrior;

	// ��¼�������ݵĵ�
	private DataPoint[] points1 = null;

	/**
	 * ���캯��HastaAndSNC�������û��������K, IFP��������ԭʼͼ�����ݵ�ӳ��������ϵͳ����ȷ����ʼ��ȫ��С�������
	 * Parameters:K:����ϵͳÿ�߰���С��������� IFP��ÿ��������ı���������С��������� line_num:ͼƬ���п�
	 * row_num:ͼƬ���п�
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
		// ʵ������������
		rawPoints = DataBase.getInstance().getPoints();
		points1 = rawPoints;
		DataBase db = new DataBase();
		// ����ȫ�ֱ���points1
		points = db.new_public_points(points1);

		// GRBֵת��ΪLUV
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
			// ������ά����ϵͳ��ά�����
			if (xmin > x)
				xmin = x;
			if (ymin > y)
				ymin = y;
			if (zmin > z)
				zmin = z;
			// ������ά����ϵͳ��ά���յ�
			if (xmax < x)
				xmax = x;
			if (ymax < y)
				ymax = y;
			if (zmax < z)
				zmax = z;
		}
		// ������ά����ϵͳ��ά�Ĳ���(����)
		dx = (xmax - xmin) / K;
		dy = (ymax - ymin) / K;
		dz = (zmax - zmin) / K;

		// Ӱ�����ӵĻ���ȷ������
		// double sigmaPre = Math.max(dx, dy);
		// sigma = Math.max(sigmaPre, dz);
		// sigma = sigma * IFP;

		// ��ÿ�����ݵ���䵽��Ӧ��С�����У�������һ����С����Ϊ��
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

			// ��¼���зǿ�����
			if (!notEmptyGrid.contains(steps[xl][yl][zl]))
				notEmptyGrid.add(steps[xl][yl][zl]);
		}

		// ȷ������С���������ʼ��
		// ͳ��С�����������
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					double xloc = xmin + dx / 2 + dx * x;
					double yloc = ymin + dy / 2 + dy * y;
					double zloc = zmin + dz / 2 + dz * z;

					Grid grid = steps[x][y][z];
					if (grid == null) {
						grid = new Grid();
						// С�������������Ϊ�伸������
						grid.setCenter(xloc, yloc, zloc);
						steps[x][y][z] = grid;
					}
					// ����С������������
					allCoords[x][y][z] = new Coord(x, y, z);
				}
			}
		}
	}

	/**
	 * ����RGBToLUV����������RGBPoints,������ԭʼͼ�����ݵ��������ά����ֵ(R,G,B)����ʽת��Ϊ(L,U,V),���ھ���ͷָ����
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

		// ����������ݵ�(L,U,V)ֵ���ı��ļ�
		String s = new String();
		try {
			File fo = new File("C:\\Users\\YinJF\\Work\\DataPointsLUV.txt");

			if (fo.exists()) {
				System.out.println("�ļ�DataPointsLUV����,ɾ����ǰ��...");
				fo.delete();
			}
			if (fo.createNewFile()) {
				System.out.println("�ļ�DataPointsLUV�����ɹ���");
			} else {
				System.out.println("�ļ�DataPointsLUV����ʧ�ܣ�");
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
	 * ����LUVToRGB����������LUVPoints,���������õ�ͼ�����ݵ��������ά����ֵ(L,U,V)ת��Ϊ(R,G,B),����ͼ�����
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
				System.out.println("�ļ�DataPointsLUVRGBC����,ɾ����ǰ��...");
				fo.delete();
			}
			if (fo.createNewFile()) {
				System.out.println("�ļ�DataPointsLUVRGBC�����ɹ���");
			} else {
				System.out.println("�ļ�DataPointsLUVRGBC����ʧ�ܣ�");
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
	 * �����߼����������ֿ鴦����process���㷨����
	 * 
	 * 1.���ú���markCandidatePoints()������ͼ�����ݵ����������ݳ��е�ƫ���ƴ�С��ϵѰ����ֵ�ļ���㣬ȷ�����г�ʼ������������
	 * 2.�Գ�ʼ������������������򣬵��ú���markCandidatePoints�����ռ�������ԭ��ϲ����ڵĳ�ʼ�������ģ��ϲ�������ĵĳ�ʼ������(
	 * ��Χ��)���غϽ���ȫ���� 3.���ú���searchforCluster(Cluster
	 * cluster)�����ݺϲ���ľ������ģ���������������������Χ�򣬲������������÷�Χ���ڵ��������ݵ������
	 * 4.���ú���transformToImagePoints(toWrite, clusters, steps, allCoords,
	 * K),�һؾ�������е�������ʧ�㲢���߽�����㰴������֪����ŷʽ������̵�ԭ�����·�������
	 * 5.���ú���smoothpix(DataPoint[] points, List<Cluster>
	 * clusters)���Է���������ƽ����ȥ�룬�޸��������ݵ�Ķ�Ӧ���� 6.���ú���initSegments(toWrite,
	 * segments)�������ݵ������,��ʼ������ͼƬ�ֿ�,ȷ��ÿ�����ʼ���,�����ĵ㼯�Ϻ����ڿ�ŵļ���
	 * 7.���ú���generateAuxiliarySeg(segments,
	 * auxiliaryNum)Ϊÿ��������ָ���ĵĶ��⸨����,���Ķ���Ӧ�Ŀ�ź��ڿ���Ϣ
	 * 8.���ú���matrixWPriorByAverageValue
	 * /matrixWPriorByCovariance/matrixWPriorBySegAverageValue(segments,
	 * auxiliaryNum)���������ֿ���������W(��ʼ����)
	 * 9.���ú���matrixDeclareAndSovle(segments)����ͼƬ���ʼ�����ʼ����ǰ��W��D
	 * ,������ֵ�����������õ������ӿ�,�����д��ȫ�ֱ���,�����еݹ�����,����ɷֿ����
	 * 10.���ú���segmentsArrange(segments,
	 * toWrite)����ͳ�Ʒֿ���,Ϊͼ�����ݵ��趨��Ų���������ĵ��������ά����
	 * 11.���ú���LUVToRGB(toWrite)��������ع�Ϊ(R,G,B)
	 * 12.��ʼ����ͼAPI�����ú���show_DataPoint_on_Image(BufferedImage image,DataPoint
	 * point,Color color)��ͼ��ʾ������
	 */
	public void process() {
		// ���г�ʼ�����������������
		ArrayList<Coord> centers = markCandidatePoints();
		// ���г�ʼ��������
		ArrayList<Cluster> initialClusters = new ArrayList<Cluster>();
		// ��¼����ʼ�������Ƿ��Ѻϲ�����-1��ʾ��δ�ϲ��������ʾ�ѱ��ϲ�
		ArrayList<Integer> initMerged = new ArrayList<Integer>();
		// �ϲ���ĳ�ʼ��������
		LinkedList<Cluster> finnalInitList = new LinkedList<Cluster>();
		// ����������ľ�������
		ArrayList<Cluster> clusters = new ArrayList<Cluster>();

		// �洢����������ͼ�����ݵ�
		DataPoint[] toWrite;

		// ��������ֵ��ʼ�����г�ʼ��������
		for (int x = 0, clusterNumber = 1; x < centers.size(); x++) {
			Cluster cluster = new Cluster(centers.get(x).getX(), centers.get(x)
					.getY(), centers.get(x).getZ(), clusterNumber++);
			initialClusters.add(cluster);
		}

		// ���������İ�����ά����ֵ��С��������
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
		// ����������ϲ���ʼ��������, ������������˳��������
		for (int i = 0; i < initialClusters.size(); i++) {
			// ÿ�κϲ���ɺ����õ����¾���ķ�Χ��
			LinkedList<Coord> finnalInitCoordList = new LinkedList<Coord>();
			if (!initMerged.contains(i)) {
				int j = 0;
				// �ϲ���ľ�������
				Cluster mergedCluster = initialClusters.get(i);
				finnalInitCoordList = initialClusters.get(i).contentList;
				while (j < initialClusters.size()) {
					Coord toVerify = initialClusters.get(j).initialClusterCenter;
					// ������ж��ľ������������Ѿ�������һ��ķ�Χ������������������ķ�Χ����ȥ�صĻ����Ͻ��кϲ�
					if (!initMerged.contains(j) && i != j
							&& finnalInitCoordList.contains(toVerify)) {
						for (Coord coord : initialClusters.get(j).contentList) {
							if (!finnalInitCoordList.contains(coord)) {
								finnalInitCoordList.add(coord);
							}
						}
						// ��ǳɹ��ϲ��ĳ�ʼ�������ģ������ظ��ϲ�
						initMerged.set(j, j);
						// ��ĳһ�ಢ�����µľ������ģ�����������(��Χ��)�����󣬽���ͷ��ʼ����ԭʼ���������б�
						j = 0;
					} else {
						// ���������þ�����������ֱ��������Χ���⣬�򻯲���
						if (i == j) {
							initMerged.set(j, j);
						}
						j++;
					}
				}
				int num = 0;
				// ȷ��ÿ��������౻����һ����������һ�࣬getAlreadysearched()����ֵ1��ʾ�ѱ�������������㽫�����벽��3
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
			// �ϲ���ĳ�ʼ��������
			Cluster inicluser = new Cluster(-1, -1, -1, -1);
			inicluser = finnalInitList.get(i);

			// ���ϲ���������ĵĳ�ʼ����Χ������������������Χ��
			inicluser = searchforCluster(inicluser);
			// �������õ�������Χ�������������������񣬼���������
			inicluser.adjustClusterCenter();
			clusters.add(inicluser);

		}
		// �����û�������ܶ���ֵ��ȥ��������
		// String input = JOptionPane.showInputDialog("�����������ܶȷ�ֵ",
		// Integer.toString(thread - 1 > 0 ? thread - 1 : 0));
		// thread = Integer.parseInt(input.trim());

		finalClusterNum = clusters.size();
		toWrite = new DataPoint[allPointsNum];

		// ��ʼ��
		for (int i = 0; i < allPointsNum; i++) {
			toWrite[i] = points[i];
		}
		// ͳ�ƾ����㷨��������ݵ������
		int toWriteS = 0;
		// ͳ�ƾ�����������д�������ݵ������
		int j = 0;

		for (int i = 0; i < clusters.size(); i++) {
			Cluster cluster = clusters.get(i);
			List<Coord> list = cluster.getCoordList();
			// Ϊ��������(����)�������
			for (int coordNum = 0; coordNum < list.size(); coordNum++) {
				allCoords[list.get(coordNum).getX()][list.get(coordNum).getY()][list
						.get(coordNum).getZ()].setClusterNum(i + 1);
			}

			int clustersize = 0;
			// ���������������о����㷨���������ݵ���������(L,U,V)������ClusterLabel
			for (Coord c : list) {
				Grid grid = steps[c.getX()][c.getY()][c.getZ()];
				ArrayList<DataPoint> plist = grid.list;
				clustersize += plist.size();
				// ���ĵ�����
				for (DataPoint p : plist) {
					p.setClusterLabel(i);
					toWrite[(p.getLineNum() - 1)] = p;
					toWriteS++;
				}
			}
			j += clustersize;
			System.out.println(" ��  " + (int) (i + 1) + " ��" + "���� "
					+ clustersize + " ����");
		}

		System.out.println("------���� "
				+ j
				+ " �����¼,"
				+ " ƽ��ǰ��ʧ�� "
				+ (points.length - j)
				+ " ��,"
				+ " ռ���� "
				+ new DecimalFormat("0.0000")
						.format((double) (points.length - j) / points.length)
				+ " ------");

		System.out.println(" -----����  " + toWriteS + " ����д��------");

		// ȥ�뺯��:1.��ʧ���һ� 2.�߽���������·��� 3.�����ݵ��ʼ�����ͼ������ص�
		toWrite = DataWriter.transformToImagePoints(toWrite, clusters, steps,
				allCoords, K);

		// �Ծ��������ݵ����ƽ��
		toWrite = smoothpix(toWrite, clusters);

		// ��ʼ������ͼƬ�ֿ�
		segments = initSegments(toWrite, segments);

		// ���ɸ����ڵ�
		generateAuxiliarySeg(segments, auxiliaryNum);

		// ��������
		matrixWPriorByAverageValue(segments, auxiliaryNum);
		// matrixWPriorByCovariance(segments, auxiliaryNum);
		// matrixWPriorBySegAverageValue(segments, auxiliaryNum);

		// ����ͼƬ�����������ֵ�����������õ������ӿ�,�����е�������,����ɷֿ����
		matrixDeclareAndSovle(segments);

		// ����ֿ���Ϣ,Ϊͼ�����ݵ��趨���
		toWrite = segmentsArrange(segments, toWrite);

		// ��������ع�Ϊ(R,G,B)
		toWrite = LUVToRGB(toWrite);

		// �����ͼAPI, ׼����ͼ
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
		g.drawString("��" + toWrite.length + "���㣬�۳�" + clusters.size() + "����",
				0, DataBase.sHei - 25);

		g.drawString(time, 0, DataBase.sHei - 10);
		g.dispose();

		g.dispose();
		DataField.updateImagePanel(img);

	}

	/**
	 * ����smoothpix������ݾ����ƽ��ȥ�� Parameters:points: �������Ѻ��ƽ�������ݼ� clusters:
	 * ������ľ��������б� Return ƽ��������ݼ���
	 * 
	 */
	public DataPoint[] smoothpix(DataPoint[] points, List<Cluster> clusters) {

		// ��ʼ��allPix
		for (int i = 0; i < points.length; i++) {
			allPix[(int) points[i].getLine()][(int) points[i].getRow()] = points[i];
		}

		for (int i = 0; i < points.length; i++) {
			// ��¼�������������Χ(�߽�)���������ݵ�ľ�����Ϣ
			int[] neighboorCluster = new int[finalClusterNum + 1];
			// ������������
			LinkedList<DataPoint> neighboor = new LinkedList<DataPoint>();
			// ������������Ĵ�С
			int NBsize = 0;
			// �����������������Χ(�߽�)���ݵ�����������һ��ı�ţ���Ҫƽ�������ݵ㽫��ƽ��������(��׼��)
			int NBMax = 0;

			int k = 0;
			for (int j = 0; j < neighboorCluster.length; j++) {
				neighboorCluster[j] = 0;
			}

			if (points[i].getSmoothed() == 0) {
				boolean tag;
				// ����õ�����
				neighboor.add(points[i]);
				do {
					tag = true;
					// ����������ж���temp��¼�����ݵ��ȫ���������ݵ�
					// LinkedList<DataPoint> temp =
					// allEightNeighborFeilds(neighboor
					// .get(k));
					// �����������ж���temp��¼�����ݵ��ȫ���������ݵ�,�����ڶԱ�����
					LinkedList<DataPoint> temp = allFourNeighborFeilds(neighboor
							.get(k));
					for (DataPoint p : temp) {
						// ���뱻���������ͬһ��������ݵ��������
						if (p.getClusterLabel() == points[i].getClusterLabel()) {
							if (!neighboor.contains(p)) {
								neighboor.add(p);
							}
							// �������������߽磬����λ��׼��
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
						// �����Ѿ�������ϲ�����������ѭ��
						tag = false;
					}
					// �����С����ƽ������40������ֹͣ������ѭ����ֹ
				} while (NBsize < smoothValve && tag == true);

				// ����Ҫƽ���ͳ�������ݵ㼰���������е�ƽ������׼�࣬NBsize <
				// 40��ʾ�������ƽ�����޵������㣬neighboor.get(0).getClusterLabel() ==
				// -1��ʾ�������о�Ϊ�����㷨��ʧ��ͼ�����ݵ�

				// ���뷽��:������ٽ���鲢
				if (NBsize < smoothValve
						|| neighboor.get(0).getClusterLabel() == -1) {
					if (NBMax != -1) {
						for (DataPoint p : neighboor) {
							p.setClusterLabel(NBMax);
						}
						// �����������ʧ����������Ϊȱʡֵ
					} else {
						for (DataPoint p : neighboor) {
							p.setClusterLabel(NBMax);
						}
					}
				}
				// ���ڲ���Ҫƽ���ĵ㣬��������Ļ�ͨ�ԣ�����Ҫ�ٴβ����´��ж�
				else {
					for (DataPoint p : neighboor) {
						p.setSmoothed(1);
					}
				}
			}
		}
		// ���ݾ������
		for (int i = 0; i < points.length; i++) {
			allPix[(int) points[i].getLine()][(int) points[i].getRow()] = points[i];
		}
		// �������д���ļ�
		// DataWriter.outPutMatrixtoFile(allPix);
		// ��������ݵ�[R,G,B]д���ļ�
		// DataWriter.outPutPointsRGBtoFile(points);
		// System.out.println("smoothCount= " + smoothCount);
		return points;
	}

	/**
	 * ����allEightNeighborFeilds����ָ�����ݵ���ԭͼ���еİ����� Parameters:point:ָ�������ݵ�
	 * Return:������㼯
	 */
	public LinkedList<DataPoint> allEightNeighborFeilds(DataPoint point) {
		LinkedList<DataPoint> allEightNeighboor = new LinkedList<DataPoint>();
		// ��������
		allEightNeighboor.add(point);
		// ����Ϸ��еĵ�
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
		// ��ӱ��еĵ�
		if (point.getRow() > 1)
			allEightNeighboor.add(allPix[(int) point.getLine()][(int) point
					.getRow() - 1]);

		if (point.getRow() < row_num)
			allEightNeighboor.add(allPix[(int) point.getLine()][(int) point
					.getRow() + 1]);
		// ����·��еĵ�
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
	 * ����allFourNeighborFeilds����ָ�����ݵ���ԭͼ���е������� Parameters:point:ָ�������ݵ�
	 * Return:������㼯
	 */
	public LinkedList<DataPoint> allFourNeighborFeilds(DataPoint point) {
		LinkedList<DataPoint> allFourNeighboor = new LinkedList<DataPoint>();
		// ��������
		allFourNeighboor.add(point);
		// ����Ϸ��ĵ�
		if (point.getLine() > 1)
			allFourNeighboor.add(allPix[(int) point.getLine() - 1][(int) point
					.getRow()]);
		// �����ߵĵ�
		if (point.getRow() > 1)
			allFourNeighboor.add(allPix[(int) point.getLine()][(int) point
					.getRow() - 1]);
		// ����ұߵĵ�
		if (point.getRow() < row_num)
			allFourNeighboor.add(allPix[(int) point.getLine()][(int) point
					.getRow() + 1]);
		// ����·��ĵ�
		if (point.getLine() < line_num)
			allFourNeighboor.add(allPix[(int) point.getLine() + 1][(int) point
					.getRow()]);

		return allFourNeighboor;
	}

	/**
	 * ����initSegments�������ݵ�ľ�������ʼ�����зֿ�,��������ı��,���ݵ㼯�Ϻ��ڿ��ż���
	 * Parameters:points:����ƽ������������� segments:��ʼ���б� Return:��ʼ�ֿ�����
	 */
	public LinkedList<Segment> initSegments(DataPoint[] points,
			LinkedList<Segment> segments) {
		// ����
		int segmentIndex = 0;

		for (int i = 0; i < points.length; i++) {
			// �����С
			int sizeofOnePiece = 0;
			int k = 0;
			if (points[i].getSegmentChecked() == 0) {
				boolean tag;
				// ͼ���е�����
				Segment oneSeg = new Segment(segmentIndex, -1);
				// ����õ�����
				oneSeg.getContentList().add(points[i]);
				points[i].setSegmentChecked(1);
				points[i].setSegmentLabel(segmentIndex);
				// �����ļ��߽�ĳ�ʼ����
				LinkedList<DataPoint> mirrorContentList = oneSeg
						.getContentList();
				LinkedList<DataPoint> mirrorMarginList = oneSeg.getMarginList();
				do {
					tag = true;
					// �����������ж���temp��¼�����ݵ��ȫ���������ݵ�,�������Ŀ�
					LinkedList<DataPoint> temp = allFourNeighborFeilds(mirrorContentList
							.get(k));
					for (DataPoint p : temp) {
						// ���뱻���������ͬһ����(��ͬһ��)�����ݵ��������
						if (p.getClusterLabel() == points[i].getClusterLabel()) {
							// ִ��Ч��
							if (p.getSegmentChecked() == 0) {
								mirrorContentList.add(p);
								p.setSegmentChecked(1);
								p.setSegmentLabel(segmentIndex);
							}
						}
						// ��ĳ�����������(�������)���е�������������߿�,��õ�Ϊ������ı߽��
						else {
							if (!mirrorMarginList.contains(mirrorContentList
									.get(k)))
								mirrorMarginList.add(mirrorContentList.get(k));
						}
					}
					k++;
					sizeofOnePiece = mirrorContentList.size();
					if (k == sizeofOnePiece) {
						// �������Ѿ�������ϲ�����������ѭ��
						tag = false;
					}
				} while (tag == true);
				// ���ô���,���޸�
				oneSeg.setContentList(mirrorContentList);
				oneSeg.setMarginList(mirrorMarginList);
				oneSeg.setNumAll();
				segments.add(oneSeg);
				segmentIndex++;
			}

		}

		// �ڶ���ɨ�裬��ʼ��ÿ��������ڿ���Ϣ
		for (int j = 0; j < segments.size(); j++) {
			// �߽缰����ĳ�ʼ����
			LinkedList<Integer> mirrorNeighborSegLabelList = segments.get(j)
					.getNeighborSegLabelList();
			LinkedList<DataPoint> mirrorMarginList = segments.get(j)
					.getMarginList();

			for (DataPoint marginPoint : mirrorMarginList) {
				LinkedList<DataPoint> temp2 = allFourNeighborFeilds(marginPoint);
				for (DataPoint neighboorPoint : temp2) {
					// ���뱻����㲻ͬ������ݵ��ż�¼����
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
	 * ����generateAuxiliarySegΪÿ�����������������鲢�޸Ŀ�����������Ϣ Parameters:oriSegments:ԭʼ�ֿ�
	 * auxNum:ÿ����չ��ﵽ���¿��� Return:���������ķֿ�����
	 */
	public void generateAuxiliarySeg(LinkedList<Segment> oriSegments, int auxNum) {
		// ��չ��������Ŀ鼯
		LinkedList<Segment> extendedSegments = new LinkedList<Segment>();
		// �����µĿ鼯,���ź��ڿ�ı�Ż������仯,��������Բ���
		for (int i = 0; i < oriSegments.size(); i++) {
			LinkedList<Integer> extendedSegmentsNeighboors = new LinkedList<Integer>();
			extendedSegmentsNeighboors = oriSegments.get(i)
					.getNeighborSegLabelList();
			int beforeExtend = extendedSegmentsNeighboors.size();
			// ÿ������ڿ���ԭ�����ڿ���������
			for (int k = 0; k < beforeExtend; k++) {
				extendedSegmentsNeighboors.add(new Integer(
						extendedSegmentsNeighboors.get(k) * 3 + 1));
				extendedSegmentsNeighboors.add(new Integer(
						extendedSegmentsNeighboors.get(k) * 3 + 2));
				extendedSegmentsNeighboors.set(k, new Integer(
						extendedSegmentsNeighboors.get(k) * 3));
			}
			// ÿ��һ��Ϊ��
			for (int j = 0; j < 3; j++) {
				Segment newSegment = new Segment(-1, -1);
				// newSegment =
				// segments.get(i);ǳ���ƽ��ܵõ�����ָ��ͬһ�飺segments.get(i)��ָ��
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
	 * ����matrixWPriorByAverageValue���ա����ݵ�-���ݿ顱��ֵ��ʽ����ȫͼƬ����������W,
	 * �����Ӿ����W�Ľ����ոþ����ж�Ӧλ�õ�ֵ���й��� Parameters:fatherSegments:ԭʼ�ֿ鼯(�ѽ��и����������)
	 * auxiliaryNum:ÿ����չ��ﵽ���¿��� Return:��ʼ�ֿ����W(segment-point AverageValue based)
	 */
	public void matrixWPriorByAverageValue(LinkedList<Segment> fatherSegments,
			int auxiliaryNum) {
		// ���������
		int sideLength = fatherSegments.size();
		// ������󲻿��ٷ�,�������������(��)�޳�
		if (sideLength > 1) {
			// �����������������
			double[][] arrayW = new double[sideLength][sideLength];
			// ��ʼȫ����0,�Խ����ϼ���ͬһ����չ�õ����ֵ��Ϊ1
			int index = sideLength / auxiliaryNum;
			for (int num = 0; num < index; num++) {
				for (int i = auxiliaryNum * num; i < auxiliaryNum * (num + 1); i++) {
					for (int j = auxiliaryNum * num; j < auxiliaryNum
							* (num + 1); j++) {
						arrayW[i][j] = 1;
					}
				}
			}
			// ���������ϵ�޸ľ����ֵ
			for (int row = 0; row < sideLength; row++) {
				// ����W��ֵW[i,j]=W[j,i]
				double Wij = 0;
				// �뵱ǰ������������ͬһ��(��ͬһ����)�����п�Ŀ��
				LinkedList<Integer> checkForNeighboors = new LinkedList<Integer>();
				// �뵱ǰ�����ڵ����п�Ŀ��
				LinkedList<Integer> allNeighboors = fatherSegments.get(row)
						.getNeighborSegLabelList();
				for (Integer segLabel : allNeighboors) {
					checkForNeighboors.add(segLabel);
				}
				// �����������
				int column = 0;
				// ���ڵ����ݵ㼯��
				LinkedList<DataPoint> checkForSeg = fatherSegments.get(row)
						.getContentList();
				// ��ǰ������
				double segCenterX = fatherSegments.get(row).getXloc();
				double segCenterY = fatherSegments.get(row).getYloc();
				double segCenterZ = fatherSegments.get(row).getZloc();

				for (Integer neighboor : checkForNeighboors) {
					// ȷ����ǰ���ڿ��ھ����е��к�
					for (int columnIndex = 0; columnIndex < fatherSegments
							.size(); columnIndex++) {
						if ((fatherSegments.get(columnIndex).getSegmentLabel() == (int) neighboor)) {
							column = columnIndex;
						}
					}
					// ���ڿ������
					double segNeighboorCenterX = segments.get((int) neighboor)
							.getXloc();
					double segNeighboorCenterY = segments.get((int) neighboor)
							.getYloc();
					double segNeighboorCenterZ = segments.get((int) neighboor)
							.getZloc();
					LinkedList<DataPoint> oneNeighboorSeg = segments.get(
							(int) neighboor).getContentList();
					Wij = 0;
					// ��ǰ�����������ĳһ���ڿ������е�ľ����ƽ����
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
					// ��ǰ���е����е�����ĳһ���ڿ����ĵľ����ƽ����
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
	 * ����matrixWPriorBySegAverageValue���ա����ݿ�-���ݿ顱��ֵ��ʽ����ȫͼƬ����������W,
	 * �����Ӿ����W�Ľ����ոþ����ж�Ӧλ�õ�ֵ���й��� Parameters:fatherSegments:ԭʼ�ֿ鼯(�ѽ��и����������)
	 * auxiliaryNum:ÿ����չ��ﵽ���¿��� Return:��ʼ�ֿ����W(segment-segment AverageValue
	 * based)
	 */
	public void matrixWPriorBySegAverageValue(
			LinkedList<Segment> fatherSegments, int auxiliaryNum) {
		// ���������
		int sideLength = fatherSegments.size();
		// ������󲻿��ٷ�,�������������(��)�޳�
		if (sideLength > 1) {
			// �����������������
			double[][] arrayW = new double[sideLength][sideLength];
			// ��ʼȫ����0,�Խ����ϼ���ͬһ����չ�õ����ֵ��Ϊ1
			int index = sideLength / auxiliaryNum;
			for (int num = 0; num < index; num++) {
				for (int i = auxiliaryNum * num; i < auxiliaryNum * (num + 1); i++) {
					for (int j = auxiliaryNum * num; j < auxiliaryNum
							* (num + 1); j++) {
						arrayW[i][j] = 1;
					}
				}
			}
			// ���������ϵ�޸ľ����ֵ
			for (int row = 0; row < sideLength; row++) {
				// ����W��ֵW[i,j]=W[j,i]
				double Wij = 0;
				// �뵱ǰ������������ͬһ��(��ͬһ����)�����п�Ŀ��
				LinkedList<Integer> checkForNeighboors = new LinkedList<Integer>();
				// �뵱ǰ�����ڵ����п�Ŀ��
				LinkedList<Integer> allNeighboors = fatherSegments.get(row)
						.getNeighborSegLabelList();
				for (Integer segLabel : allNeighboors) {
					checkForNeighboors.add(segLabel);
				}
				// �����������
				int column = 0;

				// ������
				double segCenterX = fatherSegments.get(row).getXloc();
				double segCenterY = fatherSegments.get(row).getYloc();
				double segCenterZ = fatherSegments.get(row).getZloc();

				for (Integer neighboor : checkForNeighboors) {
					// ȷ����ǰ���ڿ��ھ����е��к�
					for (int columnIndex = 0; columnIndex < fatherSegments
							.size(); columnIndex++) {
						if ((fatherSegments.get(columnIndex).getSegmentLabel() == (int) neighboor)) {
							column = columnIndex;
						}
					}
					// ���ڿ������
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
	 * ����matrixWPriorBySegAverageValue���ա�����
	 * ��-���ݵ㡱Э���ʽ����ȫͼƬ����������W,�����Ӿ����W�Ľ����ոþ����ж�Ӧλ�õ�ֵ���й���
	 * Parameters:fatherSegments:ԭʼ�ֿ鼯(�ѽ��и����������) auxiliaryNum:ÿ����չ��ﵽ���¿���
	 * Return:��ʼ�ֿ����W(segment-point Covariance based)
	 */
	public void matrixWPriorByCovariance(LinkedList<Segment> fatherSegments,
			int auxiliaryNum) {
		// ���������
		int sideLength = fatherSegments.size();
		// ��ά��λ����
		Matrix I = new Matrix(new double[][] { { 1, 0, 0 }, { 0, 1, 0 },
				{ 0, 0, 1 } });

		if (sideLength > 1) {
			// �����������������
			double[][] arrayW = new double[sideLength][sideLength];
			// ��ʼȫ����0,�Խ����ϼ���ͬһ����չ�õ����ֵ��Ϊ1
			int index = sideLength / auxiliaryNum;
			for (int num = 0; num < index; num++) {
				for (int i = auxiliaryNum * num; i < auxiliaryNum * (num + 1); i++) {
					for (int j = auxiliaryNum * num; j < auxiliaryNum
							* (num + 1); j++) {
						arrayW[i][j] = 1;
					}
				}
			}
			// ���������ϵ�޸ľ����ֵ
			for (int row = 0; row < sideLength; row++) {
				// ����W��ֵW[i,j]=W[j,i]
				double Wij = 0;
				// Э�������
				Matrix onePointCovR = new Matrix(new double[][] { { 0, 0, 0 },
						{ 0, 0, 0 }, { 0, 0, 0 } });

				Matrix CovarianceSegRow = new Matrix(new double[][] {
						{ 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } });

				// �뵱ǰ������������ͬһ��(��ͬһ����)�����п�Ŀ��
				LinkedList<Integer> checkForNeighboors = new LinkedList<Integer>();
				// �뵱ǰ�����ڵ����п�Ŀ��
				LinkedList<Integer> allNeighboors = fatherSegments.get(row)
						.getNeighborSegLabelList();
				for (Integer segLabel : allNeighboors) {
					checkForNeighboors.add(segLabel);
				}
				// �����������
				int column = 0;
				// ���ڵ����ݵ㼯��
				LinkedList<DataPoint> checkForSeg = fatherSegments.get(row)
						.getContentList();
				// ��ǰ��(L,U,V)��ֵ
				double segCenterX = fatherSegments.get(row).getXloc();
				double segCenterY = fatherSegments.get(row).getYloc();
				double segCenterZ = fatherSegments.get(row).getZloc();

				// ���쵱ǰ���Э�������Matrix CovarianceSegRow
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

				// ��ǰ���ֵ����
				Matrix RowMeanVector = new Matrix(new double[] { segCenterX,
						segCenterY, segCenterZ }, 3);

				for (Integer neighboor : checkForNeighboors) {
					// ȷ����ǰ���ڿ��ھ����е��к�
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

					// �������ڿ��Э�������Matrix CovarianceSegColumn
					// ���ڿ�ĵ㼯
					LinkedList<DataPoint> oneNeighboorSeg = segments.get(
							(int) neighboor).getContentList();
					// ���ڿ������(L,U,V)��ֵ
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

					// ���ڿ��ֵ����
					Matrix ColMeanVector = new Matrix(new double[] {
							segNeighboorCenterX, segNeighboorCenterX,
							segNeighboorCenterX }, 3);
					// ����Wij
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
	 * ����matrixDeclareAndSovle���յ�ǰ���ڵ�������������D��W,����������������ֵ,���õ������ӽڵ�,�ӽڵ���еݹ�,
	 * �������ֿ���ֵʱ�ݹ鷵�� Parameters:fatherSegments:��ǰ���ֿ鼯 Return:�ݹ��ʶ
	 */
	public int matrixDeclareAndSovle(LinkedList<Segment> fatherSegments) {
		// ���ζ��ֵ����������ӿ鼯
		LinkedList<Segment> leftChildSegments = new LinkedList<Segment>();
		LinkedList<Segment> rightChildSegments = new LinkedList<Segment>();
		// ����D,W,D-W,pow(D,-1/2)(D-W)pow(D-1/2)
		Matrix D;
		Matrix W;
		Matrix dInverseSqar;
		Matrix d_Minus_w;
		Matrix equator;

		// ���������
		int sideLength = 0;
		sideLength = fatherSegments.size();

		// ������󲻿��ٷ�,�������������(��)�޳�
		int checkForSingular = 0;
		while (checkForSingular < fatherSegments.size()) {
			boolean noNeighboorInMatrix = true;
			// �뵱ǰ�����ڵ����п�Ŀ��
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
			// �����������������
			double[][] arrayD = new double[sideLength][sideLength];
			double[][] arrayW = new double[sideLength][sideLength];
			// ��ʼȫ����0
			for (int i = 0; i < sideLength; i++)
				for (int j = 0; j < sideLength; j++) {
					arrayD[i][j] = 0;
					arrayW[i][j] = 0;
				}

			// ���ݵ�ǰ����������˳���Ӧ���ҳ�ʼW��ֵ,���쵱ǰ�����W��D
			// �����������
			int row;
			int column;
			for (row = 0; row < sideLength; row++) {
				// ����D�ĶԽ����ϵ�ֵD[k,k]
				double Dk = 0;
				// �����������
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

			// ��������=pow(D,-1/2)��eigenVectors
			eigenVectors = dInverseSqar.times(eigenVectors);

			double[] eigenValues = equator.eig().getRealEigenvalues();

			LinkedList<LinkedList<Segment>> checkForSubDivsion = computeNcutValue(
					eigenValues, eigenVectors, fatherSegments, D, W);

			// �ж���ǰ�鼯�Ƿ�����������ֵ���ֵ����,�������ȡ�ӿ鼯
			if ((fatherSegments.size() >= 2)
					&& (checkForSubDivsion.size() == 2)) {
				leftChildSegments = checkForSubDivsion.get(0);
				rightChildSegments = checkForSubDivsion.get(1);
				// �����ӿ鼯��������
				for (Segment seg : leftChildSegments) {
					seg.setGroupLabel(segmentGroupLabelIndex);
				}
				for (Segment seg : rightChildSegments) {
					int nextSegmentGroupLabelIndex = segmentGroupLabelIndex + 1;
					seg.setGroupLabel(nextSegmentGroupLabelIndex);
				}
				// �����ֽ��д��ȫ�ֱ���segments
				for (Segment seg : segments) {
					if (leftChildSegments.contains(seg)) {
						seg.setGroupLabel(segmentGroupLabelIndex);
					} else if (rightChildSegments.contains(seg)) {
						int nextSegmentGroupLabelIndex = segmentGroupLabelIndex + 1;
						seg.setGroupLabel(nextSegmentGroupLabelIndex);
					}

				}
				// ���ŵ���
				segmentGroupLabelIndex = segmentGroupLabelIndex + 2;
				// �����ӿ鼯�ݹ����
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
	 * ����computeNcutValueÿ�������еĶ��ֲ�����Ѱ�����Ż��֣��Ҳ����򷵻ؿռ�
	 * Parameters:eigenValues:�����������ֵ���� eigenVectors:�������������������
	 * fatherSegments:���鼯 D_of_FatherSeg:����D���� W_of_FatherSeg:����W���� Return:�����ӿ�
	 */
	public LinkedList<LinkedList<Segment>> computeNcutValue(
			double[] eigenValues, Matrix eigenVectors,
			LinkedList<Segment> fatherSegments, Matrix D_of_FatherSeg,
			Matrix W_of_FatherSeg) {

		LinkedList<LinkedList<Segment>> result = new LinkedList<LinkedList<Segment>>();
		// ���������Ż���
		double minNcutValue = Double.MAX_VALUE;
		// ����ֵ��С�������������¼��Ӧ�Ŀ��
		double[] segLabelSortByEigenValues = new double[eigenValues.length];
		// ���п��ܵĶ��ֲ���
		int possibleDivision = 0;
		// �ڶ�С������ֵ��Ӧ������������Ӧ�Ŀ������
		Matrix secondSmallestEigenVectors = new Matrix(eigenValues.length, 1);
		// �ڶ�С������ֵ��Ӧ����������
		Matrix secondSmallestEigenVectorsMirror = new Matrix(
				eigenValues.length, 1);

		// ����ֵ����
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

		// ����ֵ��С�������������¼��Ӧ�Ŀ��,�ҵ��ڶ�С������ֵ��Ӧ����������
		for (int i = 0; i < segLabelSortByEigenValues.length; i++) {
			// ���ж������ֵ��ȵ�����������״�����ԭ��
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

		// �����������е�ֵ��С�����˳������
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
		// ��¼���ھ����еĺ�
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
			// ĳһ�ֶ��ַ�������ӿ鼯,��ѡ��
			LinkedList<Segment> leftChildSegmentsCandidate = new LinkedList<Segment>();
			LinkedList<Segment> rightChildSegmentsCandidate = new LinkedList<Segment>();
			// ���λ��ֵ�Ncutֵ����صļ����м�ֵ
			double certainNcutValue = 0;

			double cutLR = 0;
			double assocLF = 0;
			double cutRL = 0;
			// double cutRL1 = 0;
			double assocRF = 0;
			// ���쵥�λ��ֵ��ӿ鼯
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

		// ���ַ��ٽ��ж�
		if (minNcutValue >= nCutValve) {
			result.remove(1);
			result.remove(0);
		}

		return result;
	}

	/**
	 * ����segmentsArrange�����е����ݵ㰴�����Ź鲢 Parameters:allSegmentsGrouped:���п�����Ľ����
	 * allPoints:���е㼯
	 * 
	 * Return:�ֿ����,�޸�������ά����ֵ����������ݵ㼯��
	 */

	public DataPoint[] segmentsArrange(LinkedList<Segment> allSegmentsGrouped,
			DataPoint[] allPoints) {
		// ���������ݵ�
		DataPoint[] classifiedPoints = new DataPoint[allPoints.length];
		LinkedList<Group> allGroups = new LinkedList<Group>();

		for (int i = 0; i < allSegmentsGrouped.size(); i++) {
			int included = -1;
			// �ÿ��Ƿ������ѽ���������
			for (int j = 0; j < allGroups.size(); j++) {
				if (allGroups.get(j).getGroupLabel() == allSegmentsGrouped.get(
						i).getGroupLabel()) {
					included = j;

				}
			}
			// ���������µķ���
			if (included == -1) {
				Group newGroup = new Group(allSegmentsGrouped.get(i)
						.getGroupLabel());
				// Ϊ�µ��齨�����ݵ㼯��,for{...}�ڴ˴���ָֹ�봫��(ǳ����)
				LinkedList<DataPoint> segentsContentMirrow = new LinkedList<DataPoint>();
				for (int all = 0; all < allSegmentsGrouped.get(i)
						.getContentList().size(); all++) {
					segentsContentMirrow.add(allSegmentsGrouped.get(i)
							.getContentList().get(all));
				}
				newGroup.setContentList(segentsContentMirrow);
				newGroup.setNumAll();
				allGroups.add(newGroup);
				// ����������з���
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

		// ���ݷ������޸����ݵ�������ռ�ֵ
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
	 * ����ָ����������������ľ���ذ�������������(��Χ��) cluster��ָ���������� return �����ľ����
	 */
	public Cluster searchforCluster(Cluster cluster) {
		// ����������ķ�Χ��
		LinkedList<Coord> list = new LinkedList<Coord>();
		list = cluster.contentList;

		while (true) {
			int loopS = list.size();
			// �����ռ���������·���
			List<Coord> upDown = checkUpAndDown(list);
			if (upDown.size() != 0) {
				for (Coord c : upDown) {
					if (!list.contains(c))
						list.add(c);
				}
				upDown.clear();
			}
			// �����ռ���������ҷ���
			List<Coord> leftRight = checkLeftAndRight(list);
			if (leftRight.size() != 0) {
				for (Coord c : leftRight) {
					if (!list.contains(c))
						list.add(c);
				}
				leftRight.clear();
			}
			// �����ռ������ǰ����
			List<Coord> frontBack = checkFrontAndBack(list);
			for (Coord c : frontBack) {
				if (!list.contains(c))
					list.add(c);
			}
			frontBack.clear();
			int loopE = list.size();
			// ������ķ�Χ��������ʱ����������
			if (loopE == loopS)
				break;
		}
		cluster.contentList = list;
		return cluster;
	}

	/**
	 * �ڿռ���������ҷ��������������������������� list������ĳ�ʼ��Χ�� return result�����ҷ�������������Χ��
	 * 
	 */
	List<Coord> checkLeftAndRight(List<Coord> list) {
		List<Coord> result = new LinkedList<Coord>();
		for (Coord gnumber : list) {
			int x = gnumber.getX();
			int y = gnumber.getY();
			int z = gnumber.getZ();
			// �߽�����������
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
	 * �ڿռ���������·��������������������������� list������ĳ�ʼ��Χ�� return result�����·�������������Χ��
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
	 * �ڿռ������ǰ���������������������������� list������ĳ�ʼ��Χ�� return result��ǰ��������������Χ��
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
	 * ����󷽵���������󷽵�grid��ƫ��ֵ�ݶȴ��ڵ�ǰ�����ƫ��ֵ�ݶȣ��ͷ���true�����򷵻�false x ��ǰ�����X���� y
	 * ��ǰ�����Y���� z ��ǰ�����Z���� list ��ǰ�����������������б�(����Χ��) return �󷽵������Ƿ��������
	 */
	private boolean checkLeftGrid(int x, int y, int z, List<Coord> list,
			List<Coord> result) {
		if (x == 0)
			return false;
		if (steps[x - 1][y][z].getAlreadysearched() != 0)
			return false;
		// ����������ݶ�ֵ
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
	 * ����ҷ�����������ҷ���grid��ƫ��ֵ�ݶȴ��ڵ�ǰ�����ƫ��ֵ�ݶȣ��ͷ���true�����򷵻�false x ��ǰ�����X���� y
	 * ��ǰ�����Y���� z ��ǰ�����Z���� list ��ǰ�����������������б�(����Χ��) return �ҷ��������Ƿ��������
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
	 * ����Ϸ�����������Ϸ���grid��ƫ��ֵ�ݶȴ��ڵ�ǰ�����ƫ��ֵ�ݶȣ��ͷ���true�����򷵻�false x ��ǰ�����X���� y
	 * ��ǰ�����Y���� z ��ǰ�����Z���� list ��ǰ�����������������б�(����Χ��) return �Ϸ��������Ƿ��������
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
	 * ����·�����������·���grid��ƫ��ֵ�ݶȴ��ڵ�ǰ�����ƫ��ֵ�ݶȣ��ͷ���true�����򷵻�false x ��ǰ�����X���� y
	 * ��ǰ�����Y���� z ��ǰ�����Z���� list ��ǰ�����������������б�(����Χ��) return �·��������Ƿ��������
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
	 * ���ǰ�����������ǰ����grid��ƫ��ֵ�ݶȴ��ڵ�ǰ�����ƫ��ֵ�ݶȣ��ͷ���true�����򷵻�false x ��ǰ�����X���� y
	 * ��ǰ�����Y���� z ��ǰ�����Z���� list ��ǰ�����������������б�(����Χ��) return ǰ���������Ƿ��������
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
	 * ���󷽵���������󷽵�grid��ƫ��ֵ�ݶȴ��ڵ�ǰ�����ƫ��ֵ�ݶȣ��ͷ���true�����򷵻�false x ��ǰ�����X���� y
	 * ��ǰ�����Y���� z ��ǰ�����Z���� list ��ǰ�����������������б�(����Χ��) return �󷽵������Ƿ��������
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
	 * ���������ƫ���ƴ�С��ϵ��������ʼ�������� return ��ʼ�����������꼯
	 */
	ArrayList<Coord> markCandidatePoints() {
		// �����С�����xƫ����yƫ����zƫ��
		countPotential();
		// �ظ��������ҳ���ʼ�������ĵ㼯��
		ArrayList<Coord> xlist = getCanPointsOnX();
		ArrayList<Coord> ylist = getCanPointsOnY();
		ArrayList<Coord> zlist = getCanPointsOnZ();

		ArrayList<Coord> centers = new ArrayList<Coord>();
		// �󼯺ϵĽ����õ���ʼ��������
		for (Coord c : xlist) {
			if (ylist.contains(c) && zlist.contains(c)) {
				if (steps[c.getX()][c.getY()][c.getZ()].pnumber > 0)
					centers.add(c);
			}
		}
		return centers;
	}

	/**
	 * ��X�᷽��������ѡ�������� return X�᷽���ʼ�����������꼯
	 */
	private ArrayList<Coord> getCanPointsOnX() {
		ArrayList<Coord> points = new ArrayList<Coord>();
		// �ȽϿռ����������Xƫ������ѡȡ��ѡ��������
		for (int y = 0; y < K; y++) {
			for (int z = 0; z < K; z++) {
				for (int x = 1; x < K; x++) {
					double dpx1 = steps[x - 1][y][z].getXdp();
					double dpx2 = steps[x][y][z].getXdp();
					if (dpx1 == 0 && dpx2 < 0)
						// dpx1�Ǻ�ѡ��������
						points.add(new Coord(x - 1, y, z));
					else if (dpx1 > 0 && dpx2 == 0)
						// dpx2�Ǻ�ѡ��������
						points.add(new Coord(x, y, z));
					// dpx1,dpx2���Ǻ�ѡ��������
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
	 * ��Y�᷽��������ѡ�������� return X�᷽���ʼ�����������꼯
	 */
	private ArrayList<Coord> getCanPointsOnY() {
		ArrayList<Coord> points = new ArrayList<Coord>();
		// �ȽϿռ����������Yƫ������ѡȡ��ѡ��������
		for (int z = 0; z < K; z++) {
			for (int x = 0; x < K; x++) {
				for (int y = 1; y < K; y++) {
					double dpy1 = steps[x][y - 1][z].getYdp();
					double dpy2 = steps[x][y][z].getYdp();
					if (dpy1 == 0 && dpy2 < 0)
						// dpy1�Ǻ�ѡ��������
						points.add(new Coord(x, y - 1, z));
					else if (dpy1 > 0 && dpy2 == 0)
						// dpy2�Ǻ�ѡ��������
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
	 * ��Z�᷽��������ѡ�������� return X�᷽���ʼ�����������꼯
	 */
	private ArrayList<Coord> getCanPointsOnZ() {
		ArrayList<Coord> points = new ArrayList<Coord>();
		// �ȽϿռ����������Zƫ������ѡȡ��ѡ��������
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 1; z < K; z++) {
					double dpy1 = steps[x][y][z - 1].getZdp();
					double dpy2 = steps[x][y][z].getZdp();
					if (dpy1 == 0 && dpy2 < 0)
						// dpy1�Ǻ�ѡ��������
						points.add(new Coord(x, y, z - 1));
					else if (dpy1 > 0 && dpy2 == 0)
						// dpy2�Ǻ�ѡ��������
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
	 * ����ÿ��С�����xƫ����yƫ����zƫ��
	 */
	private void countPotential() {
		// ����С�����ʼ��������
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
					// С�������Ե��ۼ�
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
						// �մ����������������Ϊ��������
						BigGrid grid = new BigGrid(xmin + kx * IFP * dx + IFP
								* dx / 2, ymin + ky * IFP * dy + IFP * dy / 2,
								zmin + kz * IFP * dz + IFP * dz / 2, 0, 0, pall);
						grids[kx][ky][kz] = grid;
						allBigGrid.add(grids[kx][ky][kz]);
					}

				}
			}
		}

		// �ǿ�С��������꼯
		ArrayList<Coord> notEmptyGridCoord = new ArrayList<Coord>();
		for (int y = 0; y < K; y++) {
			for (int z = 0; z < K; z++) {
				for (int x = 0; x < K; x++) {
					if (steps[x][y][z].pnumber >= 1)
						notEmptyGridCoord.add(new Coord(x, y, z));
				}
			}
		}
		// �ǿմ���������꼯
		ArrayList<Coord> notEmptyBigGridCoord = new ArrayList<Coord>();
		for (int y = 0; y < BK; y++) {
			for (int z = 0; z < BK; z++) {
				for (int x = 0; x < BK; x++) {
					if (grids[x][y][z].pnumber >= 1)
						notEmptyBigGridCoord.add(new Coord(x, y, z));
				}
			}
		}

		// ���������,С�������Ϣ���ı��ļ�
		String s = new String();
		try {
			File fo = new File("C:\\Users\\YinJF\\Work\\CenterOfGravity.txt");
			File fb = new File(
					"C:\\Users\\YinJF\\Work\\CenterOfGravityBigGrid.txt");
			if (fo.exists()) {
				System.out.println("�ļ�og����,ɾ����ǰ��...");
				fo.delete();
			}
			if (fb.exists()) {
				System.out.println("�ļ�bg����,ɾ����ǰ��...");
				fb.delete();
			}
			if (fo.createNewFile()) {
				System.out.println("�ļ�og�����ɹ���");
			} else {
				System.out.println("�ļ�og����ʧ�ܣ�");
			}
			if (fb.createNewFile()) {
				System.out.println("�ļ�bg�����ɹ���");
			} else {
				System.out.println("�ļ�bg����ʧ�ܣ�");
			}

			BufferedWriter outputfb = new BufferedWriter(new FileWriter(fb));
			// �������������
			for (BigGrid bigGrid : allBigGrid) {
				System.out.println("------���д�����: ["
						+ (int) Math.floor((bigGrid.getXloc() - xmin)
								/ (dx * IFP))
						+ "]["
						+ (int) Math.floor((bigGrid.getYloc() - ymin)
								/ (dy * IFP))
						+ "]["
						+ (int) Math.floor((bigGrid.getZloc() - zmin)
								/ (dz * IFP)) + "]" + " �е�"
						+ bigGrid.getWeight() + " ��------");
				s = Math.round(bigGrid.getXloc() * 100) / 100.0 + " "
						+ Math.round(bigGrid.getYloc() * 100) / 100.0 + " "
						+ Math.round(bigGrid.getZloc() * 100) / 100.0 + " "
						+ +bigGrid.pnumber + "\r\n";
				// �����������
				outputfb.write(s);
			}
			outputfb.close();
			// С�������ļ��������ݵ������
			BufferedWriter outputfo = new BufferedWriter(new FileWriter(fo));
			for (Grid steps : notEmptyGrid) {
				System.out.println("------�ǿ�С���� : ["
						+ (int) Math.floor((steps.getXloc() - xmin) / dx)
						+ "]["
						+ (int) Math.floor((steps.getYloc() - ymin) / dy)
						+ "]["
						+ (int) Math.floor((steps.getZloc() - zmin) / dz) + "]"
						+ " �е�" + steps.getWeight() + " ��------");
				s = Math.round(steps.getXloc() * 100) / 100.0 + " "
						+ Math.round(steps.getYloc() * 100) / 100.0 + " "
						+ Math.round(steps.getZloc() * 100) / 100.0 + " "
						+ steps.getWeight() + "\r\n";
				// �����������
				outputfo.write(s);
			}
			outputfo.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		// �������С����Ȩ�����ı��ļ�
		String sw = new String();
		try {
			File fo = new File("C:\\Users\\YinJF\\Work\\GridWeight.txt");

			if (fo.exists()) {
				System.out.println("�ļ�GridWeight����,ɾ����ǰ��...");
				fo.delete();
			}
			if (fo.createNewFile()) {
				System.out.println("�ļ�GridWeight�����ɹ���");
			} else {
				System.out.println("�ļ�GridWeight����ʧ�ܣ�");
			}
			BufferedWriter outputfo = new BufferedWriter(new FileWriter(fo));

			// ����xƫ����yƫ����zƫ��
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

						// �����ݳ��ƺ����Ķ��壬�����С�������ֵ��ƫ����
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
									// ��������
									double dist2 = dx * dx + dy * dy + dz * dz;
									// �ռ����
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
	 * ��Coord����ռ�������������
	 */
	public class Coord {
		// X��Y��Z������ֵ
		private int x;
		private int y;
		private int z;

		// �����������������
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
		// �ж��������������Ƿ����
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
	 * ��BigGrid����ռ���������
	 */
	private class BigGrid {
		// �������X��Y��Z������ֵ
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		private double spacialX = -1;

		private double spacialY = -1;

		// �������ڰ��������ݵ���
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
	 * ��Cluster��������������
	 */
	public class Cluster {
		// ������
		int clusterLabel;

		// ��������
		double xloc = -1;

		double yloc = -1;

		double zloc = -1;

		// �������������ֵͳ��
		double xAll = 0;

		double yAll = 0;

		double zAll = 0;

		double distanceSQ = 0;

		// �����а��������ݵ�����
		double numAll = 0;

		// ���෶Χ��(������)�������꼯
		LinkedList<Coord> contentList;

		// ������������
		Coord clusterCenter;

		// ��ʼ������������
		Coord initialClusterCenter;

		// �Ҿ����ʼ���ĵ�������
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

		// �ϲ��������������
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
	 * ��Segment����ͼ��ֿ����
	 */
	private class Segment {
		// ����
		int segmentLabel;

		// ����
		int groupLabel;

		// �����������ռ���������
		double xloc = -1;

		double yloc = -1;

		double zloc = -1;

		// �����ռ������ֵͳ��
		double xAll = 0;

		double yAll = 0;

		double zAll = 0;;

		// ���а��������ݵ�����
		double numAll = 0;

		// �������ݵ㼯��
		LinkedList<DataPoint> contentList;

		// ���ڱ߽�㼯��
		LinkedList<DataPoint> marginList;

		// ���ڿ�ı�ż���
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
	 * ��Group����ͼ�������,����ĵļ���
	 */
	private class Group {
		// ����
		int groupLabel;

		// ס���������ռ���������
		double xloc = -1;

		double yloc = -1;

		double zloc = -1;

		// �����ռ������ֵͳ��
		double xAll = 0;

		double yAll = 0;

		double zAll = 0;

		// ���а��������ݵ�����
		double numAll = 0;

		// �������ݵ㼯��
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
	 * ��Grid����ռ�С�������
	 */
	public class Grid {
		// �������������ֵͳ��
		private double xall = 0;

		private double yall = 0;

		private double zall = 0;

		private double spacialXall = 0;

		private double spacialYall = 0;

		// С�����X��Y��Z������ֵ,����������ƽ��ֵ
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// С����ռ�����ƽ��ֵ
		private double spacialX = -1;

		private double spacialY = -1;

		// ����ƫ����
		private double xdp = -1;

		private double ydp = -1;

		private double zdp = -1;

		// С������������ݵ�����
		private int pnumber = 0;

		// С����������������ݵ㼯
		private ArrayList<DataPoint> list = new ArrayList<DataPoint>();

		// С�����Ƿ�������
		private int alreadysearched = 0;

		// С�������ֵ
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
