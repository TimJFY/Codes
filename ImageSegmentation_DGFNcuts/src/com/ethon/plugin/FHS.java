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
 * FHS ��������ͼ��ָ��㷨
 * 
 * @author liwei
 * 
 */

public class FHS {

	// ͼƬ���п���п�
	int line_num, row_num;

	// С����������յ�
	public double xmin, ymin, zmin, xmax, ymax, zmax;

	// FHS�㷨��ʼ����ʱ��
	private long start;

	// FHS�㷨���н���ʱ��
	private long end;

	// С�������͵ľ���
	public Cells[][][] cells;

	// ����ϵͳ�����зǿ�С���񼯺�
	private ArrayList<Cells> notEmptyCells = new ArrayList<Cells>();

	// ����ϵͳ�зǿ�����ĸ���
	private int notEmptyCellsNum = 0;

	// ��ʼ�ܶȺ�����ÿ��С�����Ӧһ����ʼ�ܶȺ�����
	public double[][][] initial;

	// ����Ӧ�ܶȺ�����ÿ��С�����Ӧһ������Ӧ�ܶȺ�����
	public double[][][] adaptive;

	// ��������(���ڼ���С����ı߳�)
	public double[][][] varFactor;

	// �ɱ�������ڼ�������Ӧ�ܶȵ�һ������ֵ��
	public double[][][] varBandwidth;

	// h(x)��d�η�
	public double[][][] varHd;

	// �۳˻������ڼ����ʼ�ܶȺ����ļ���ƽ��ֵ��
	private BigDecimal multiply;
	
	// �۳˻������ڼ�������Ӧ�ܶȺ����ļ���ƽ��ֵ��
	private BigDecimal multiply2;
	
	// ��ʼ�ܶȺ����ļ���ƽ��ֵ
	private double geoMean = 0;

	// ����Ӧ���ܶȺ����ļ���ƽ��ֵ
	private double geoMean2 = 0;

	// ����ˮƽ
	private double noiseLevel = 0;

	// ����ϵͳÿ�߰���С���������
	private int K;

	// �̶�����h0(���û�����)
	private int h0;

	// ro���ӵĳ�ʼֵ�Ƕ��٣�����
	private double ro = 1.0;

	// �̶�����(Ҳ����ÿ��С����ı߳�)
	private double width;

	// ����ԭʼͼ�����ݵ�
	DataPoint[] RGBPoints;

	// ���о���ת��ΪLUV����ϵ�����ݵ�
	DataPoint[] LUVPoints;

	// �þ����У��������ݵ㰴����ԭͼ���е�����(x,y)��¼�ڶ�Ӧ��λ�ã�ԭͼ��Ϊline_num*row_num
	DataPoint[][] allPix;

	// ��¼�������ݵĵ�
	private DataPoint[] points1 = null;

	// ԭʼͼ�����ݵ������
	private int allPointsNum = 0;

	// ��ʼ������Ŀ
	private int initalClusterNum;

	// �м������Ŀ
	private int middleClusterNum;
	
	// ���վ�����Ŀ
	private int finalClusterNum;

	// ��ʼ�ֿ���Ŀ
	private int initSegmentNum;

	// ���շֿ���Ŀ
	private int finalSegmentNum;

	// ƽ����ֵ
	private double smoothValve;

	// ���ƽ����ͼ�������еĿ�,Ҳ��ͼ�񰴿黮�ֵ����ս��
	LinkedList<Segment> segments = new LinkedList<Segment>();

	// ����������ı��
	private int segmentGroupLabelIndex = 0;

	/**
	 * ���캯��FHS�� �����û��������h0��������ԭʼͼ�����ݵ�ӳ��������ϵͳ����ȷ����ʼ��ȫ��С������� Parameters:
	 * h0:�̶������û����룩 line_num:ͼƬ���п� row_num:ͼƬ���п� smoothValve:ƽ����ֵ
	 */
	public FHS(int h0, int line_num, int row_num) {

		// ��ʼ��
		start = System.currentTimeMillis();

		this.h0 = h0;
		this.line_num = line_num;
		this.row_num = row_num;

		this.smoothValve = 40;

		allPix = new DataPoint[line_num + 1][row_num + 1];

		// ��ʼ�����ݵ㣨ԭʼͼ�����ݵ�ļ��� ��
		RGBPoints = DataBase.getInstance().getPoints();
		points1 = RGBPoints;

		// ���������ݳ�
		DataBase db = new DataBase();
		// ����ȫ�ֱ���points1
		LUVPoints = db.new_public_points(points1);

		// RGBֵת��ΪLUV������浽LUVPoint������
		LUVPoints = RGBToLUV(LUVPoints);

		// ͳ��LUV�ռ����ص�ĸ�����������LUVPoints���鳤�ȣ�
		allPointsNum = LUVPoints.length;

		// Ϊxmin��xmin��zmin��xmax��ymax��zmax��ʼ��
		xmin = Double.MAX_VALUE;
		ymin = Double.MAX_VALUE;
		zmin = Double.MAX_VALUE;

		xmax = Double.MIN_VALUE;
		ymax = Double.MIN_VALUE;
		zmax = Double.MIN_VALUE;

		// ����ÿ�����ݵ�
		for (DataPoint p : LUVPoints) {
			// ���ݵ���LUV�ռ��x��y��z����ֵ
			double x = p.getCoord_X();
			double y = p.getCoord_Y();
			double z = p.getCoord_Z();
			
			//4��23�����������
			//p.setClusterLabel(-1);
			
			// ������ά����ϵͳ��ά����㣨xmin��ymin��zmin��
			if (xmin > x)
				xmin = x;
			if (ymin > y)
				ymin = y;
			if (zmin > z)
				zmin = z;

			// ������ά����ϵͳ��ά���յ㣨xmax��ymax��zmax��
			if (xmax < x)
				xmax = x;
			if (ymax < y)
				ymax = y;
			if (zmax < z)
				zmax = z;
		}// foreachѭ������

		// ���㴰��Ҳ����С����ı߳���
		this.width = h0 / 3.0;

		System.out.println("С����Ŀ��width  = " + width);
		double dx = xmax - xmin;
		double dy = ymax - ymin;
		double dz = zmax - zmin;
		// �Ƚ�dx��dy��dz��ȡ���������Ǹ�
		double minusMax = 0;
		double temp = dx > dy ? dx :dy;
		minusMax = temp > dz ? temp : dz;
		// ��������ϵͳÿ�߰���С���������
		/**
		 * ����public static double ceil(double a) ������С�ģ���ӽ��������double
		 * ֵ����ֵ���ڵ��ڲ�����������ĳ��������
		 */
		//this.K = (int) Math.ceil((xmax - xmin) / width);
		this.K = (int) Math.ceil(minusMax / width);
		System.out.println(" ÿ�����������С����ĸ���K = " + K);
		
		// С�������;����ʼ��
		cells = new Cells[K][K][K];

		// ��ÿ�����ݵ���䵽��Ӧ��С�����У�������һ����С����Ϊ��
		for (DataPoint p : LUVPoints) {
			// ����:ÿһ�����ݵ��������ͬ����
			// ���ݵ������ֵ
			double x = p.getCoord_X();
			double y = p.getCoord_Y();
			double z = p.getCoord_Z();

			/**
			 * public static double floor(double a) �������ģ���ӽ��������double
			 * ֵ����ֵС�ڵ��ڲ�����������ĳ�������� cells[xl][yl][zl] Ϊ���ݵ㣨x,y,z����Ӧ������
			 */
			// ���ݵ������������
			int xl = (int) Math.floor((x - xmin) / width);
			int yl = (int) Math.floor((y - ymin) / width);
			int zl = (int) Math.floor((z - zmin) / width);
			//���⣺��Ȼ������zl > 11�����
			System.out.println("cells[" + xl + "][" + yl + "][" + zl +"]");

			if (xl == K)
				xl = K - 1;
			if (yl == K)
				yl = K - 1;
			if (zl == K)
				zl = K - 1;

			// �����ݵ���ӵ�С�������֮��
			if (cells[xl][yl][zl] == null)
				//cells[xl][yl][zl] = new Cells();
				cells[xl][yl][zl] = new Cells(xl, yl, zl);

			cells[xl][yl][zl].add(p);

			// ��¼���зǿ�С����
			/**
			 * public boolean contains(Object o) ������б��а���ָ����Ԫ�أ��򷵻� true��
			 */

			/**
			 * public boolean add(E e) ��ָ����Ԫ����ӵ����б��β����
			 */
			if (!notEmptyCells.contains(cells[xl][yl][zl]))
				notEmptyCells.add(cells[xl][yl][zl]);

		}// foreachѭ������
		System.out.println("���ݵ���䵽��ӦС������ɣ�");
		
		// ȷ������С���������ʼ��
		// ͳ��С���������
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					// С����ļ�����������
					double xloc = xmin + width / 2 + width * x;
					double yloc = ymin + width / 2 + width * y;
					double zloc = zmin + width / 2 + width * z;

					Cells cell = cells[x][y][z];
					// ��ʼ����Щ�յ�С����
					if (cell == null) {						
						cell = new Cells();
						// С�������������Ϊ�伸������
						cell.setCenter(xloc, yloc, zloc);
						cells[x][y][z] = cell;
					}
				}
			}
		}// forѭ������
		System.out.println("С�����ʼ����ɣ�");		
		
		// ��ʼ����ʼ�ܶȺ���������Ӧ�ܶȺ���
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
		//��ʼ��multiply
		multiply = new BigDecimal(Double.toString(1.0));
		System.out.println("��ʼ���� multiply = " + multiply);
		
		/**
		 * ��3�� ��ʼ�ܶȹ���
		 */
		// ��ʼ�ܶȹ���
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {

					// ��ʼ�ܶȹ���,ÿ��С������һ����ʼ�ܶ�
					// ÿ��С����ĳ�ʼ�ܶ� = ÿ��С�����ڰ��������ݵ����/LUV�ռ����ݵ������
					System.out.println("��С���������ݵ�ĸ�����"
							+ cells[x][y][z].getpNum());
					System.out.println("LUV��ɫ�ռ����ݵ��ܸ�����" + LUVPoints.length);
					initial[x][y][z] = (double) (cells[x][y][z].getpNum() )/ LUVPoints.length;
					System.out.println("��ʼ�ܶ�ֵ��" + initial[x][y][z]);

					// ��ʼ�ܶ��۳�
					// ע�⣺ȥ����Щ�յ�С����ֻͳ����Щ�ǿյ�С��������
					//if (!(Math.abs(initial[x][y][z] - 0.0) < 1e-6))
					if((initial[x][y][z] - 0.0) != 0)
					{
						// ÿ����һ����ʼ�ܶ�ֵ��Ϊ0�����񣬷ǿ�����������+1
						notEmptyCellsNum++;// ͳ�Ʒǿ�����ĸ���
						System.out.println("�ǿյ���������� "+notEmptyCellsNum);						
					}
				}
			}
		}// forѭ������
		
		System.out.println("�ǿ�С����ĸ�����" + notEmptyCellsNum);
		
		for(int x = 0; x < K; x++){
			for(int y = 0; y < K; y++){
				for(int z = 0; z < K; z++){
					if((initial[x][y][z] - 0.0) != 0){
						
						//double���ͱ�����ת��ΪString���ͣ�Ȼ��ͨ��BigDecimal�Ĺ��캯��תΪBigDecimal���ͣ����м���
						/*
						 * BigDecimal bd = new BigDecimal(string);
						 * ���BigDecimal�÷��ĵ�
						 */
						double v1 = initial[x][y][z];
						double multiplyDou = 1.0;
						BigDecimal initial = new BigDecimal(Double.toString(v1));
						//�ȼ���ÿһ����ʼ�ܶȵ�notEmptyCellsNum�η���ֵ
						multiplyDou = StrictMath.pow(initial.doubleValue(), 1.0/notEmptyCellsNum);
						//�۳˼���,�õ�����ƽ��ֵ
						multiply = new BigDecimal(Double.toString(multiplyDou));
						multiply = multiply.multiply(initial);
						
						System.out.println("ÿһ���۳˵Ľ��multiply  = " + multiply.doubleValue());
					}//if����
					
				}
			}
		}//forѭ������	

		System.out.println("multiply  = " + multiply.doubleValue());
		
		// �����ʼ�ܶȺ����ļ���ƽ��ֵ
		geoMean = multiply.doubleValue();
		System.out.println("��ʼ�ܶȺ����ļ���ƽ��ֵ�ǣ�geoMean = " + geoMean);
		
		/**
		 * ��4�� ����Ӧ���ܶȹ���
		 */
		// ����Ӧ���ܶȹ���
		//ͳ�Ƹ���
		int counters = 1;
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					counters++;
					double plus = 0;// �ۼӺ�

					//ֻ���㵱ǰС������Χ��26��С��������
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
					
					// ������xΪ��׼���x������Ӧ���ܶȺ���ֵ(�����⣺�����λ������������ĸ������Ƿǿ�����ĸ�������)					
					adaptive[x][y][z] = plus / notEmptyCellsNum;
					System.out.println("����Ӧ�ܶȺ���: " + adaptive[x][y][z]);
				}
			}
		}// ����forѭ������

		
		//��ʼ��multiply2
		multiply2 = new BigDecimal(Double.toString(1.0));
		System.out.println("��ʼ���� multiply2 = " + multiply2);
		
		// ��������Ӧ�ܶȺ���ֵ
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					//ֻͳ�Ʒǿ�����
					if ((adaptive[x][y][z] - 0.0) != 0) {
						
						double v2 = adaptive[x][y][z];
						
						BigDecimal adaptive = new BigDecimal(Double.toString(v2));
						
						double multiplyDou = 1.0;
						//�ȼ���ÿһ����ʼ�ܶȵ�notEmptyCellsNum�η���ֵ				
						multiplyDou = StrictMath.pow(adaptive.doubleValue(), 1.0 / notEmptyCellsNum);
						
						//�۳˼���,�õ�����ƽ��ֵ
						multiply2 = new BigDecimal(Double.toString(multiplyDou));
						multiply2 = multiply2.multiply(adaptive);
						
						System.out.println("ÿһ���۳˵Ľ��multiply2  = " + multiply2.doubleValue());												
					}// if����
				}
			}
		}// forѭ������

		// ����Ӧ�ܶȺ����ļ���ƽ��ֵ
		geoMean2 = multiply2.doubleValue();
		System.out.println("����Ӧ�ܶȺ����ļ���ƽ��ֵ :" + geoMean2);

	}// ���캯������
	
	
	

	/**
	 * ��Cells
	 * ����ռ�С�������
	 * 
	 */
	public class Cells {

		// X��Y��Z������ֵ
		private int x;

		private int y;

		private int z;

		// �����������ֵͳ�ƣ��ۼӺ�
		private double xall = 0;

		private double yall = 0;

		private double zall = 0;

		// ��Ӧ��ͼ�������ֵ�ۼӺ�
		private double spacialXall = 0;

		private double spacialYall = 0;

		// С��������x��y��z�������ֵ�����������ݵ������ƽ��ֵ��
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// ��Ӧ��ͼ�������ƽ��ֵ(���ݵ��Ӧ��ͼ���x��y����ֵ������ƽ��ֵ)
		private double spacialX = -1;

		private double spacialY = -1;

		// x��y��z �������������ƫ����
		private double xdp = -1;

		private double ydp = -1;

		private double zdp = -1;

		// С�����а��������ݵ�ĸ���
		private int pNum;

		// �������������ľ�����
		private int clusterNum;

		// С����ĳ�ʼ�ܶ�ֵ
		private double density1 = -1;

		// С���������Ӧ�ܶ�ֵ
		private double density2 = -1;

		// С�������ֵ
		private double potential = -1;

		// С����������������ݵ㼯
		private ArrayList<DataPoint> list = new ArrayList<DataPoint>();

		// С�����Ƿ�������
		private int alreadysearched = 0;

		// С�����Ƿ񱻺ϲ���
		private int alreadyMerged = 0;
		
		// �������Ĺ��캯��
		public Cells(int x, int y, int z) {
			this.setX(x);
			this.setY(y);
			this.setZ(z);
		}

		// ���������Ĺ��캯��
		public Cells() {
			
		}

		// �����ݵ��������Ӧ��С����֮��
		public void add(DataPoint p) {
			// ��Ӧ��С���������ݵ�ĸ���+1
			pNum++;

			xall += p.getCoord_X();
			yall += p.getCoord_Y();
			zall += p.getCoord_Z();

			// ���ݵ���ͼƬ�е�����(x,y)=[line,row]
			spacialXall += p.getLine();
			spacialYall += p.getRow();

			// �����ݵ���뵽������
			list.add(p);
		}

		
		// ����С������������(xloc, yloc, zloc)
		public void setCenter(double xloc, double yloc, double zloc) {
			this.xloc = xloc;
			this.yloc = yloc;
			this.zloc = zloc;
		}

		// getter and setter ����
		
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

		//С������������ֵ,ȡ�������ĵ�����ֵ��Ϊ�ö��������ֵ
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

		// ���·�����֪����ʲô��
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
	 * ��Cluster
	 * �������������� һ������������ɸ�С����
	 * 
	 */
	public class Cluster {
		// �������ı��
		private int clusterLabel;

		// ��������x��y��z��������ֵ
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// ��ͬ��������ڵĶ�����ص��z��y��z��������ۼӺ�
		private double xAll = 0;

		private double yAll = 0;

		private double zAll = 0;

		// ���෶Χ���ڰ��������ݵ��ܸ���
		private int pointsAll = 0;

		// ���෶Χ���������꼯
		LinkedList<Cells> contentList;

		// ���༸���������񣨸�С��������ֵ��Ϊ�þ��������ֵ��
		Cells clusterCenter;

		// ��ʼ��������С����
		Cells initialClusterCenter;

		//����ı���
		private double distanceSQ;
		
		// �������������캯��

		// Ѱ�ҳ�ʼ��������С������Χ�����������ڵ�����С����
		public Cluster(int x, int y, int z, int clusterLabel) {
			// ��ʼ�������ĵ�С������cell[x][y][z]
			this.initialClusterCenter = new Cells(x, y, z);

			// �������ı��
			this.clusterLabel = clusterLabel;

			// ��ʼ���෶Χ����С�������꼯
			contentList = new LinkedList<Cells>();

			// ����ʼ��������С������ӵ���ʼ���෶Χ��
			// ��ʼ���෶Χ��ֻ��һ��Ԫ�أ���һ����Ϊ��׼��С����
			contentList.add(initialClusterCenter);

		}// ���캯������

		//�޲������캯��
		public Cluster() {

		}
		
		//��һ�������Ĺ��캯��
		public Cluster(int clusterLabel) {
			this.clusterLabel = clusterLabel;
			contentList = new LinkedList<Cells>();
		}
		
		public Cluster(int x, int y, int z) {
			// ��ʼ�������ĵ�С������cell[x][y][z]
			this.initialClusterCenter = new Cells(x, y, z);

			// ��ʼ���෶Χ����С�������꼯
			contentList = new LinkedList<Cells>();

			// ����ʼ��������С������ӵ���ʼ���෶Χ��
			// ��ʼ���෶Χ��ֻ��һ��Ԫ�أ���һ����Ϊ��׼��С����
			contentList.add(initialClusterCenter);

		}// ���캯������
		
		
		
		// getter and setter����
		
		public Cells getClusterCenter() {
			return clusterCenter;
		}

		public double getDistanceSQ() {
			return distanceSQ;
		}

		//�Ѿ������޸�
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

		// ͳ�ƾ��෶Χ���ڰ���������С���������ݵ��x����ֵ���ۼӺ�
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

		// ͳ�ƾ��෶Χ���ڰ���������С���������ݵ��y����ֵ���ۼӺ�
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

		// ͳ�ƾ��෶Χ���ڰ���������С���������ݵ��z����ֵ���ۼӺ�
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

		// setPointsAll()�������޸ģ�ͳ�ƾ��෶Χ�������ݵ�ĸ���������С����ĸ�����
		public void setPointsAll() {
			for (Cells c : contentList) {
				// ����Ѱ�ҵ����෶Χ���ÿ��С����ͨ��С��������getpNum()�������ÿ��С�����ڵ����ݵ�ĸ���
				pointsAll += cells[c.getX()][c.getY()][c.getZ()].getpNum();
			}
		}		

		//���෶Χ�������е�С����
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

		// ����֮�������෶Χ����ص�
		public void adjustClusterCenter() {
			this.setPointsAll();			
			
			/*
			 * public static double floor(double a)
			 * �������ģ���ӽ��������double ֵ����ֵС�ڵ��ڲ�����������ĳ������
			 * return :�����ӽ�������󣩸���ֵ����ֵС�ڵ��ڸò�����������ĳ������
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
		
		//����ĺ���
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
		
		//������Ϊ����һ��С�����ʱ��,��ΰѸ�С������ӵ����෶Χ����
		public void addCell(Cells c){
			//���С����Ϊ�գ���ӵ����෶Χ��Ҳ����˵����ʹ�ǿյ�С����Ҳ����ĳ�����෶Χ���ڣ�
			if(xloc == -1 && yloc == -1 && zloc ==-1){
				xloc = cells[c.getX()][c.getY()][c.getZ()].getXloc();
				yloc = cells[c.getX()][c.getY()][c.getZ()].getYloc();
				zloc = cells[c.getX()][c.getY()][c.getZ()].getZloc();
			}
			contentList.add(c);
		}
		
		//������Ϊ���С�����ʱ����ΰѸ�С������ӵ����෶Χ����
		public void addCells(LinkedList<Cells> cs){
			for(Cells c : cs){
				contentList.add(c);//һ�ΰ�ÿһ��С������ӵ����෶Χ������ϵͳ��
			}
		}
		
	}

	/**
	 * ��Segment
	 * ����ͼ��ֿ���� ÿһ��������ɸ����ݵ�
	 * 
	 */
	private class Segment {
		// ����
		private int segmentLabel;

		// �ÿ�������ı��
		private int groupLabel = -1;

		// �����������ռ���������
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// �����ռ������ֵͳ�ƣ��ۼӺͣ�
		private double xAll = 0;

		private double yAll = 0;

		private double zAll = 0;;

		// ���а��������ݵ�����
		private double numAll = 0;

		// �������ݵ㼯��
		LinkedList<DataPoint> contentList;

		// ���ڱ߽�㼯��
		LinkedList<DataPoint> marginList;

		// ���ڿ�ı�ż���
		LinkedList<Integer> neighborSegLabelList;

		//getter and setter����
		
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

		// �������Ĺ��캯��
		public Segment(int segmentLabel, int groupLabel) {
			this.segmentLabel = segmentLabel;
			this.groupLabel = groupLabel;
			contentList = new LinkedList<DataPoint>();
			marginList = new LinkedList<DataPoint>();
			neighborSegLabelList = new LinkedList<Integer>();
		}
		// ���������Ĺ��캯��
		public Segment(int segmentLabel) {
			this.segmentLabel = segmentLabel;
			contentList = new LinkedList<DataPoint>();
			marginList = new LinkedList<DataPoint>();
			neighborSegLabelList = new LinkedList<Integer>();
		}

	}

	/**
	 * ��Group
	 * ����ͼ�������,����ĵļ��� ÿһ������ܶ��ͼ���
	 * 
	 */
	private class Group {
		// ����
		private int groupLabel;

		// ס���������ռ���������
		private double xloc = -1;

		private double yloc = -1;

		private double zloc = -1;

		// �����ռ������ֵͳ��
		private double xAll = 0;

		private double yAll = 0;

		private double zAll = 0;

		// ���а��������ݵ�����
		private double numAll = 0;

		// �������ݵ㼯��
		LinkedList<DataPoint> contentList;

		//getter and setter ����
		
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

		// �������Ĺ��캯��
		public Group(int groupLabel) {
			this.groupLabel = groupLabel;
			contentList = new LinkedList<DataPoint>();
		}
		
		
		//���·���ò��û�õ�
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
	 * ����RGBToLUV��
	 * ��������RGBPoints,������ԭʼͼ�����ݵ��������ά����ֵ(R,G,B)����ʽת��Ϊ(L,U,V),���ھ���ͷָ����
	 * Parameters: RGBPoints ����RGB���ݵ� 
	 * return DataPoint[]
	 */
	public DataPoint[] RGBToLUV(DataPoint[] RGBPoints) {

		double L = 0;
		double U = 0;
		double V = 0;

		// ��ά���ݴ��ת������
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
			// ���µõ���LUVֵ�浽ԭ���ص����
			RGBPoints[i].setCoord_X(L);
			RGBPoints[i].setCoord_Y(U);
			RGBPoints[i].setCoord_Z(V);

		}
		return RGBPoints;
	}

	
	/**
	 * ����LUVToRGB����������LUVPoints,���������õ�ͼ�����ݵ��������ά����ֵ(L,U,V)ת��Ϊ(R,G,B),����ͼ�����
	 * Parameters: LUVPoints ����LUV���ݵ� 
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
			// ���µõ���(R,G,B)ֵ�浽ԭLUVPoints[]���ص����
			LUVPoints[i].setCoord_X(R);
			LUVPoints[i].setCoord_Y(G);
			LUVPoints[i].setCoord_Z(B);

		}

		return LUVPoints;
	}

	
	/**
	 * ��Χ������߼�������process 
	 * 1����ʼ�ܶȹ��ƣ� 
	 * 2������Ӧ���ܶȹ��ƣ�
	 * 3����Χ����� ��
	 * 4������Χ�������ӳ�䵽ͼ��ռ䣻
	 */
	public void process() {		
		/**
		 * ��5�� ��Χ�����
		 * 1����������ˮƽnoiseLevel;
		 * 2��ȷ����ʼ��������С����ȷ��ĳ����׼С���񣬱Ƚ�����Χ����С�����ݶȵĲ�ֵ�������ݶȲ�ֵ������һ��������Ϊ��׼���񣬼����ж�
		 * ֱ������������ݶȲ�ֵ�����Ǹ�<0.01����ô�������С�������ĳ����ʼ��������initialClusterCenter��
		 * 3�����ݳ�ʼ�������ģ����������������������򣬲������������÷�Χ���ڵ��������ݵ�����ţ�
		 * 4���жϳ�ʼ��������С���������Ӧ���ܶȺ���adaptive[x][y][z] < noiseLevel�Ļ�����þ��������Ϊ������
		 */

		double constant1 = 0.1;
		noiseLevel = constant1 * geoMean2;
		System.out.println("����ˮƽ noiseLevel = " + noiseLevel);
		
		/**
		 * ���̣�
		 * 1������ȷ��ĳ��С������Ϊ��׼С���� ���⣺���������forѭ������ÿһ��С����Ļ�����ô
		 * 2���ɻ�׼С����Ѱ����һ��nextС���񣻵���searchNextCell(Cells c)��
		 * 3����nextС������Ϊ�µĻ�׼С�������Ѱ����һ��nextС����
		 * 4��ֱ��nextС�����Max < 0.01��Ѱ�ҽ�������ǰС������Ϊ��������centerС����
		 * 5���ҵ���centerΪ���ĵ����е�С����path�����ҵ�������С������Ϊһ����������
		 * 6��Ϊ�����������д���š�
		 */
		
		// �����С�����xƫ����yƫ����zƫ��
		partialDerivative();

		//System.out.println("����ƫ��������");
		
		// �洢��ʼ��������С���񼯺�(ÿһ�������Ӧһ������С����)
		ArrayList<Cells> centers = new ArrayList<Cells>();
		
		//��ʼ���༯��
		ArrayList<Cluster> initialClusters = new ArrayList<Cluster>();
		
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					
					// alreadysearched == 0 ��ʾû�б���������
					if(cells[x][y][z].getAlreadysearched() == 0){						
						
						//�½�һ���������
						Cluster cluster = new Cluster(-1, -1, -1);						
						
						// ��������
						Cells center = new Cells(-1, -1, -1);

						// �洢·���ϵ����е�С�������
						LinkedList<Cells> path = new LinkedList<Cells>();
						//LinkedList<Cells> contentList = new LinkedList<Cells>();
						//System.out.println("��ӡ·������(�������֮ǰ): path.size() = " + path.size());
						
						// ��������
						path.add(new Cells(x, y, z));

						//System.out.println("��ӡ·�����ȣ��������֮��: path.size() = " + path.size());
						System.out.println("cells[" + path.get(0).getX() + "][" + path.get(0).getY() + "][" + path.get(0).getZ() +"]");
						
						// �ɵĻ�׼���ʼ��(�������Ĺ��캯��)
						Cells oldBase = new Cells(-1, -1, -1);
						//oldBase = null;
						
						// �µĻ�׼���ʼ��
						Cells newBase =new Cells(-1, -1, -1);
						//newBase = null;
						
						// ������
						int counter = 0;
						
						//����
						double peak = 1.0;
						
						//��ʼ������
						double[] v = { -1, -1, -1, -1, -1, -1 };	
						do {							
							oldBase = path.get(counter);

							// ��׼С������������ֵ
							int x1 = oldBase.getX();
							int y1 = oldBase.getY();
							int z1 = oldBase.getZ();

							System.out.println("��׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							//��׼��������ȱ��Ϊ1
							cells[x1][y1][z1].setAlreadysearched(1);
							
							v[0] = gapLeftCell(x1, y1, z1);
							v[1] = gapRightCell(x1, y1, z1);
							v[2] = gapUpCell(x1, y1, z1);
							v[3] = gapDownCell(x1, y1, z1);
							v[4] = gapFrontCell(x1, y1, z1);
							v[5] = gapBackCell(x1, y1, z1);
							
							//���δ�ӡ��������ֵ
							for(int i = 0; i < v.length; i++){
								System.out.println("v[" + i +"] = "+ v[i]);
							}
							
							// v[0]���
							if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
									&& v[0] >= v[4] && v[0] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[0] != -1){
									if(v[0] < peak){										
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
									} else {										
										newBase = new Cells(x1 - 1, y1, z1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ (x1 - 1) +"][" + y1 + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1 - 1][y1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}// if����

							
							// v[1]���
							if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
									&& v[1] >= v[4] && v[1] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[1] != -1){
									if(v[1] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//��ΰ�һ�����󸳸���һ������
										//newBase = cells[x1 + 1][y1][z1];
										newBase = new Cells(x1 + 1, y1, z1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ (x1 + 1) +"][" + y1 + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1 + 1][y1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[2]���
							if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
									&& v[2] >= v[4] && v[2] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[2] != -1){
									if(v[2] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1 - 1][z1];
										newBase = new Cells(x1, y1 - 1, z1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + (y1 - 1) + "][" + z1 + "]");
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1][y1 - 1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}							
							

							// v[3]���
							if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
									&& v[3] >= v[4] && v[3] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[3] != -1){
									if(v[3] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1 + 1][z1];
										newBase = new Cells(x1, y1 + 1, z1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + (y1 + 1) + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1][y1 + 1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[4]���
							if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
									&& v[4] >= v[3] && v[4] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[4] != -1){
									if(v[4] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1][z1 - 1];
										newBase = new Cells(x1, y1, z1 - 1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + (z1 - 1) + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1][y1][z1 - 1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[5]���
							if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
									&& v[5] >= v[3] && v[5] >= v[4]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[5] != -1){
									if(v[5] < peak) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1][z1 + 1];
										newBase = new Cells(x1, y1, z1 + 1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + (z1 + 1) + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1][y1][z1 + 1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}	
							}
							
							//�����׼������Χ����С���񶼱��������ˣ�������Ϊһ��
							if(v[0] == -1 && v[1] == -1 && v[2] == -1 && 
									v[3] == -1 && v[4] == -1 && v[4] == -1){
								center = new Cells(x1, y1, z1);
								System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
								newBase = null;
							}

						} while (counter <= path.size() && newBase != null && oldBase !=newBase);

						//��������ͳ�ƹ���
						if (!centers.contains(center)){					
							cluster = new Cluster(center.getX(), center.getY(), center.getZ());
							cluster.setContentList(path);
							System.out.println("�þ��෶Χ����С����ĸ�����" + cluster.getContentList().size());
							centers.add(center);
							initialClusters.add(cluster);
							System.out.println("���ĵ�����Ϊ��" + centers.size());
							System.out.println("��ʼ���༯���ھ��������Ϊ��" + initialClusters.size());							
						}//�ڲ�if����
					}//���if����
					
				}
			}
		}// forѭ������
		
		System.out.println("�ж��ٸ���ʼ������������" + centers.size());
		
		// ��ʼ���༯���о���ĸ���
		initalClusterNum = 1;

		// �洢����������ͼ�����ݵ�
		DataPoint[] toWrite;
		
		//��ʼ���༯���о���ĸ���
		initalClusterNum = initialClusters.size();

		System.out.println("��ʼ���༯���о���ĸ�����" + initalClusterNum);
		
		//Ϊÿһ����ʼ�����趨һ����ţ���Ŵ�1��ʼ
		for(int x = 0, clusterNumber = 1; x < initalClusterNum; x++){
			initialClusters.get(x).setClusterLabel(clusterNumber++);
		}

		// ���ϵ����������վ��༯��
		LinkedList<Cluster> finnalInitList = new LinkedList<Cluster>();
		
		// ���վ��༯���о���ĸ���
		finalClusterNum = 1;
		
		for (int i = 0; i < initalClusterNum; i++) {
			// �����ʼ��������С���������Ӧ���ܶȺ���adaptive[x][y][z] < noiseLevel�Ļ�����þ��������Ϊ����

			// ��ȡĳ������ľ�������С�����x��y��z����
			int x = initialClusters.get(i).initialClusterCenter.getX();
			int y = initialClusters.get(i).initialClusterCenter.getY();
			int z = initialClusters.get(i).initialClusterCenter.getZ();
			double a = adaptive[x][y][z];
			System.out.println("adaptive[x][y][z] = " + a);
			
			if (adaptive[x][y][z] >= noiseLevel) {
				finnalInitList.add(initialClusters.get(i));
			}

		}// forѭ������

		// ���վ���ĸ���
		finalClusterNum = finnalInitList.size();
		System.out.println("���վ��༯���о���ĸ�����" + finalClusterNum);
		
		// allPointsNum��ԭʼͼ�����ݵ�ĸ���
		toWrite = new DataPoint[allPointsNum];
		System.out.println("ԭʼͼ�����ݵ������(toWrite.length)  :" + toWrite.length);		
	
		// ��ʼ��
		for (int i = 0; i < allPointsNum; i++) {
			toWrite[i] = LUVPoints[i];
		}
		
		// ͳ�ƾ����㷨��������ݵ������
		int toWriteS = 0;

		// ͳ�ƾ�����������д�������ݵ������
		int j = 0;

		for (int i = 0; i < finalClusterNum; i++) {
			//System.out.println("��"+ (i + 1) + "�������������: ");
			Cluster cluster = finnalInitList.get(i);// ��ȡ��(i+1)������
						
			// ���෶Χ��(������)�������꼯,Ҳ���������������Χ�����е�С�������ļ���
			LinkedList<Cells> list = cluster.getContentList();
			
			System.out.println("��"+ (i + 1) + "�����෶Χ���ڰ��� " + list.size() +"��С����");
			
			// Ϊ������෶Χ��������С����(����)�������,��Ŵ�1��ʼ
			// coordNumΪС���������
			for (int coordNum = 0; coordNum < list.size(); coordNum++) {
				cells[list.get(coordNum).getX()][list.get(coordNum).getY()][list
						.get(coordNum).getZ()].setClusterNum(i + 1);//����С�������ţ�clusterNum
			}

			int clusterSize = 0;// ���������ݵ����

			// ���������������о����㷨���������ݵ���������(L,U,V)������ClusterLabel
			for (Cells c : list) {    // �������෶Χ�������е�С����				
				Cells grid = cells[c.getX()][c.getY()][c.getZ()];					
				ArrayList<DataPoint> plist = grid.list;// �洢ÿ��С���������ݵ�ļ���				
				clusterSize += plist.size();// ÿ��С���������ݵ�����ۼӺ� = ���������ݵ�ĸ���
				
				// ���ĵ�����
				for (DataPoint p : plist) {
					// ���⣺??
					p.setClusterLabel(i);// ���������ݵ�ľ����ǩ���������ڵľ���ı����ͬ�����ʼ�����Ų�ͬ
					toWrite[(p.getLineNum() - 1)] = p;
					toWriteS++;
				}// �ڲ�foreachѭ������
				
			}// ���foreachѭ������

			j += clusterSize;
			System.out.println(" ��  " + (int) (i + 1) + " ��" + "���� "
					+ clusterSize + " �����ݵ�");
			cluster.setPointsAll(clusterSize);
		}
		
		
		System.out.println("------���� "
				+ j
				+ " �����¼,"
				+ " ƽ��ǰ��ʧ�� "
				+ (LUVPoints.length - j)
				+ " ��,"
				+ " ռ���� "
				+ new DecimalFormat("0.0000")
						.format((double) (LUVPoints.length - j)
								/ LUVPoints.length) + " ------");

		System.out.println(" -----����  " + toWriteS + " ����д��------");

		/**
		 * ��6�� ����Χ�������ӳ�䵽ͼ��ռ� 0�����ú���transformToImagePoints��
		 * �һؾ�������е�������ʧ�㲢���߽�����㰴������֪����ŷʽ������̵�ԭ�����·������ţ�
		 * 1�����ú���smoothpix�Է���������ƽ����ȥ�룬�޸��������ݵ�Ķ�Ӧ���ţ�
		 * 2�����ú���initSegments���������ݵ������,��ʼ������ͼƬ�ֿ�,ȷ��ÿ�����ʼ���,�����ĵ㼯�Ϻ����ڿ�ŵļ��ϣ�
		 * 3�����ú���segmentsArrange����ͳ�Ʒֿ���,Ϊͼ�����ݵ��趨��Ų���������ĵ��������ά���ꣻ
		 * 4�����ú���LUVToRGB��������ع�Ϊ(R,G,B)��
		 * 5����ʼ����ͼAPI�����ú���show_DataPoint_on_Image(BufferedImage image,DataPoint
		 * point,Color color)��ͼ��ʾ������
		 */


		// ȥ�뺯��:1.��ʧ���һ� 2.�߽���������·��� 3.�����ݵ��ʼ�����ͼ������ص�
		//toWrite = transformToimagePoints(LUVPoints, finnalInitList);

		// ��ǰ������ӵ�е��������ݵ�(������ʧ��)
		DataPoint[] imagePoints = new DataPoint[toWrite.length];
		
		// ���ݵ�����
		int clusterLabel = -1;
		
		for(int i = 0; i < toWrite.length; i++){
			// ������δ��ʧ��
			// ������ݵ㲻����ʧ�㣬��ô��ø����ݵ�����
			// ע������ʧ��������-1
			if (toWrite[i].getClusterLabel() != -1) {
				clusterLabel = toWrite[i].getClusterLabel();
			} 
			//�����ݵ��ʼ�����ͼ������ص�
			imagePoints[i] = new DataPoint(toWrite[i].getLine(),
					toWrite[i].getRow(), toWrite[i].getCoord_X(),
					toWrite[i].getCoord_Y(), toWrite[i].getCoord_Z(), 1, 1,
					toWrite[i].getColor(), toWrite[i].getLineNum(), clusterLabel,
					toWrite[i].getSegmentLabel(), toWrite[i].getSmoothed(),
					toWrite[i].getCombined(), toWrite[i].getSegmentChecked());
		}
		System.out.println("imagePoints.length = " + imagePoints.length);
				
		// �Ծ��������ݵ����ƽ��
		toWrite = smoothpix(toWrite, finnalInitList);

		System.out.println("ƽ��������ɣ���");
		// ��ʼ������ͼƬ�ֿ�
		segments = initSegments(toWrite, segments);
		
		DataPoint[] classifiedPoints = new DataPoint[toWrite.length];
		
		// ���ݷ������޸����ݵ�������ռ�ֵ���������ݵ������ֵ�����ļ�����������ֵ��ͬ
		for(Segment s : segments){
			// ָ����ļ�����������ֵ
			double centerL = s.getXloc();
			double centerU = s.getYloc();
			double centerV = s.getZloc();
			// s�ڵ��������ݵ��LUV�����ռ�ֵ����ͬ������s�ļ�����������ֵ
			for (DataPoint point : s.getContentList()) {
				point.setCoord_X(centerL);
				point.setCoord_Y(centerU);
				point.setCoord_Z(centerV);
				classifiedPoints[point.getLineNum() - 1] = point;
			}
		}
		System.out.println("���շֿ�ĸ����� = " + segments.size());
		// ����ֿ���Ϣ,Ϊͼ�����ݵ��趨���
		// ���ط��������ݵ㼯��
		//toWrite = segmentsArrange(segments, toWrite);

		// ��������ع�Ϊ(R,G,B)������GRB����ֵ�����ݵ㼯��
		toWrite = LUVToRGB(toWrite);
		System.out.println("toWrite= " + toWrite.length);
		// �����ͼAPI, ׼����ͼ
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

		g.drawString("��" + toWrite.length + "���㣬�۳�" + finnalInitList.size()
				+ "����", 0, DataBase.sHei - 45);

		g.drawString(time, 0, DataBase.sHei - 5);
		g.dispose();

		g.dispose();
		DataField.updateImagePanel(img);

	}// pocess��������

	
	/**
	 * ����segmentsArrange 
	 * �����е����ݵ㰴�����Ź鲢 
	 * Parameters:
	 * allSegmentsGrouped:���п�����Ľ���� 
	 * allPoints:���е㼯
	 * Return:�ֿ����,�޸�������ά����ֵ����������ݵ㼯��
	 */

	public DataPoint[] segmentsArrange(LinkedList<Segment> allSegmentsGrouped,
			DataPoint[] allPoints) {
		// ���������ݵ㼯��
		DataPoint[] classifiedPoints = new DataPoint[allPoints.length];
		//�����鼯��
		LinkedList<Group> allGroups = new LinkedList<Group>();

		// allSegmentsGrouped.size()�����ڿ�ĸ���
		for (int i = 0; i < allSegmentsGrouped.size(); i++) {
			int included = -1;// ��־λ������ÿ��Ѿ���ĳһ�����У���Ѹ÷�����븳����

			// �жϸÿ��Ƿ������ѽ��������У�ͨ������жϣ�
			for (int j = 0; j < allGroups.size(); j++) {
				// ����ÿ��Ѿ���ĳ�������У����Ӧ�������ͬ
				if (allGroups.get(j).getGroupLabel() == allSegmentsGrouped.get(
						i).getGroupLabel()) {
					included = j;
				}
			}// for����

			/*
			 * ����Ĳ����ǰѿ���뵽��Ӧ�ķ�����
			 */
			// ����ÿ鲻��ĳ�������У������µķ���
			if (included == -1) {
				// ���ʣ����Ѿ�������ˣ���ô��û����Ӧ�ķ����أ���������һ�����Ͻ����·��飬����Ҫ��������ظ��Ĺ���
				Group newGroup = new Group(allSegmentsGrouped.get(i)
						.getGroupLabel());
				// Ϊ�µ��齨�����ݵ㼯��,for{...}�ڴ˴���ָֹ�봫��(ǳ����)
				LinkedList<DataPoint> segentsContentMirrow = new LinkedList<DataPoint>();// ����µķ��������е����ݵ�
				for (int all = 0; all < allSegmentsGrouped.get(i)
						.getContentList().size(); all++) {
					segentsContentMirrow.add(allSegmentsGrouped.get(i)
							.getContentList().get(all));
				}
				newGroup.setContentList(segentsContentMirrow);
				newGroup.setNumAll();
				allGroups.add(newGroup);// ���µķ�����뵽���鼯����
				// ����������з���
			} else {
				Group alreadyExistedGroup = allGroups.get(included);// ��ȡ�ÿ����ڵ���
				LinkedList<DataPoint> newContentList = alreadyExistedGroup
						.getContentList();// ��ȡ�����������ݵ�

				/*
				 * public boolean addAll(Collection<? extends E> c) ���ָ��
				 * collection �е�����Ԫ�ص����б�Ľ�β��˳����ָ�� collection �ĵ�����������ЩԪ�ص�˳��
				 * ���ָ���� collection �ڲ��������б��޸ģ���˲�������Ϊ�ǲ�ȷ���ġ�
				 */
				newContentList.addAll(allSegmentsGrouped.get(i)
						.getContentList());
				alreadyExistedGroup.setContentList(newContentList);// �趨�ÿ�����������ݵ㼯��
				alreadyExistedGroup.setNumAll();
				/*
				 * public E set(int index, E element) �����б���ָ��λ�õ�Ԫ���滻Ϊָ����Ԫ�ء�
				 * index - Ҫ�滻��Ԫ�ص����� element - Ҫ��ָ��λ�ô洢��Ԫ�� return ��ǰ��ָ��λ�õ�Ԫ��
				 */
				allGroups.set(included, alreadyExistedGroup);// ��ָ�����޸�֮�󣬷Ż�ԭ����λ��
			}
		}// ���forѭ������

		// ���ݷ������޸����ݵ�������ռ�ֵ���������ݵ������ֵ�����ļ�����������ֵ��ͬ
		for (Group group : allGroups) {
			// ָ����ļ�����������ֵ
			double centerL = group.getXloc();
			double centerU = group.getYloc();
			double centerV = group.getZloc();

			// group�ڵ��������ݵ��LUV�����ռ�ֵ����ͬ������group�ļ�����������ֵ
			for (DataPoint point : group.getContentList()) {
				point.setCoord_X(centerL);
				point.setCoord_Y(centerU);
				point.setCoord_Z(centerV);
				classifiedPoints[point.getLineNum() - 1] = point;
			}
		}

		// ���շֿ����Ŀ
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
	 * ����outputGraphStructuretoFile ������ͼ״�ṹ�����ı�,�˴�����Ϊ���п飬������������������֮�����һ����
	 */
	public static void outputGraphStructuretoFile(LinkedList<Segment> segments) {
		String singleLine = new String();
		try {

			File f = new File("D:\\Graduation Project\\src\\0331\\"
					+ "ClusterForTD\\GraphStructure.txt");

			if (f.exists()) {
				System.out.println("GraphStructure���ڣ�ɾ��");
				f.delete();
			} else {
				System.out.println("GraphStructure�����ڣ����ڴ���...");
				if (f.createNewFile()) {
					System.out.println("GraphStructure�����ɹ���");
				} else {
					System.out.println("GraphStructure����ʧ�ܣ�");
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
	 * ����initSegments 
	 * �������ݵ�ľ�������ʼ�����зֿ�,��������ı��,���ݵ㼯�Ϻ��ڿ��ż��� 
	 * Parameters:
	 * points:����ƽ������������� 
	 * segments:��ʼ���б� 
	 * Return:��ʼ�ֿ�����
	 */
	public LinkedList<Segment> initSegments(DataPoint[] points,
			LinkedList<Segment> segments) {
		// ���ţ���������
		int segmentIndex = 0;

		for (int i = 0; i < points.length; i++) {
			// �����С
			int sizeofOnePiece = 0;
			int k = 0;
			
			if (points[i].getSegmentChecked() == 0) {//����Ƿ��Ѿ����ֿ�
				boolean tag;
				// ͼ���е�����
				Segment oneSeg = new Segment(segmentIndex);

				// ����õ�����
				oneSeg.getContentList().add(points[i]);// �������ݵ�ļ���
				points[i].setSegmentChecked(1);// ��ǵ�i�����ݵ��Ѿ��������
				points[i].setSegmentLabel(segmentIndex);// ��i�����ݵ������

				// �����ļ��߽�ĳ�ʼ����
				// ���������ݵ㼯��
				LinkedList<DataPoint> mirrorContentList = oneSeg.getContentList();
				System.out.println("��ʼ���������ݵ㼯�������ݵ������" + mirrorContentList.size());
				// �߽����ݵ㼯��
				LinkedList<DataPoint> mirrorMarginList = oneSeg.getMarginList();
				do {
					tag = true;
					// ����������ж���temp��¼�����ݵ��ȫ���������ݵ�,�������Ŀ�
					//LinkedList<DataPoint> temp = allFourNeighborFeilds(mirrorContentList.get(k));//���ص�k�����ݵ�
					LinkedList<DataPoint> temp = allEightNeighborFeilds(mirrorContentList.get(k));//���ص�k�����ݵ�
					System.out.println("temp.size() =" + temp.size());
					
					for (DataPoint p : temp) {
						//System.out.println("p.getClusterLabel() =" + p.getClusterLabel());
						//System.out.println("points[i].getClusterLabel() =" + points[i].getClusterLabel());
						// ���뱻���������ͬһ����(��ͬһ��)�����ݵ��������
						if (p.getClusterLabel() == points[i].getClusterLabel() ) {
							// ֻ����Щû�б��ֿ�����ݵ�
							if (p.getSegmentChecked() == 0) {
								mirrorContentList.add(p);
								p.setSegmentChecked(1);
								p.setSegmentLabel(segmentIndex);
							}
						}
						// ��ĳ�����������(�������)���е�������������߿�,��õ�Ϊ������ı߽��
						else {
							if (!mirrorMarginList.contains(mirrorContentList.get(k)))
								// ����Щ��׼�������Χ�����׼�㲻ͬ����ĵ㣬���뵽�߽�㼯��
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
				//System.out.println("�ֿ�֮�����������ݵ㼯�������ݵ������" + mirrorContentList.size());
				oneSeg.setMarginList(mirrorMarginList);
				//System.out.println("�ֿ�֮��߽����ݵ㼯�������ݵ������" + mirrorMarginList.size());
				oneSeg.setNumAll();
				segments.add(oneSeg);
				//System.out.println("�ֿ������" + segments.size());
				segmentIndex++;//���Ų�������
			}//if����

		}//forѭ������

		System.out.println("�ֿ������" + segments.size());
		
		// �ڶ���ɨ�裬��ʼ��ÿ��������ڿ���Ϣ
		for (int j = 0; j < segments.size(); j++) {
			// �߽缰����ĳ�ʼ����
			// ���ڿ�ı�ż���
			LinkedList<Integer> mirrorNeighborSegLabelList = segments.get(j)
					.getNeighborSegLabelList();
			
			// ���������ݿ���������ݵ㼯��
			LinkedList<DataPoint> mirrorMarginList = segments.get(j)
					.getMarginList();

			//���������������ݵ�
			for (DataPoint marginPoint : mirrorMarginList) {
				//ͳ���������ݵ����������ݵ㼯��temp2
				//LinkedList<DataPoint> temp2 = allFourNeighborFeilds(marginPoint);
				LinkedList<DataPoint> temp2 = allEightNeighborFeilds(marginPoint);
				//�������������ݵ㼯��
				for (DataPoint neighboorPoint : temp2) {
					// ���뱻����㲻ͬ������ݿ�ż�¼����
					if (neighboorPoint.getSegmentLabel() != marginPoint
							.getSegmentLabel()) {
						// ��ȡ�뱻������Ų�ͬ����Щ���
						Integer neighboorSegLab = new Integer(
								neighboorPoint.getSegmentLabel());
						if (!mirrorNeighborSegLabelList
								.contains(neighboorSegLab)
								|| mirrorNeighborSegLabelList.size() == 0) {
							// �뱻����㲻ͬ������ݿ�ŵļ���
							mirrorNeighborSegLabelList.add(neighboorSegLab);							
						}
					}
				}
			}
			segments.get(j).setNeighborSegLabelList(mirrorNeighborSegLabelList);
			System.out.println("���" + j + "�����ڿ�ı���У�"+ 
					mirrorNeighborSegLabelList.size() + "��");
		}
		return segments;
	}

	
	/**
	 * ����transformToImagePoints 
	 * �����ݵ��ʼ�����ͼ������ص� 
	 * Parameters:
	 * points:��ǰͼ��LUV����ϵ���е����ݵ� 
	 * clusters�����༯�� 
	 * cells��С������� 
	 * return ���ͼ������ص㼯��
	 */
	public static DataPoint[] transformToimagePoints(DataPoint[] points,
			LinkedList<Cluster> clusters) {

		// ��ǰ������ӵ�е��������ݵ�(������ʧ��)
		DataPoint[] imagePoints = new DataPoint[points.length];

		// ���ݵ�����
		int clusterLabel;

		for (int i = 0; i < points.length; i++) {
			
			// ������δ��ʧ��
			// ������ݵ㲻����ʧ�㣬��ô��ø����ݵ�����(��ʧ��û����ţ�clusterLabel )
			// ע������ʧ��������-1
			if (points[i].getClusterLabel() != -1) {
				clusterLabel = points[i].getClusterLabel();
			} 
			
			// ��������ʧ��
			// ������ݵ㱾�������ʧ�㣬��ô��������ݵ㵽ÿ�����༸����������ŷʽ���ξ��룬�һ���ʧ��
			else { 

				// ÿ��������������ֵ(clusterCenterX, clusterCenterY, clusterCenterZ)
				double clusterCenterX = 0;
				double clusterCenterY = 0;
				double clusterCenterZ = 0;

				// ��ʧ����������ļ�����ŷʽ����
				double minDistance = Double.MAX_VALUE;
				//System.out.println("minDistance = " +minDistance);
				
				// �鲢���,��������ʧ����������ļ��ŷʽ������̵ľ�������
				int minDistanceCluster = 0;

				//�������еľ���
				for (int j = 0; j < clusters.size(); j++) {
					double dis = 0;
					
					// �þ��༸������С���񣨲���ɽ��ⶥ��С����Ҳ����˵���Ǿ������ģ�
					Cluster cluster = clusters.get(j);
					System.out.println("��" + (j + 1) + "�����෶Χ�����У� " +cluster.contentList.size()+ "��С����" +
							"�У� " +cluster.getPointsAll()+ " �����ݵ�" );
					//ȷ���˸þ���ļ�������С����
					cluster.adjustClusterCenter();
					//��ȡ��������С����
					Cells geoCenter = cluster.getClusterCenter();

					// ��ȡ���༸����������ֵ�����������ֵ��Ϊ���������ֵ
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

					// �������(ŷʽ����)
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
						minDistanceCluster = j;// j�Ǿ������
					}

				}// forѭ������
				//���������Ѱ�Ҿ�����̵ľ�������
				clusterLabel = minDistanceCluster;
			}// else����
			//���������Ϊͼ����û��һ�����ݵ����һ�����
			System.out.println("��" + i + "�����ݵ������ǣ�" +clusterLabel + "�� Ϊ��"+ i +"�����ݵ���������ɣ���");
			
			//�����ݵ��ʼ�����ͼ������ص�
			imagePoints[i] = new DataPoint(points[i].getLine(),
					points[i].getRow(), points[i].getCoord_X(),
					points[i].getCoord_Y(), points[i].getCoord_Z(), 1, 1,
					points[i].getColor(), points[i].getLineNum(), clusterLabel,
					points[i].getSegmentLabel(), points[i].getSmoothed(),
					points[i].getCombined(), points[i].getSegmentChecked());
		}// forѭ������������˶�ÿ�����ݵ�Ĳ���

		return imagePoints;
	}

	/**
	 * ����smoothpix 
	 * ������ݾ����ƽ��ȥ�� 
	 * Parameters: 
	 * points: �����������ƽ�������ݵ㼯 
	 * clusters: ������ľ����б� 
	 * Return ƽ��������ݼ���
	 */
	public DataPoint[] smoothpix(DataPoint[] points,
			LinkedList<Cluster> clusters) {
		smoothValve = 40;
		// ��ʼ��allPix,ԭʼͼ�������ֵ
		for (int i = 0; i < points.length; i++) {
			allPix[(int) points[i].getLine()][(int) points[i].getRow()] = points[i];
		}

		for (int i = 0; i < points.length; i++) {

			// ��¼�������������Χ(�߽�)���������ݵ�ľ�����Ϣ
			int[] neighboorCluster = new int[finalClusterNum + 1];

			// �������������������ݵ�ļ���
			LinkedList<DataPoint> neighboor = new LinkedList<DataPoint>();

			// ������������Ĵ�С(�������������������ݵ�ĸ���)
			int NBsize = 0;

			// �����������������Χ(�߽�)���ݵ�����������һ��ı�ţ���Ҫƽ�������ݵ㽫��ƽ��������(��׼��)
			int NBMax = 0;

			int k = 0;// �������

			// ��ʼ���������������Χ���������ݵ�ľ�����Ϣ����ʼ���������,neighboorCluster.length��ʾ�������
			for (int j = 0; j < neighboorCluster.length; j++) {
				neighboorCluster[j] = 0;// ��j����������ʼ��Ϊ0
			}

			if (points[i].getSmoothed() == 0) {// ��������ݵ�û�б�ƽ��
				boolean tag;
				// ����õ�����
				neighboor.add(points[i]);
				do {
					tag = true;
					// ���ҵ�����ı��Ϊk�����ݵ㣬����allFourNeighborFeilds����������ָ�����ݵ����������ݵ㼯��
					//LinkedList<DataPoint> temp = allFourNeighborFeilds(neighboor.get(k));// ���� �������ж���temp��¼�����ݵ��ȫ���������ݵ�
					LinkedList<DataPoint> temp = allEightNeighborFeilds(neighboor.get(k));// ���� �������ж���temp��¼�����ݵ��ȫ���������ݵ�

					//���α�����������ÿһ�����ݵ� 
					for (DataPoint p : temp) {
						// ���뱻���������ͬһ��������ݵ��������(���ݾ����ǩ�ж�)
						if (p.getClusterLabel() == points[i].getClusterLabel()) {
							if (!neighboor.contains(p)) {// ���ָ�����ݵ�����򲻰������������ݵ㼯���еĵ㣬��ô���������뵽ָ���������
								neighboor.add(p);
							}
						}
						// ������������ݵ������İ��������ݵ㲻����ͬһ������Ļ��������������߽磬����λ��׼��
						else {
							int clusterLabel = p.getClusterLabel();
							if (clusterLabel != -1) {
								clusterLabel++;
							} else {
								clusterLabel = 0;
							}
							neighboorCluster[clusterLabel]++;
							// �����������������Χ(�߽�)���ݵ�����������һ��ı�ţ���Ҫƽ�������ݵ㽫��ƽ��������(��׼��)
							if (NBMax == -1) {
								NBMax++;
							}
							if (neighboorCluster[clusterLabel] > neighboorCluster[NBMax]) {
								clusterLabel--;
								NBMax = clusterLabel;
							}
						}// else����
					}// foreachѭ������
					k++;
					// �������������������ݵ�ĸ���
					NBsize = neighboor.size();
					System.out.println("�������������������ݵ�ĸ����� " + NBsize);
					if (k == NBsize) {
						// �����Ѿ�������ϲ�����������ѭ��
						tag = false;
					}
					// �����С����ƽ������40������ֹͣ������ѭ����ֹ
				} while (NBsize < smoothValve && tag == true);

				/*
				 * ���뷽��:������ٽ���鲢 ����Ҫƽ���ͳ�������ݵ㼰���������е�ƽ������׼�� ��Ҫƽ���ͳ���ĵ� 
				 * 1��NBsize < 40 ��ʾ�������ƽ�����޵������㣻 
				 * 2��neighboor.get(0).getClusterLabel() == -1��ʾ�������о�Ϊ�����㷨��ʧ��ͼ�����ݵ㣻
				 */
				if (NBsize < smoothValve
						|| neighboor.get(0).getClusterLabel() == -1) {
					// �����׼�����Ų�����-1
					if (NBMax != -1) {
						for (DataPoint p : neighboor) {
							// ���������������������ݵ���������Ϊ�������������Χ(�߽�)���ݵ�����������һ��ı��
							p.setClusterLabel(NBMax);
						}
					}
					// �����������ʧ����������Ϊȱʡֵ
					else {
						for (DataPoint p : neighboor) {
							p.setClusterLabel(NBMax);
						}
					}
				}// ���if����

				// ���ڲ���Ҫƽ���ĵ㣬��������Ļ�ͨ�ԣ�����Ҫ�ٴβ����´��ж�
				else {
					for (DataPoint p : neighboor) {
						p.setSmoothed(1);// ����Ҫƽ�������ݵ㣬����Ϊ����ƽ��
					}
				}
			}
		}// �����forѭ������

		// ԭͼ�����ݾ������
		for (int i = 0; i < points.length; i++) {
			allPix[(int) points[i].getLine()][(int) points[i].getRow()] = points[i];
		}

		return points;
	}

	/**
	 * ����allEightNeighborFeilds ����ָ�����ݵ���ԭͼ���еİ����� Parameters: point:ָ�������ݵ�
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
	 * ����allFourNeighborFeilds 
	 * ����ָ�����ݵ���ԭͼ���е������� 
	 * Parameters: point:ָ�������ݵ�
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
	 * ����initialClusterCenter 
	 * ������ʼ�������� 
	 * Parameters: 
	 * c ָ����ĳ��С������� 
	 * return ��ʼ��������С����center����
	 * author weil
	 */
	public Cells initialClusterCenter(Cells c) {
		// �����С�����xƫ����yƫ����zƫ��
		partialDerivative();

		// ��������
		Cells center = new Cells();

		// �洢·���ϵ����е�С�������
		LinkedList<Cells> path = new LinkedList<Cells>();

		// ��������
		path.add(c);

		// �ɵĻ�׼��
		Cells oldBase = new Cells();

		// �µĻ�׼��
		Cells newBase = new Cells();

		// ������
		int counter = 0;

		do {
			oldBase = path.get(counter);

			// ��׼С������������ֵ
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
			
			// v[0]���
			if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
					&& v[0] >= v[4] && v[0] >= v[5]) {
				/*
				 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
				 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
				 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
				 * == null)ʱ������
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
						// ��ɽ�������Ѿ���������С������Ϊ����������
						cells[x1 -1][y1][z1].setAlreadysearched(1);
						counter++;
					}
				}

			}// if����

			// v[1]���
			if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
					&& v[1] >= v[4] && v[1] >= v[5]) {
				/*
				 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
				 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
				 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
				 * == null)ʱ������
				 */
				if (v[1] < 0.5) {
					//center = oldBase;
					center = new Cells(x1, y1, z1);
					newBase = null;
				} else {
					//��ΰ�һ�����󸳸���һ������
					//newBase = cells[x1 + 1][y1][z1];
					newBase = new Cells(x1 + 1, y1, z1);
					if (!path.contains(newBase)) {
						path.add(newBase);
						// ��ɽ�������Ѿ���������С������Ϊ����������
						cells[x1 + 1][y1][z1].setAlreadysearched(1);
						counter++;
					}
				}

			}

			// v[2]���
			if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
					&& v[2] >= v[4] && v[2] >= v[5]) {
				/*
				 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
				 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
				 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
				 * == null)ʱ������
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
						// ��ɽ�������Ѿ���������С������Ϊ����������
						cells[x1][y1 - 1][z1].setAlreadysearched(1);
						counter++;
					}
				}

			}

			// v[3]���
			if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
					&& v[3] >= v[4] && v[3] >= v[5]) {
				/*
				 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
				 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
				 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
				 * == null)ʱ������
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
						// ��ɽ�������Ѿ���������С������Ϊ����������
						cells[x1][y1 + 1][z1].setAlreadysearched(1);
						counter++;
					}
				}

			}

			// v[4]���
			if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
					&& v[4] >= v[3] && v[4] >= v[5]) {
				/*
				 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
				 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
				 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
				 * == null)ʱ������
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
						// ��ɽ�������Ѿ���������С������Ϊ����������
						cells[x1][y1][z1 - 1].setAlreadysearched(1);
						counter++;
					}
				}

			}

			// v[5]���
			if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
					&& v[5] >= v[3] && v[5] >= v[4]) {
				/*
				 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
				 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
				 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
				 * == null)ʱ������
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
						// ��ɽ�������Ѿ���������С������Ϊ����������
						cells[x1][y1][z1 + 1].setAlreadysearched(1);
						counter++;
					}
				}

			}

		} while (counter <= path.size() && newBase != null);
		
		return center;
	}
	
	

	/**
	 * ����initialCluster 
	 * ������ʼ�������� 
	 * Parameters: null 
	 * return ���г�ʼ��������С����center���󼯺�
	 * author weil
	 */
	public ArrayList<Cells> initialCluster() {
		// �����С�����xƫ����yƫ����zƫ��
		partialDerivative();

		System.out.println("����ƫ��������");
		
		// �洢��ʼ��������С���񼯺�
		ArrayList<Cells> centers = new ArrayList<Cells>();
		
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					// ��������
					Cells center = new Cells();

					// �洢·���ϵ����е�С�������
					LinkedList<Cells> path = new LinkedList<Cells>();
					
					System.out.println("��ӡ·������(�������֮ǰ): path.size() = " + path.size());
					
					// ��������
					path.add(new Cells(x, y, z));

					System.out.println("��ӡ·�����ȣ��������֮��: path.size() = " + path.size());
					System.out.println("cells[" + path.get(0).getX() + "][" + path.get(0).getY() + "][" + path.get(0).getZ() +"]");
					
					// �ɵĻ�׼���ʼ��
					Cells oldBase = new Cells();
					oldBase = null;
					
					// �µĻ�׼���ʼ��
					Cells newBase =new Cells();
					newBase = null;
					
					// ������
					int counter = 0;

					//alreadysearched == 0 ��ʾû�б���������
					if(cells[x][y][z].getAlreadysearched() == 0){
						double[] v = { -1, -1, -1, -1, -1, -1 };
						
						//LinkedList<Cells> contentList = new LinkedList<Cells>();
						
						do {
							
							oldBase = path.get(counter);

							// ��׼С������������ֵ
							int x1 = oldBase.getX();
							int y1 = oldBase.getY();
							int z1 = oldBase.getZ();

							System.out.println("��׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							cells[x1][y1][z1].setAlreadysearched(1);
							
							v[0] = gapLeftCell(x1, y1, z1);
							v[1] = gapRightCell(x1, y1, z1);
							v[2] = gapUpCell(x1, y1, z1);
							v[3] = gapDownCell(x1, y1, z1);
							v[4] = gapFrontCell(x1, y1, z1);
							v[5] = gapBackCell(x1, y1, z1);
							//���δ�ӡ��������ֵ
							for(int i = 0; i < v.length; i++){
								System.out.println("v[" + i +"] = "+ v[i]);
							}
							
							// v[0]���
							if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
									&& v[0] >= v[4] && v[0] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[0] != -1){
									if(v[0] < 1.0){
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
									} else {
										//newBase = cells[x1 - 1][y1][z1];
										newBase = new Cells(x1 - 1, y1, z1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ (x1 - 1) +"][" + y1 + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1 -1][y1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}// if����

							
							// v[1]���
							if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
									&& v[1] >= v[4] && v[1] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[1] != -1){
									if(v[1] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//��ΰ�һ�����󸳸���һ������
										//newBase = cells[x1 + 1][y1][z1];
										newBase = new Cells(x1 + 1, y1, z1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ (x1 + 1) +"][" + y1 + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1 + 1][y1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[2]���
							if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
									&& v[2] >= v[4] && v[2] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[2] != -1){
									if(v[2] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1 - 1][z1];
										newBase = new Cells(x1, y1 - 1, z1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + (y1 - 1) + "][" + z1 + "]");
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1][y1 - 1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}							
							

							// v[3]���
							if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
									&& v[3] >= v[4] && v[3] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[3] != -1){
									if(v[3] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1 + 1][z1];
										newBase = new Cells(x1, y1 + 1, z1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + (y1 + 1) + "][" + z1 + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1][y1 + 1][z1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[4]���
							if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
									&& v[4] >= v[3] && v[4] >= v[5]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[4] != -1){
									if(v[4] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1][z1 - 1];
										newBase = new Cells(x1, y1, z1 - 1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + (z1 - 1) + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1][y1][z1 - 1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}
							}
							

							// v[5]���
							if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
									&& v[5] >= v[3] && v[5] >= v[4]) {
								/*
								 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
								 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
								 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
								 * == null)ʱ������
								 */
								if(v[5] != -1){
									if(v[5] < 1.0) {
										//center = oldBase;
										center = new Cells(x1, y1, z1);
										System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
										newBase = null;
										
									} else {
										//newBase = cells[x1][y1][z1 + 1];
										newBase = new Cells(x1, y1, z1 + 1);
										System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + (z1 + 1) + "]");
										
										if (!path.contains(newBase)) {
											path.add(newBase);
											// ��ɽ�������Ѿ���������С������Ϊ����������
											cells[x1][y1][z1 + 1].setAlreadysearched(1);
											counter++;
											System.out.println("counter = " + counter);
											System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
										}
									}
								}	
							}

						} while (counter <= path.size() && newBase != null && oldBase !=newBase);

						//��������ͳ�ƹ���
						if (!centers.contains(center)) {					

							centers.add(center);
							System.out.println("���ĵ�����Ϊ��" + centers.size());
						}//�ڲ�if����
					}//���if����
					
				}
			}
		}// forѭ������

		return centers;
	}

	
	/**
	 * ����searchPathCells 
	 * ÿһ������·��������С������� 
	 * Parameters: c ָ��·���ϵ�ĳ��С�������
	 * return ָ��С�����������·���ϵ����е�С������� 
	 * author weil
	 */
	public LinkedList<Cells> searchPathCells(Cells c) {

		// �洢·���ϵ����е�С�������
		LinkedList<Cells> path = new LinkedList<Cells>();

		// ��������
		path.add(c);

		// ·���ϵľ�������С�������
		Cells center = new Cells();
		
		// �ɵĻ�׼���ʼ��
		Cells oldBase = new Cells();

		// �µĻ�׼���ʼ��
		Cells newBase =new Cells();

		// ������
		int counter = 0;

		//alreadysearched == 0 ��ʾû�б���������
		if(cells[c.getX()][c.getY()][c.getZ()].getAlreadysearched() == 0){
			double[] v = { -1, -1, -1, -1, -1, -1 };
			
			do {
				oldBase = path.get(counter);

				// ��׼С������������ֵ
				int x1 = oldBase.getX();
				int y1 = oldBase.getY();
				int z1 = oldBase.getZ();

				System.out.println("��׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
				cells[x1][y1][z1].setAlreadysearched(1);
				
				v[0] = gapLeftCell(x1, y1, z1);
				v[1] = gapRightCell(x1, y1, z1);
				v[2] = gapUpCell(x1, y1, z1);
				v[3] = gapDownCell(x1, y1, z1);
				v[4] = gapFrontCell(x1, y1, z1);
				v[5] = gapBackCell(x1, y1, z1);
				//���δ�ӡ��������ֵ
				for(int i = 0; i < v.length; i++){
					System.out.println("v[" + i +"] = "+ v[i]);
				}
				
				// v[0]���
				if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
						&& v[0] >= v[4] && v[0] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[0] != -1){
						if(v[0] < 1.0){
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1 - 1][y1][z1];
							newBase = new Cells(x1 - 1, y1, z1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ (x1 - 1) +"][" + y1 + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1 -1][y1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}// if����

				
				// v[1]���
				if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
						&& v[1] >= v[4] && v[1] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[1] != -1){
						if(v[1] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//��ΰ�һ�����󸳸���һ������
							//newBase = cells[x1 + 1][y1][z1];
							newBase = new Cells(x1 + 1, y1, z1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ (x1 + 1) +"][" + y1 + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1 + 1][y1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[2]���
				if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
						&& v[2] >= v[4] && v[2] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[2] != -1){
						if(v[2] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1 - 1][z1];
							newBase = new Cells(x1, y1 - 1, z1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + (y1 - 1) + "][" + z1 + "]");
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1][y1 - 1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}							
				

				// v[3]���
				if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
						&& v[3] >= v[4] && v[3] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[3] != -1){
						if(v[3] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1 + 1][z1];
							newBase = new Cells(x1, y1 + 1, z1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + (y1 + 1) + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1][y1 + 1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[4]���
				if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
						&& v[4] >= v[3] && v[4] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[4] != -1){
						if(v[4] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1][z1 - 1];
							newBase = new Cells(x1, y1, z1 - 1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + (z1 - 1) + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1][y1][z1 - 1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[5]���
				if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
						&& v[5] >= v[3] && v[5] >= v[4]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[5] != -1){
						if(v[5] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1][z1 + 1];
							newBase = new Cells(x1, y1, z1 + 1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + (z1 + 1) + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1][y1][z1 + 1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}	
				}

			} while (counter <= path.size() && newBase != null && oldBase !=newBase);
		}

		return path;
	}
	
	
	/**
	 * ����searchforCluster
	 * ����ָ�����෶Χ���ڵ�����С����
	 * Parameters: 
	 * cluster��ָ����������
	 * return : ���෶Χ��������С����
	 * 
	 */
	
	public LinkedList<Cells> searchforCluster(Cluster c){

		//��ȡ�þ���ľ�������
		Cells center = c.initialClusterCenter;
		System.out.println("x = " + center.getX() +", y = " + center.getY() + ", z = " + center.getZ());
		
		// �洢·���ϵ����е�С�������
		LinkedList<Cells> path = new LinkedList<Cells>();
		
		System.out.println("��ӡ·������(�������֮ǰ): path.size() = " + path.size());
		
		// ��������
		path.add(new Cells(center.getX(), center.getY(), center.getZ()));

		System.out.println("��ӡ·�����ȣ��������֮��: path.size() = " + path.size());
		System.out.println("cells[" + path.get(0).getX() + "][" + path.get(0).getY() + "][" + path.get(0).getZ() +"]");
		
		// �ɵĻ�׼���ʼ��
		Cells oldBase = new Cells();
		oldBase = null;
		
		// �µĻ�׼���ʼ��
		Cells newBase =new Cells();
		newBase = null;
		
		// ������
		int counter = 0;

		//alreadysearched == 0 ��ʾû�б���������
		if(cells[center.getX()][center.getY()][center.getZ()].getAlreadysearched() == 0){
			double[] v = { -1, -1, -1, -1, -1, -1 };
				
			do {
				oldBase = path.get(counter);

				// ��׼С������������ֵ
				int x1 = oldBase.getX();
				int y1 = oldBase.getY();
				int z1 = oldBase.getZ();

				System.out.println("��׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
				cells[x1][y1][z1].setAlreadysearched(1);
				
				v[0] = gapLeftCell(x1, y1, z1);
				v[1] = gapRightCell(x1, y1, z1);
				v[2] = gapUpCell(x1, y1, z1);
				v[3] = gapDownCell(x1, y1, z1);
				v[4] = gapFrontCell(x1, y1, z1);
				v[5] = gapBackCell(x1, y1, z1);
				//���δ�ӡ��������ֵ
				for(int i = 0; i < v.length; i++){
					System.out.println("v[" + i +"] = "+ v[i]);
				}
				
				// v[0]���
				if (v[0] >= v[1] && v[0] >= v[2] && v[0] >= v[3]
						&& v[0] >= v[4] && v[0] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[0] != -1){
						if(v[0] < 1.0){
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
						} else {
							//newBase = cells[x1 - 1][y1][z1];
							newBase = new Cells(x1 - 1, y1, z1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ (x1 - 1) +"][" + y1 + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1 -1][y1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}// if����

				
				// v[1]���
				if (v[1] >= v[0] && v[1] >= v[2] && v[1] >= v[3]
						&& v[1] >= v[4] && v[1] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[1] != -1){
						if(v[1] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//��ΰ�һ�����󸳸���һ������
							//newBase = cells[x1 + 1][y1][z1];
							newBase = new Cells(x1 + 1, y1, z1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ (x1 + 1) +"][" + y1 + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1 + 1][y1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[2]���
				if (v[2] >= v[0] && v[2] >= v[1] && v[2] >= v[3]
						&& v[2] >= v[4] && v[2] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[2] != -1){
						if(v[2] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1 - 1][z1];
							newBase = new Cells(x1, y1 - 1, z1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + (y1 - 1) + "][" + z1 + "]");
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1][y1 - 1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}							
				

				// v[3]���
				if (v[3] >= v[0] && v[3] >= v[1] && v[3] >= v[2]
						&& v[3] >= v[4] && v[3] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[3] != -1){
						if(v[3] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1 + 1][z1];
							newBase = new Cells(x1, y1 + 1, z1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + (y1 + 1) + "][" + z1 + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1][y1 + 1][z1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[4]���
				if (v[4] >= v[0] && v[4] >= v[1] && v[4] >= v[2]
						&& v[4] >= v[3] && v[4] >= v[5]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[4] != -1){
						if(v[4] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1][z1 - 1];
							newBase = new Cells(x1, y1, z1 - 1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + (z1 - 1) + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1][y1][z1 - 1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}
				}
				

				// v[5]���
				if (v[5] >= v[0] && v[5] >= v[1] && v[5] >= v[2]
						&& v[5] >= v[3] && v[5] >= v[4]) {
					/*
					 * ����ע������Ѱ�ҵ���������С������ݶȲ�ֵ < 0.01
					 * ʱ����Ϊ��oldBase��·���Ķ��壨����������С����
					 * ��ô��ɽ���̲�Ҫ������ȥ�ˣ�Ҳ����˵newBase = null;�����ڼ��if(oldBase
					 * == null)ʱ������
					 */
					if(v[5] != -1){
						if(v[5] < 1.0) {
							//center = oldBase;
							center = new Cells(x1, y1, z1);
							System.out.println("������ɽ·���Ķ���Ϊ:cells["+ x1 +"][" + y1 + "][" + z1 + "]");
							newBase = null;
							
						} else {
							//newBase = cells[x1][y1][z1 + 1];
							newBase = new Cells(x1, y1, z1 + 1);
							System.out.println("�µĻ�׼��Ϊ:cells["+ x1 +"][" + y1 + "][" + (z1 + 1) + "]");
							
							if (!path.contains(newBase)) {
								path.add(newBase);
								// ��ɽ�������Ѿ���������С������Ϊ����������
								cells[x1][y1][z1 + 1].setAlreadysearched(1);
								counter++;
								System.out.println("counter = " + counter);
								System.out.println("����·���ϵ�С�������:path.size() =  " + path.size());
							}
						}
					}	
				}

			} while (counter <= path.size() && newBase != null && oldBase !=newBase);
		}
					
		return path;		
	}
	
	// ��������������Ǽ����������ݶȲ�ֵ��
	/**
	 * ����gapLeftCell 
	 * ���㵱ǰС�����������С����Ĳ�ֵ 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return ��ǰС�����������С����Ĳ�ֵ 
	 * author weil
	 */
	public double gapLeftCell(int x, int y, int z) {

		if (x == 0)// ����ߵ�����û����ߵ������ˣ�
			return -1;

		// ��ߵ�С�����Ѿ����Ϊ�������������򷵻�-1
		if (cells[x - 1][y][z].getAlreadysearched() != 0) {
			return -1;
		}

		// ��������С������ݶȲ�ֵ
		double xd = 0;// x����ƫ������ֵ��ƽ��

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
	 * ����gapRightCell 
	 * ���㵱ǰС���������ұ�С����Ĳ�ֵ 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y����
	 * z ��ǰС�����z���� 
	 * return ��ǰС���������ұ�С����Ĳ�ֵ 
	 * author weil
	 */
	public double gapRightCell(int x, int y, int z) {

		if (x == K - 1)// ���ұߵ�����û���ұߵ������ˣ�
			return -1;

		// �ұߵ�С�����Ѿ����Ϊ�������������򷵻�-1
		if (cells[x + 1][y][z].getAlreadysearched() != 0) {
			return -1;
		}

		// ��������С������ݶȲ�ֵ
		double xd = 0;// x����ƫ������ֵ��ƽ��

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
	 * ����gapUpCell 
	 * ���㵱ǰС���������ϱ�С����Ĳ�ֵ 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z����
	 * return ��ǰС���������ϱ�С����Ĳ�ֵ 
	 * author weil
	 */
	public double gapUpCell(int x, int y, int z) {

		if (y == 0)// ���ϱߵ�����û���ϱߵ������ˣ�
			return -1;

		// �ϱߵ�С�����Ѿ����Ϊ�������������򷵻�-1
		if (cells[x][y - 1][z].getAlreadysearched() != 0) {
			return -1;
		}

		// ��������С������ݶȲ�ֵ
		double xd = 0;// x����ƫ������ֵ��ƽ��

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
	 * ����gapDownCell 
	 * ���㵱ǰС���������±�С����Ĳ�ֵ 
	 * Parameters: x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return ��ǰС���������±�С����Ĳ�ֵ 
	 * author weil
	 */
	public double gapDownCell(int x, int y, int z) {

		if (y == K - 1)// ��ǰ�ߵ�����û��ǰ�ߵ������ˣ�
			return -1;

		// ǰ�ߵ�С�����Ѿ����Ϊ�������������򷵻�-1
		if (cells[x][y + 1][z].getAlreadysearched() != 0) {
			return -1;
		}

		// ��������С������ݶȲ�ֵ
		double xd = 0;// x����ƫ������ֵ��ƽ��

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
	 * ����gapFrontCell 
	 * ���㵱ǰС��������ǰ��С����Ĳ�ֵ 
	 * Parameters: x ��ǰС�����x����(����) 
	 * y ��ǰС�����y����
	 * z ��ǰС�����z���� 
	 * return ��ǰС��������ǰ��С����Ĳ�ֵ 
	 * author weil
	 */
	public double gapFrontCell(int x, int y, int z) {

		if (z == 0)// ��ǰ�ߵ�����û��ǰ�ߵ������ˣ�
			return -1;

		// ǰ�ߵ�С�����Ѿ����Ϊ�������������򷵻�-1
		if (cells[x][y][z - 1].getAlreadysearched() != 0) {
			return -1;
		}

		// ��������С������ݶȲ�ֵ
		double xd = 0;// x����ƫ������ֵ��ƽ��

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
	 * ����gapBackCell 
	 * ���㵱ǰС����������С����Ĳ�ֵ 
	 * Parameters: x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return ��ǰС����������С����Ĳ�ֵ 
	 * author weil
	 */
	public double gapBackCell(int x, int y, int z) {

		if (z == K - 1)// ���ߵ�����û�к�ߵ������ˣ�
			return -1;

		// ��ߵ�С�����Ѿ����Ϊ�������������򷵻�-1
		if (cells[x][y][z + 1].getAlreadysearched() != 0) {
			return -1;
		}

		// ��������С������ݶȲ�ֵ
		double xd = 0;// x����ƫ������ֵ��ƽ��

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
	 * ����partialDerivative() 
	 * ����ÿ��С�����xƫ����yƫ����zƫ�� 
	 * return NULL
	 * author weil
	 */
	public void partialDerivative() {

		// �ǿ�С��������꼯
		ArrayList<Cells> notEmptyGridCoord = new ArrayList<Cells>();
		for (int y = 0; y < K; y++) {
			for (int z = 0; z < K; z++) {
				for (int x = 0; x < K; x++) {
					if (cells[x][y][z].pNum >= 1)
						notEmptyGridCoord.add(new Cells(x, y, z));
				}
			}
		}

		// �������С��������Ӧ���ܶȺ������ı��ļ�GridWeight.txt
		String s = new String();
		try {
			File fo = new File(
					"D:\\Graduation Project\\src\\0331\\GridWeight.txt");

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

						// С�������������
						double xloc = cells[x][y][z].getXloc();
						double yloc = cells[x][y][z].getYloc();
						double zloc = cells[x][y][z].getZloc();

						// С�����Ӧ�Ķ�άͼ�����������ֵ
						double spacialX = cells[x][y][z].getSpacialX();
						double spacialY = cells[x][y][z].getSpacialY();

						// x��y��z�����ƫ����
						double dpx = 0;
						double dpy = 0;
						double dpz = 0;
						double potential = 0;
						double gridWeight = 0;

						// �����ݳ��ƺ����Ķ��壬�����С�������ֵ��ƫ����
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

									// ��������
									double dist2 = dx * dx + dy * dy + dz * dz;

									// �ռ����
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

						// �趨ƫ����ֵ
						cells[x][y][z].setPotential(potential);
						cells[x][y][z].setXdp(dpx);
						cells[x][y][z].setYdp(dpy);
						cells[x][y][z].setZdp(dpz);
						//System.out.println("x����ƫ��ֵ dpx = "+dpx);
						//System.out.println("y����ƫ��ֵ dpy = "+dpy);
						//System.out.println("z����ƫ��ֵ dpz = "+dpz);
					}
				}
			}

			outputfo.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	//�����Ǽ��㵱ǰС��������Χ����ľ���
	
	/**
	 * ����disLeftCell 
	 * ���㵱ǰС�����������С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftCell(int x, int y, int z) {

		if (x == 0)// ����ߵ�����û����ߵ������ˣ�
			return 0;

		double d = 0;
		double xd = 0;
		double yd = 0;
		double zd = 0;
		double Z = 0;//����ֵ

		xd = StrictMath.pow(
				cells[x][y][z].getXloc() - cells[x - 1][y][z].getXloc(), 2.0);

		yd = StrictMath.pow(
				cells[x][y][z].getYloc() - cells[x - 1][y][z].getYloc(), 2.0);

		zd = StrictMath.pow(
				cells[x][y][z].getZloc() - cells[x - 1][y][z].getZloc(), 2.0);

		d = StrictMath.pow((xd + yd + zd), 1.0 / 2);
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y][z] - 0.0) != 0){
			varFactor[x - 1][y][z] = ro	* StrictMath.pow(geoMean / initial[x - 1][y][z], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y][z]);
			
			if (d <= varFactor[x - 1][y][z] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y][z] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y][z],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y][z] = " + varBandwidth[x - 1][y][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y][z].getXloc()) 
										/ varBandwidth[x - 1][y][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y][z].getYloc())
										/ varBandwidth[x - 1][y][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y][z].getZloc())
										/ varBandwidth[x - 1][y][z]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y][z] = StrictMath.pow(varBandwidth[x - 1][y][z], 3.0);
				Z = (initial[x - 1][y][z] * kernel) / varHd[x - 1][y][z];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y][z]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disRightCell 
	 * ���㵱ǰС���������ұ�С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightCell(int x, int y, int z) {

		if (x == K - 1) // ���ұߵ�����û���ұߵ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y][z] - 0.0) != 0){
			varFactor[x + 1][y][z] = ro	* StrictMath.pow(geoMean / initial[x + 1][y][z], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y][z]);
			
			if (d <= varFactor[x + 1][y][z] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y][z] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y][z],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y][z] = " + varBandwidth[x + 1][y][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y][z].getXloc()) 
										/ varBandwidth[x + 1][y][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y][z].getYloc())
										/ varBandwidth[x + 1][y][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y][z].getZloc())
										/ varBandwidth[x + 1][y][z]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y][z] = StrictMath.pow(varBandwidth[x + 1][y][z], 3.0);
				Z = (initial[x + 1][y][z] * kernel) / varHd[x + 1][y][z];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y][z]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disUpCell 
	 * ���㵱ǰС���������ϱ�С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disUpCell(int x, int y, int z) {

		if (y == 0)// ���ϱߵ�����û���ϱߵ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x][y - 1][z] - 0.0) != 0){
			varFactor[x][y - 1][z] = ro	* StrictMath.pow(geoMean / initial[x][y - 1][z], 1.0 / 2);
			System.out.println("��������:" + varFactor[x][y - 1][z]);
			
			if (d <= varFactor[x][y - 1][z] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x][y - 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x][y - 1][z],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x][y - 1][z] = " + varBandwidth[x][y - 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y - 1][z].getXloc()) 
										/ varBandwidth[x][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y - 1][z].getYloc())
										/ varBandwidth[x][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y - 1][z].getZloc())
										/ varBandwidth[x][y - 1][z]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x][y - 1][z] = StrictMath.pow(varBandwidth[x][y - 1][z], 3.0);
				Z = (initial[x][y - 1][z] * kernel) / varHd[x][y - 1][z];
				
				System.out.println(" h(x)��d�η� : " + varHd[x][y - 1][z]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disDownCell 
	 * ���㵱ǰС���������±�С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disDownCell(int x, int y, int z) {

		if (y == K - 1)// ���±ߵ�����û���±ߵ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x][y + 1][z] - 0.0) != 0){
			varFactor[x][y + 1][z] = ro	* StrictMath.pow(geoMean / initial[x][y + 1][z], 1.0 / 2);
			System.out.println("��������:" + varFactor[x][y + 1][z]);
			
			if (d <= varFactor[x][y + 1][z] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x][y + 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x][y + 1][z],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x][y + 1][z] = " + varBandwidth[x][y + 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y + 1][z].getXloc()) 
										/ varBandwidth[x][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y + 1][z].getYloc())
										/ varBandwidth[x][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y + 1][z].getZloc())
										/ varBandwidth[x][y + 1][z]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x][y + 1][z] = StrictMath.pow(varBandwidth[x][y + 1][z], 3.0);
				Z = (initial[x][y + 1][z] * kernel) / varHd[x][y + 1][z];
				
				System.out.println(" h(x)��d�η� : " + varHd[x][y + 1][z]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disFrontCell 
	 * ���㵱ǰС��������ǰ��С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disFrontCell(int x, int y, int z) {

		if (z == 0)// ��ǰ�ߵ�����û��ǰ�ߵ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x][y][z - 1] - 0.0) != 0){
			varFactor[x][y][z - 1] = ro	* StrictMath.pow(geoMean / initial[x][y][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x][y][z - 1]);
			
			if (d <= varFactor[x][y][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x][y][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x][y][z - 1] = " + varBandwidth[x][y][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y][z - 1].getXloc()) 
										/ varBandwidth[x][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y][z - 1].getYloc())
										/ varBandwidth[x][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y][z - 1].getZloc())
										/ varBandwidth[x][y][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x][y][z - 1] = StrictMath.pow(varBandwidth[x][y][z - 1], 3.0);
				Z = (initial[x][y][z - 1] * kernel) / varHd[x][y][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x][y][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disBackCell 
	 * ���㵱ǰС����������С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disBackCell(int x, int y, int z) {

		if (z == K - 1)// ��ǰ�ߵ�����û��ǰ�ߵ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x][y][z + 1] - 0.0) != 0){
			varFactor[x][y][z + 1] = ro	* StrictMath.pow(geoMean / initial[x][y][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x][y][z + 1]);
			
			if (d <= varFactor[x][y][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x][y][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x][y][z + 1] = " + varBandwidth[x][y][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y][z + 1].getXloc()) 
										/ varBandwidth[x][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y][z + 1].getYloc())
										/ varBandwidth[x][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y][z + 1].getZloc())
										/ varBandwidth[x][y][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x][y][z + 1] = StrictMath.pow(varBandwidth[x][y][z + 1], 3.0);
				Z = (initial[x][y][z + 1] * kernel) / varHd[x][y][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x][y][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disLeftUpCell 
	 * ���㵱ǰС�����������ϱ�С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftUpCell(int x, int y, int z) {

		if (x == 0 || y == 0)// �����ϵ�����û�����ϵ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y - 1][z] - 0.0) != 0){
			varFactor[x - 1][y - 1][z] = ro	* StrictMath.pow(geoMean / initial[x - 1][y - 1][z], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y - 1][z]);
			
			if (d <= varFactor[x - 1][y - 1][z] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y - 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y - 1][z],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y - 1][z] = " + varBandwidth[x - 1][y - 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y - 1][z].getXloc()) 
										/ varBandwidth[x - 1][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y - 1][z].getYloc())
										/ varBandwidth[x - 1][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y - 1][z].getZloc())
										/ varBandwidth[x - 1][y - 1][z]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y - 1][z] = StrictMath.pow(varBandwidth[x - 1][y - 1][z], 3.0);
				Z = (initial[x - 1][y - 1][z] * kernel) / varHd[x - 1][y - 1][z];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y - 1][z]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disLeftDownCell 
	 * ���㵱ǰС�����������ϱ�С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftDownCell(int x, int y, int z) {

		if (x == 0 || y == K - 1)// �����µ�����û�����µ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y + 1][z] - 0.0) != 0){
			varFactor[x - 1][y + 1][z] = ro	* StrictMath.pow(geoMean / initial[x - 1][y + 1][z], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y + 1][z]);
			
			if (d <= varFactor[x - 1][y + 1][z] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y + 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y + 1][z],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y + 1][z] = " + varBandwidth[x - 1][y + 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y + 1][z].getXloc()) 
										/ varBandwidth[x - 1][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y + 1][z].getYloc())
										/ varBandwidth[x - 1][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y + 1][z].getZloc())
										/ varBandwidth[x - 1][y + 1][z]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y + 1][z] = StrictMath.pow(varBandwidth[x - 1][y + 1][z], 3.0);
				Z = (initial[x - 1][y + 1][z] * kernel) / varHd[x - 1][y + 1][z];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y + 1][z]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disLeftFrontCell 
	 * ���㵱ǰС����������ǰ��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftFrontCell(int x, int y, int z) {

		if (x == 0 || z == 0)// ����ǰ������û����ǰ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y][z - 1] - 0.0) != 0){
			varFactor[x - 1][y][z - 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y][z - 1]);
			
			if (d <= varFactor[x - 1][y][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y][z - 1] = " + varBandwidth[x - 1][y][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y][z - 1].getXloc()) 
										/ varBandwidth[x - 1][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y][z - 1].getYloc())
										/ varBandwidth[x - 1][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y][z - 1].getZloc())
										/ varBandwidth[x - 1][y][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y][z - 1] = StrictMath.pow(varBandwidth[x - 1][y][z - 1], 3.0);
				Z = (initial[x - 1][y][z - 1] * kernel) / varHd[x - 1][y][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disLeftBackCell 
	 * ���㵱ǰС������������С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftBackCell(int x, int y, int z) {

		if (x == 0 || z == K - 1)// ����������û�����������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y][z + 1] - 0.0) != 0){
			varFactor[x - 1][y][z + 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y][z + 1]);
			
			if (d <= varFactor[x - 1][y][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y][z + 1] = " + varBandwidth[x - 1][y][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y][z + 1].getXloc()) 
										/ varBandwidth[x - 1][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y][z + 1].getYloc())
										/ varBandwidth[x - 1][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y][z + 1].getZloc())
										/ varBandwidth[x - 1][y][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y][z + 1] = StrictMath.pow(varBandwidth[x - 1][y][z + 1], 3.0);
				Z = (initial[x - 1][y][z + 1] * kernel) / varHd[x - 1][y][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disRightUpCell 
	 * ���㵱ǰС�����������ϱ�С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightUpCell(int x, int y, int z) {

		if (x == K - 1 || y == 0)// �����ϵ�����û�����ϵ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y - 1][z] - 0.0) != 0){
			varFactor[x + 1][y - 1][z] = ro	* StrictMath.pow(geoMean / initial[x + 1][y - 1][z], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y - 1][z]);
			
			if (d <= varFactor[x + 1][y - 1][z] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y - 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y - 1][z],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y - 1][z] = " + varBandwidth[x + 1][y - 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y - 1][z].getXloc()) 
										/ varBandwidth[x + 1][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y - 1][z].getYloc())
										/ varBandwidth[x + 1][y - 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y - 1][z].getZloc())
										/ varBandwidth[x + 1][y - 1][z]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y - 1][z] = StrictMath.pow(varBandwidth[x + 1][y - 1][z], 3.0);
				Z = (initial[x + 1][y - 1][z] * kernel) / varHd[x + 1][y - 1][z];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y - 1][z]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disRightDownCell 
	 * ���㵱ǰС�����������±�С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightDownCell(int x, int y, int z) {

		if (x == K - 1 || y == K - 1)// �����µ�����û�����µ������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y + 1][z] - 0.0) != 0){
			varFactor[x + 1][y + 1][z] = ro	* StrictMath.pow(geoMean / initial[x + 1][y + 1][z], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y + 1][z]);
			
			if (d <= varFactor[x + 1][y + 1][z] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y + 1][z] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y + 1][z],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y + 1][z] = " + varBandwidth[x + 1][y + 1][z]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y + 1][z].getXloc()) 
										/ varBandwidth[x + 1][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y + 1][z].getYloc())
										/ varBandwidth[x + 1][y + 1][z]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y + 1][z].getZloc())
										/ varBandwidth[x + 1][y + 1][z]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y + 1][z] = StrictMath.pow(varBandwidth[x + 1][y + 1][z], 3.0);
				Z = (initial[x + 1][y + 1][z] * kernel) / varHd[x + 1][y + 1][z];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y + 1][z]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disRightFrontCell 
	 * ���㵱ǰС����������ǰ��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightFrontCell(int x, int y, int z) {

		if (x == K - 1 || z == 0)// ����ǰ������û����ǰ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y][z - 1] - 0.0) != 0){
			varFactor[x + 1][y][z - 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y][z - 1]);
			
			if (d <= varFactor[x + 1][y][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y][z - 1] = " + varBandwidth[x + 1][y][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y][z - 1].getXloc()) 
										/ varBandwidth[x + 1][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y][z - 1].getYloc())
										/ varBandwidth[x + 1][y][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y][z - 1].getZloc())
										/ varBandwidth[x + 1][y][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y][z - 1] = StrictMath.pow(varBandwidth[x + 1][y][z - 1], 3.0);
				Z = (initial[x + 1][y][z - 1] * kernel) / varHd[x + 1][y][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disRightBackCell 
	 * ���㵱ǰС���������Һ��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightBackCell(int x, int y, int z) {

		if (x == K - 1 || z == K- 1)// ���Һ������û���Һ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y][z + 1] - 0.0) != 0){
			varFactor[x + 1][y][z + 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y][z + 1]);
			
			if (d <= varFactor[x + 1][y][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y][z + 1] = " + varBandwidth[x + 1][y][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y][z + 1].getXloc()) 
										/ varBandwidth[x + 1][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y][z + 1].getYloc())
										/ varBandwidth[x + 1][y][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y][z + 1].getZloc())
										/ varBandwidth[x + 1][y][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y][z + 1] = StrictMath.pow(varBandwidth[x + 1][y][z + 1], 3.0);
				Z = (initial[x + 1][y][z + 1] * kernel) / varHd[x + 1][y][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	
	/**
	 * ����disUpFrontCell 
	 * ���㵱ǰС����������ǰ��С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disUpFrontCell(int x, int y, int z) {

		if (y == 0 || z == 0)// ����ǰ������û����ǰ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x][y - 1][z - 1] - 0.0) != 0){
			varFactor[x][y - 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x][y - 1][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x][y - 1][z - 1]);
			
			if (d <= varFactor[x][y - 1][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x][y - 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y - 1][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x][y - 1][z - 1] = " + varBandwidth[x][y - 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y - 1][z - 1].getXloc()) 
										/ varBandwidth[x][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y - 1][z - 1].getYloc())
										/ varBandwidth[x][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y - 1][z - 1].getZloc())
										/ varBandwidth[x][y - 1][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x][y - 1][z - 1] = StrictMath.pow(varBandwidth[x][y - 1][z - 1], 3.0);
				Z = (initial[x][y - 1][z - 1] * kernel) / varHd[x][y - 1][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x][y - 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disUpBackCell 
	 * ���㵱ǰС���������Ϻ��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disUpBackCell(int x, int y, int z) {

		if (y == 0 || z == K - 1)// ���Ϻ������û���Ϻ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x][y - 1][z + 1] - 0.0) != 0){
			varFactor[x][y - 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x][y - 1][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x][y - 1][z + 1]);
			
			if (d <= varFactor[x][y - 1][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x][y - 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y - 1][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x][y - 1][z + 1] = " + varBandwidth[x][y - 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y - 1][z + 1].getXloc()) 
										/ varBandwidth[x][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y - 1][z + 1].getYloc())
										/ varBandwidth[x][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y - 1][z + 1].getZloc())
										/ varBandwidth[x][y - 1][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x][y - 1][z + 1] = StrictMath.pow(varBandwidth[x][y - 1][z + 1], 3.0);
				Z = (initial[x][y - 1][z + 1] * kernel) / varHd[x][y - 1][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x][y - 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disDownFrontCell 
	 * ���㵱ǰС����������ǰ��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disDownFrontCell(int x, int y, int z) {

		if (y == K - 1 || z == 0)// ����ǰ������û����ǰ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x][y + 1][z - 1] - 0.0) != 0){
			varFactor[x][y + 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x][y + 1][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x][y + 1][z - 1]);
			
			if (d <= varFactor[x][y + 1][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x][y + 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y + 1][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x][y + 1][z - 1] = " + varBandwidth[x][y + 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y + 1][z - 1].getXloc()) 
										/ varBandwidth[x][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y + 1][z - 1].getYloc())
										/ varBandwidth[x][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y + 1][z - 1].getZloc())
										/ varBandwidth[x][y + 1][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x][y + 1][z - 1] = StrictMath.pow(varBandwidth[x][y + 1][z - 1], 3.0);
				Z = (initial[x][y + 1][z - 1] * kernel) / varHd[x][y + 1][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x][y + 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disDownBackCell 
	 * ���㵱ǰС���������º��С���������Ӧ�ܶ�  
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disDownBackCell(int x, int y, int z) {

		if (y == K - 1 || z == K - 1)// ���º������û���º�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x][y + 1][z + 1] - 0.0) != 0){
			varFactor[x][y + 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x][y + 1][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x][y + 1][z + 1]);
			
			if (d <= varFactor[x][y + 1][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x][y + 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x][y + 1][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x][y + 1][z + 1] = " + varBandwidth[x][y + 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x][y + 1][z + 1].getXloc()) 
										/ varBandwidth[x][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x][y + 1][z + 1].getYloc())
										/ varBandwidth[x][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x][y + 1][z + 1].getZloc())
										/ varBandwidth[x][y + 1][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x][y + 1][z + 1] = StrictMath.pow(varBandwidth[x][y + 1][z + 1], 3.0);
				Z = (initial[x][y + 1][z + 1] * kernel) / varHd[x][y + 1][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x][y + 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disLeftUpFrontCell 
	 * ���㵱ǰС������������ǰ��С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftUpFrontCell(int x, int y, int z) {

		if (x == 0 || y == 0 || z == 0)// ������ǰ������û������ǰ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y - 1][z - 1] - 0.0) != 0){
			varFactor[x - 1][y - 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y - 1][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y - 1][z - 1]);
			
			if (d <= varFactor[x - 1][y - 1][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y - 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y - 1][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y - 1][z - 1] = " + varBandwidth[x - 1][y - 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y - 1][z - 1].getXloc()) 
										/ varBandwidth[x - 1][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y - 1][z - 1].getYloc())
										/ varBandwidth[x - 1][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y - 1][z - 1].getZloc())
										/ varBandwidth[x - 1][y - 1][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y - 1][z - 1] = StrictMath.pow(varBandwidth[x - 1][y - 1][z - 1], 3.0);
				Z = (initial[x - 1][y - 1][z - 1] * kernel) / varHd[x - 1][y - 1][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y - 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disLeftUpBackCell 
	 * ���㵱ǰС�����������Ϻ��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftUpBackCell(int x, int y, int z) {

		if (x == 0 || y == 0 || z == K - 1)// �����Ϻ������û�����Ϻ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y - 1][z + 1] - 0.0) != 0){
			varFactor[x - 1][y - 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y - 1][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y - 1][z + 1]);
			
			if (d <= varFactor[x - 1][y - 1][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y - 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y - 1][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y - 1][z + 1] = " + varBandwidth[x - 1][y - 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y - 1][z + 1].getXloc()) 
										/ varBandwidth[x - 1][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y - 1][z + 1].getYloc())
										/ varBandwidth[x - 1][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y - 1][z + 1].getZloc())
										/ varBandwidth[x - 1][y - 1][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y - 1][z + 1] = StrictMath.pow(varBandwidth[x - 1][y - 1][z + 1], 3.0);
				Z = (initial[x - 1][y - 1][z + 1] * kernel) / varHd[x - 1][y - 1][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y - 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disLeftDownFrontCell 
	 * ���㵱ǰС������������ǰ��С���������Ӧ�ܶ�  
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftDownFrontCell(int x, int y, int z) {

		if (x == 0 || y == K - 1 || z == 0)// ������ǰ������û������ǰ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y + 1][z - 1] - 0.0) != 0){
			varFactor[x - 1][y + 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y + 1][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y + 1][z - 1]);
			
			if (d <= varFactor[x - 1][y + 1][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y + 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y + 1][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y + 1][z - 1] = " + varBandwidth[x - 1][y + 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y + 1][z - 1].getXloc()) 
										/ varBandwidth[x - 1][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y + 1][z - 1].getYloc())
										/ varBandwidth[x - 1][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y + 1][z - 1].getZloc())
										/ varBandwidth[x - 1][y + 1][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y + 1][z - 1] = StrictMath.pow(varBandwidth[x - 1][y + 1][z - 1], 3.0);
				Z = (initial[x - 1][y + 1][z - 1] * kernel) / varHd[x - 1][y + 1][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y + 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	
	/**
	 * ����disLeftDownBackCell 
	 * ���㵱ǰС�����������º��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disLeftDownBackCell(int x, int y, int z) {

		if (x == 0 || y == K - 1 || z == K - 1)// �����º������û�����º�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x - 1][y + 1][z + 1] - 0.0) != 0){
			varFactor[x - 1][y + 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x - 1][y + 1][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x - 1][y + 1][z + 1]);
			
			if (d <= varFactor[x - 1][y + 1][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x - 1][y + 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x - 1][y + 1][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x - 1][y + 1][z + 1] = " + varBandwidth[x - 1][y + 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x - 1][y + 1][z + 1].getXloc()) 
										/ varBandwidth[x - 1][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x - 1][y + 1][z + 1].getYloc())
										/ varBandwidth[x - 1][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x - 1][y + 1][z + 1].getZloc())
										/ varBandwidth[x - 1][y + 1][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x - 1][y + 1][z + 1] = StrictMath.pow(varBandwidth[x - 1][y + 1][z + 1], 3.0);
				Z = (initial[x - 1][y + 1][z + 1] * kernel) / varHd[x - 1][y + 1][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x - 1][y + 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	
	
	/**
	 * ����disRightUpFrontCell 
	 * ���㵱ǰС������������ǰ��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightUpFrontCell(int x, int y, int z) {

		if (x == K - 1 || y == 0 || z == 0)// ������ǰ������û������ǰ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y - 1][z - 1] - 0.0) != 0){
			varFactor[x + 1][y - 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y - 1][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y - 1][z - 1]);
			
			if (d <= varFactor[x + 1][y - 1][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y - 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y - 1][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y - 1][z - 1] = " + varBandwidth[x + 1][y - 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y - 1][z - 1].getXloc()) 
										/ varBandwidth[x + 1][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y - 1][z - 1].getYloc())
										/ varBandwidth[x + 1][y - 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y - 1][z - 1].getZloc())
										/ varBandwidth[x + 1][y - 1][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y - 1][z - 1] = StrictMath.pow(varBandwidth[x + 1][y - 1][z - 1], 3.0);
				Z = (initial[x + 1][y - 1][z - 1] * kernel) / varHd[x + 1][y - 1][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y - 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disRightUpBackCell 
	 * ���㵱ǰС�����������Ϻ��С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightUpBackCell(int x, int y, int z) {

		if (x == K - 1 || y == 0 || z == K - 1)// �����Ϻ������û�����Ϻ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y - 1][z + 1] - 0.0) != 0){
			varFactor[x + 1][y - 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y - 1][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y - 1][z + 1]);
			
			if (d <= varFactor[x + 1][y - 1][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y - 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y - 1][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y - 1][z + 1] = " + varBandwidth[x + 1][y - 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y - 1][z + 1].getXloc()) 
										/ varBandwidth[x + 1][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y - 1][z + 1].getYloc())
										/ varBandwidth[x + 1][y - 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y - 1][z + 1].getZloc())
										/ varBandwidth[x + 1][y - 1][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y - 1][z + 1] = StrictMath.pow(varBandwidth[x + 1][y - 1][z + 1], 3.0);
				Z = (initial[x + 1][y - 1][z + 1] * kernel) / varHd[x + 1][y - 1][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y - 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	/**
	 * ����disRightDownFrontCell 
	 * ���㵱ǰС������������ǰ��С���������Ӧ�ܶ�
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightDownFrontCell(int x, int y, int z) {

		if (x == K - 1 || y == K - 1 || z == 0)// ������ǰ������û������ǰ�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y + 1][z - 1] - 0.0) != 0){
			varFactor[x + 1][y + 1][z - 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y + 1][z - 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y + 1][z - 1]);
			
			if (d <= varFactor[x + 1][y + 1][z - 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y + 1][z - 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y + 1][z - 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y + 1][z - 1] = " + varBandwidth[x + 1][y + 1][z - 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y + 1][z - 1].getXloc()) 
										/ varBandwidth[x + 1][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y + 1][z - 1].getYloc())
										/ varBandwidth[x + 1][y + 1][z - 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y + 1][z - 1].getZloc())
										/ varBandwidth[x + 1][y + 1][z - 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y + 1][z - 1] = StrictMath.pow(varBandwidth[x + 1][y + 1][z - 1], 3.0);
				Z = (initial[x + 1][y + 1][z - 1] * kernel) / varHd[x + 1][y + 1][z - 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y + 1][z - 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	
	
	
	/**
	 * ����disRightDownBackCell 
	 * ���㵱ǰС�����������º��С���������Ӧ�ܶ� 
	 * Parameters: 
	 * x ��ǰС�����x����(����) 
	 * y ��ǰС�����y���� 
	 * z ��ǰС�����z���� 
	 * return Z:�ۼ�֮ǰ�Ľ��
	 * author weil
	 */
	public double disRightDownBackCell(int x, int y, int z) {

		if (x == K - 1 || y == K - 1 || z == K - 1)// �����º������û�����º�������ˣ�
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
		
		// ����varFactor����
		//ֻ����ǿ�С����
		if((initial[x + 1][y + 1][z + 1] - 0.0) != 0){
			varFactor[x + 1][y + 1][z + 1] = ro	* StrictMath.pow(geoMean / initial[x + 1][y + 1][z + 1], 1.0 / 2);
			System.out.println("��������:" + varFactor[x + 1][y + 1][z + 1]);
			
			if (d <= varFactor[x + 1][y + 1][z + 1] * width) {

				// ŷʽ����ƽ����
				double dis = 0;
				// �ɱ����
				varBandwidth[x + 1][y + 1][z + 1] = h0 * StrictMath.pow(geoMean	/ initial[x + 1][y + 1][z + 1],	1.0 / 2);
				System.out.println("�ɱ���� varBandwidth[x + 1][y + 1][z + 1] = " + varBandwidth[x + 1][y + 1][z + 1]);
				
				dis = StrictMath.pow(((cells[x][y][z].getXloc() - cells[x + 1][y + 1][z + 1].getXloc()) 
										/ varBandwidth[x + 1][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getYloc() - cells[x + 1][y + 1][z + 1].getYloc())
										/ varBandwidth[x + 1][y + 1][z + 1]), 2.0)
						+ StrictMath.pow(((cells[x][y][z].getZloc() - cells[x + 1][y + 1][z + 1].getZloc())
										/ varBandwidth[x + 1][y + 1][z + 1]), 2.0);
									
				System.out.println("ŷʽ���룺" + dis);
				
				double kernel = 0;

				// ��������||x||<=1
				if (StrictMath.pow(dis, 1.0 / 2) <= 1) {

					kernel = 1 - dis;					

				} else {

					kernel = 0;
				}
				
				System.out.println(" k = " + kernel);
				// h(x)��d�η�
				varHd[x + 1][y + 1][z + 1] = StrictMath.pow(varBandwidth[x + 1][y + 1][z + 1], 3.0);
				Z = (initial[x + 1][y + 1][z + 1] * kernel) / varHd[x + 1][y + 1][z + 1];
				
				System.out.println(" h(x)��d�η� : " + varHd[x + 1][y + 1][z + 1]);
				System.out.println(" Z = " + Z);
				
			}//�ڲ�if����	
		}//���if����
		else {
			Z = 0;
		}
			
		return Z;
	}
	

}
