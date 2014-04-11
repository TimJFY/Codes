package com.ethon.plugin;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.util.LinkedList;

import com.ethon.DataBase;
import com.ethon.model.DataPoint;
import com.ethon.tools.ImageGenerator;
import com.ethon.ui.DataField;

/**
 * 均值偏移算法
 * @author QiLiang Chen
 * 
 */
public class MeanShift {

	// 记录开始时间
	private long startTime;

	// 记录所有数据的点
	private DataPoint[] points = null;
	private DataPoint[] points_pre = null;
	private DataPoint[] points1 = null;
	DataPoint[][] points_xy;
	DataPoint[][] points_xy_pre;

	// 记录数据点的总数
	// private int allPointsNum;
	int line_num, row_num;
	
	//均值偏移的hs值和hr值
	double hs = 7.0, hr = 8.0;
	
	//记录被平滑的块数的点数的阈值
	int combine_class_point=60;

	// 构造方法,初始化图像，RGB值转换成了LUV值
	public MeanShift(int line_num, int row_num,double hs, double hr, int soomth_value) {

		startTime = System.currentTimeMillis();

		points1 = DataBase.getInstance().getPoints();

		// 记录动态的行与列数
		this.line_num = line_num;
		this.row_num = row_num;
		this.hr = hr;
		this.hs = hs;
		this.combine_class_point = soomth_value;

		System.out.println(points1.length);
		DataBase db = new DataBase();
		points = db.new_public_points(points1);
		points_pre = db.new_public_points(points1);
		// allPointsNum = points.length;

		RgbToLuv(points);
		RgbToLuv(points_pre);
	}

	/**
	 * process方法，经过自己的优化
	 * process方法，mean shift逐渐聚类(1.取出luv3维数据，2.将所有的点存入矩阵 3.初始化j=1;Yi,1 = Xi ;
	 * 3.通过方程式计算与旁边的差值，并转换其luv值4.把luv转换回rgb值)
	 * @author QiLiang Chen
	 * 
	 */
	public void process() {

		startTime = System.currentTimeMillis();
		points_xy = new DataPoint[line_num+1][row_num+1];
		points_xy_pre = new DataPoint[line_num+1][row_num+1];
		

		// 将points所有点放到二维数组points_xy中，以二维坐标轴作为基点

		for (int i = 1; i <= line_num; i++)
			for (int j = 1; j <= row_num; j++)
				if ((i-1) * row_num + j-1 < points.length) {
					points_xy[i][j] = points[(i-1) * row_num + j-1];
					points_xy_pre[i][j] = points[(i-1) * row_num + j-1];
					// System.out.println("正常" + i + "and" + j);
				} else {
					points_xy[i][j] = points[points.length - 1];
					points_xy_pre[i][j] = points[points.length - 1];
				}
		
		smoothing(hs);
		

/*		
		// 调用均值偏移滤波方法进行滤波处理
		double move_range = 15;	
		if(hs<move_range)
		{
			smoothing((int)hs+1);			//调用单个方向的均值偏移，比较范围为hsXhs
			
		}
		else
		{
			smoothing(move_range);			//调用单个方向的均值偏移，比较范围为15X15
		}
*/		
		// 调用classify方法进行分类，然后combineCalsses方法将类进行分块然后去噪（除去点数量较少的块）
		int cluster_num = classify(points_xy, line_num, row_num, hs, hr); // 逐层遍历分类的方法
		int class_num = combineCalsses(points_xy, cluster_num, line_num,
				row_num, combine_class_point);

		System.out.println("6.将points_xy_pre数组中的点存入points");
		// 将points_xy_pre数组中的点存入points
		for (int k = 1; k <= line_num; k++)
			for (int l = 1; l <= row_num; l++)
				if ((k-1) * row_num + l-1 < points.length) {
					points[(k-1) * row_num + l-1] = points_xy[k][l];
					// System.out.println("正常" +k+" and "+l );
				} else {

					points[points.length - 1] = points_xy[k][l];
					// System.out.println("不正常" +k+" and "+l );
				}

		// 将luv装换成xyz然后转换成rgb
		System.out.println("将luv装换成xyz然后转换成rgb");
		LuvToRgb(points);

		// 将数据点初始化输出图像的像素点
		for (DataPoint p : points) {
			Color pointColor = new Color(255, 0, 0);
			if (p.getCoord_X() < 0)
				p.setCoord_X(0.0);
			if (p.getCoord_Y() < 0)
				p.setCoord_Y(0.0);
			if (p.getCoord_Z() < 0)
				p.setCoord_Z(0.0);
			if (p.getCoord_X() > 255)
				p.setCoord_X(255.0);
			if (p.getCoord_Y() > 255)
				p.setCoord_Y(255.0);
			if (p.getCoord_Z() > 255)
				p.setCoord_Z(255.0);
			pointColor = new Color((int) p.getCoord_X(), (int) p.getCoord_Y(),
					(int) p.getCoord_Z());
			p.setColor(pointColor);
		}

		// 激活绘图API, 准备绘图
		System.out.println("激活绘图API, 准备绘图");
		BufferedImage img = ImageGenerator.drawImage(null, DataBase.sLen, DataBase.sHei);
		Graphics g = img.getGraphics();
		long end = System.currentTimeMillis();
		String time = Long.toString(end - startTime) + "ms";
		for (DataPoint p : points) {
			ImageGenerator.show_DataPoint_on_Image(img, p, p.getColor());
		}
		g.setColor(Color.black);
		g.drawString("MS  HS=" + hs + ",HR=" + hr + ",T=" + combine_class_point
				+ " ", 0, DataBase.sHei - 40);
		g.drawString("共" + points.length + "个点，聚成" + class_num + "个块", 0,
				DataBase.sHei - 25);

		g.drawString(time, 0, DataBase.sHei - 10);
		g.dispose();

		g.dispose();
		DataField.updateImagePanel(img);
	}

	// 将rgb装换成xyz然后转换成luv
	public void RgbToLuv(DataPoint[] points) {
		for (DataPoint p : points) {
			// 取出其中的RGB颜色和实际坐标
			double r = p.getCoord_X();
			double g = p.getCoord_Y();
			double b = p.getCoord_Z();
			// System.out.println(r+"  "+g+"   "+b);
			double x, y, z;
			double point_l, point_u, point_v;

			// RGB转换成CIE XYZ再转换成CIE LUV并且存回Datapoints中

			r /= 255.0f;
			g /= 255.0f;
			b /= 255.0f;

			x = 0.412453 * r + 0.357580 * g + 0.180423 * b;
			y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
			z = 0.019334 * r + 0.119193 * g + 0.950227 * b;

			double un = 0, vn = 0, xn = 0.950456, yn = 1, zn = 1.088754;
			un = 4 * xn / (xn + 15 * yn + 3 * zn);
			vn = 9 * yn / (xn + 15 * yn + 3 * zn);

			double u1, v1;
			if ((x + 15 * y + 3 * z) != 0) {
				u1 = 4 * x / (x + 15 * y + 3 * z);
				v1 = 9 * y / (x + 15 * y + 3 * z);
			} else {
				u1 = 0;
				v1 = 0;
			}

			if (y > 0.008856) {
				point_l = 116 * Math.pow(y, 1.0 / 3.0) - 16;
			} else {
				point_l = 116 * (7.787 * y + 16 / 116) - 16;
			}
			point_u = 13 * point_l * (u1 - un);
			point_v = 13 * point_l * (v1 - vn);
			// System.out.println(point_l+"  "+point_u+"   "+point_v);
			// 将luv值覆盖之前的rgb值
			p.setCoord_X(point_l);
			p.setCoord_Y(point_u);
			p.setCoord_Z(point_v);

			// System.out.println(p.getCoord_X()+"   "+p.getCoord_Y()+"   "+
			// p.getCoord_Z());
		}
	}

	// 将luv装换成xyz然后转换成rgb
	public void LuvToRgb(DataPoint[] points) {

		double un = 0, vn = 0, xn = 0.950456, yn = 1, zn = 1.088754;
		un = 4 * xn / (xn + 15 * yn + 3 * zn);
		vn = 9 * yn / (xn + 15 * yn + 3 * zn);

		double l, u, v, x, y, z, r, g, b;
		double u1, v1;

		for (DataPoint p : points) {
			l = p.getCoord_X();
			u = p.getCoord_Y();
			v = p.getCoord_Z();
			// System.out.println(l+"  "+u+"   "+v);
			if (l != 0) {
				u1 = u / (13 * l) + un;
				v1 = v / (13 * l) + vn;
			} else {
				u1 = un;
				v1 = vn;
			}

			if ((l + 16) / 116 > 0.206893) {
				y = Math.pow((l + 16) / 116, 3);
			} else {
				y = ((l + 16) / 116 - 16 / 116) / 7.787;
			}
			if (v1 != 0) {
				x = y * 9 * u1 / (4 * v1);
				z = y * (3 - 0.75 * u1 - 5 * v1) / v1;
			} else {
				x = 0;
				z = 0;
			}

			r = 3.240479 * x - 1.537150 * y - 0.498535 * z;
			g = -0.969256 * x + 1.875992 * y + 0.041556 * z;
			b = 0.055648 * x - 0.204043 * y + 1.057311 * z;

			r = 255 * r;
			g = 255 * g;
			b = 255 * b;
			// System.out.println(r+"  "+g+"   "+b);
			p.setCoord_X(r);
			p.setCoord_Y(g);
			p.setCoord_Z(b);
		}

	}

	// 计算两点XY的距离
	public double compute_distance_xy(DataPoint point1, DataPoint point2) {

		// 计算两点XY的距离
		double d_line = point1.getLine() - point2.getLine();
		double d_row = point1.getRow() - point2.getRow();
		double value = d_line * d_line + d_row * d_row;

		return value;

	}

	// 计算两点LUV的距离
	public double compute_distance_luv(DataPoint point1, DataPoint point2) {

		// 计算两点LUV的距离
		double d_point_l = point1.getCoord_X() - point2.getCoord_X();
		double d_point_u = point1.getCoord_Y() - point2.getCoord_Y();
		double d_point_v = point1.getCoord_Z() - point2.getCoord_Z();
		double value = d_point_l * d_point_l + d_point_u * d_point_u
				+ d_point_v * d_point_v;

		return value;

	}

	// 计算两点LUV的距离
	public double compute_distance(DataPoint point1, DataPoint point2) {

		// 计算两点LUV的距离
		double d_line = point1.getLine() - point2.getLine();
		double d_row = point1.getRow() - point2.getRow();
		double d_point_l = point1.getCoord_X() - point2.getCoord_X();
		double d_point_u = point1.getCoord_Y() - point2.getCoord_Y();
		double d_point_v = point1.getCoord_Z() - point2.getCoord_Z();
		double value = d_line * d_line + d_row * d_row + d_point_l * d_point_l
				+ d_point_u * d_point_u + d_point_v * d_point_v;

		return value;

	}
	
	/**
	 * 滤波的方法
	 * @author QiLiang Chen
	 * 
	 */
	public void smoothing( double hs) 
	{
		// 这里是滤波的原始方法

		int test_count; // 测试用的变量
		int i, j = 1;
		DataPoint point = new DataPoint(0, 0, 0, 0, 0, 1, 1,
				new Color(0, 0, 0), 0, -1, -1, 0, 0, 0); // 记录要比较的点，即Xi
		DataPoint compute_point = new DataPoint(0, 0, 0, 0, 0, 1, 1, new Color(
				0, 0, 0), 0, -1, -1, 0, 0, 0); // 计算g(x)式的中间产值点
		for (i = 1; i <= line_num; i++) {
			System.out.println("现在滤波处理的行数：" + i);
			for (j = 1; j <= row_num; j++) {
				test_count = 0; // 测试用的变量
				double line = 0, row = 0, point_l = 0, point_u = 0, point_v = 0;
				double distance; // 记录两点的距离
				Boolean key = true;
				point = points_xy[i][j];
				int count = 0; // 记录有多少个点是实际有效比较的点
				do {

					int virtual_i = i - 1, virtual_j = j;
					int compare_num_x = 0, compare_num_y = 0;

					while (virtual_i > 1 && compare_num_x < 15) {
						virtual_j = j - 1;
						compare_num_y = 0;
						while (virtual_j >= 1 && compare_num_y < 15) {

							if ((points_xy[virtual_i][virtual_j].getLine() - point
									.getLine()) < hs
									&& (points_xy[virtual_i][virtual_j]
											.getRow() - point.getRow()) < hs) {
								double value_xy = compute_distance_xy(point,
										points_xy[virtual_i][virtual_j]);
								if (value_xy < hs * hs) {
									double value_luv = compute_distance_luv(
											point,
											points_xy[virtual_i][virtual_j]);
									if (value_luv < hr * hr) {
										line += points_xy[virtual_i][virtual_j]
												.getLine();
										row += points_xy[virtual_i][virtual_j]
												.getRow();
										point_l += points_xy[virtual_i][virtual_j]
												.getCoord_X();
										point_u += points_xy[virtual_i][virtual_j]
												.getCoord_Y();
										point_v += points_xy[virtual_i][virtual_j]
												.getCoord_Z();
										count++;
									}
								}
							}
							virtual_j--;
							compare_num_y++;
						}

						// 再次初始化坐标
						virtual_j = j;
						compare_num_y = 0;
						// System.out.println(i +" and " + j);

						while (virtual_j <= row_num && compare_num_y < 16) {
							if ((points_xy[virtual_i][virtual_j].getLine() - point
									.getLine()) < hs
									&& (points_xy[virtual_i][virtual_j]
											.getRow() - point.getRow()) < hs) {
								double value_xy = compute_distance_xy(point,
										points_xy[virtual_i][virtual_j]);
								if (value_xy < hs * hs) {
									double value_luv = compute_distance_luv(
											point,
											points_xy[virtual_i][virtual_j]);
									if (value_luv < hr * hr) {
										line += points_xy[virtual_i][virtual_j]
												.getLine();
										row += points_xy[virtual_i][virtual_j]
												.getRow();
										point_l += points_xy[virtual_i][virtual_j]
												.getCoord_X();
										point_u += points_xy[virtual_i][virtual_j]
												.getCoord_Y();
										point_v += points_xy[virtual_i][virtual_j]
												.getCoord_Z();
										count++;
									}
								}
							}
							virtual_j++;
							compare_num_y++;
						}
						virtual_i--;
						compare_num_x++;
					}

					// 再次初始化虚拟比较的坐标
					virtual_i = i;
					compare_num_x = 0;
					compare_num_y = 0;

					// System.out.println(i +" and " + j +"    "+ virtual_j);
					while (virtual_i <= line_num && compare_num_x < 16) {
						virtual_j = j - 1;
						compare_num_y = 0;
						while (virtual_j > 1 && compare_num_y < 15) {
							// System.out.println(i +" and " + j +"    "+
							// virtual_j);
							if ((points_xy[virtual_i][virtual_j].getLine() - point
									.getLine()) < hs
									&& (points_xy[virtual_i][virtual_j]
											.getRow() - point.getRow()) < hs) {
								double value_xy = compute_distance_xy(point,
										points_xy[virtual_i][virtual_j]);
								if (value_xy < hs * hs) {
									double value_luv = compute_distance_luv(
											point,
											points_xy[virtual_i][virtual_j]);
									if (value_luv < hr * hr) {
										line += points_xy[virtual_i][virtual_j]
												.getLine();
										row += points_xy[virtual_i][virtual_j]
												.getRow();
										point_l += points_xy[virtual_i][virtual_j]
												.getCoord_X();
										point_u += points_xy[virtual_i][virtual_j]
												.getCoord_Y();
										point_v += points_xy[virtual_i][virtual_j]
												.getCoord_Z();
										count++;
									}
								}
							}
							virtual_j--;
							compare_num_y++;
						}

						// 再次初始化坐标
						virtual_j = j;
						compare_num_y = 0;

						while (virtual_j <= row_num && compare_num_y < 16) {
							if ((points_xy[virtual_i][virtual_j].getLine() - point
									.getLine()) < hs
									&& (points_xy[virtual_i][virtual_j]
											.getRow() - point.getRow()) < hs) {
								double value_xy = compute_distance_xy(point,
										points_xy[virtual_i][virtual_j]);
								if (value_xy < hs * hs) {
									double value_luv = compute_distance_luv(
											point,
											points_xy[virtual_i][virtual_j]);
									if (value_luv < hr * hr) {
										line += points_xy[virtual_i][virtual_j]
												.getLine();
										row += points_xy[virtual_i][virtual_j]
												.getRow();
										point_l += points_xy[virtual_i][virtual_j]
												.getCoord_X();
										point_u += points_xy[virtual_i][virtual_j]
												.getCoord_Y();
										point_v += points_xy[virtual_i][virtual_j]
												.getCoord_Z();
										count++;
									}
								}
							}

							virtual_j++;
							compare_num_y++;
						}
						virtual_i++;
						compare_num_x++;
					}

					if (count != 0) {
						compute_point.setLine((line / count)); // 整数与double问题
						compute_point.setRow((row / count));
						compute_point.setCoord_X((point_l / count));
						compute_point.setCoord_Y((point_u / count));
						compute_point.setCoord_Z((point_v / count));

						distance = compute_distance(point, compute_point);

						point.setLine(compute_point.getLine());
						point.setRow(compute_point.getRow());
						point.setCoord_X(compute_point.getCoord_X());
						point.setCoord_Y(compute_point.getCoord_Y());
						point.setCoord_Z(compute_point.getCoord_Z());

						if (distance <= 0.0001) {
							key = false;
						}

					}
					test_count++;
					if (test_count > 500) {
						break;
					}
				} while (key);
			}
		}
	}

	/**
	 * 逐层遍历分类的方法
	 * @author QiLiang Chen
	 * 
	 */
	public int classify(DataPoint[][] points_xy, int line, int row, double hs,
			double hr) {

		DataPoint point = new DataPoint(0, 0, 0, 0, 0, 1, 1,
				new Color(0, 0, 0), 0, -1, -1, 0, 0, 0); // 记录要比较的点，即Xi
		// 计算g(x)式的中间产值点

		// 这里是第一步，初始化第一个点，设置为第一类，并且设置为已搜索
		int cluster_num = 0; // 用于记录类别数
		double distance_xy = 0, distance_luv = 0;

		int range = 15;
		// 用于存放有关系的点的类别的数组
		int[] relation_points = new int[200000];
		for (int i = 0; i < 200000; i++) {
			relation_points[i] = -1;
		}
		int now_num = 0;
		int relation_point_num = 0;

		relation_points[0] = 0;

		for (int i = 1; i <= line; i++) {
			System.out.println("现在是分类处理的行数：" + i);
			for (int j = 1; j <= row; j++) {

				if (points_xy[i][j].getClusterLabel() == -1) {
					points_xy[i][j].setClusterLabel(cluster_num);
					cluster_num++;
				} else {
					continue;
				}

				for (int n = 0; n < 200000; n++) {
					relation_points[n] = -1;
				}
				relation_points[0] = i * row + j;

				// point = points_xy[i][j];
				relation_point_num = 0; // 初始化有关系的点数
				now_num = 0;

				do {
					int num = relation_points[now_num];
					int v_i, v_j;
					if(num % row == 0)
					{
						v_i = num / row - 1;
						v_j = row;
					}
					else
					{
						v_i = num / row;
						v_j = num % row;
					}
					point = points_xy[v_i][v_j];
					now_num++;
					int virtual_i = v_i - 1, virtual_j = v_j;
					int compare_num_x = 0, compare_num_y = 0;

					while (virtual_i > 1 && compare_num_x < range) {
						virtual_j = v_j - 1;
						compare_num_y = 0;
						while (virtual_j >= 1 && compare_num_y < range) {
							if ((points_xy[virtual_i][virtual_j].getLine() - point
									.getLine()) < hs
									&& (points_xy[virtual_i][virtual_j]
											.getRow() - point.getRow()) < hs) {
								distance_xy = compute_distance_xy(point,
										points_xy[virtual_i][virtual_j]);
								if (distance_xy < hs * hs) {
									distance_luv = compute_distance_luv(point,
											points_xy[virtual_i][virtual_j]);
									if (distance_luv < hr * hr) {
										if (points_xy[virtual_i][virtual_j]
												.getClusterLabel() == -1) {
											points_xy[virtual_i][virtual_j]
													.setClusterLabel(point
															.getClusterLabel());
											relation_point_num++;
											relation_points[relation_point_num] = virtual_i
													* row + virtual_j;
										}
									}
								}
							}
							virtual_j--;
							compare_num_y++;
						}

						// 再次初始化坐标
						virtual_j = v_j;
						compare_num_y = 0;
						// System.out.println(i +" and " + j);

						while (virtual_j <= row && compare_num_y < range + 1) {
							if ((points_xy[virtual_i][virtual_j].getLine() - point
									.getLine()) < hs
									&& (points_xy[virtual_i][virtual_j]
											.getRow() - point.getRow()) < hs) {
								distance_xy = compute_distance_xy(point,
										points_xy[virtual_i][virtual_j]);
								if (distance_xy < hs * hs) {
									distance_luv = compute_distance_luv(point,
											points_xy[virtual_i][virtual_j]);
									if (distance_luv < hr * hr) {
										if (points_xy[virtual_i][virtual_j]
												.getClusterLabel() == -1) {
											points_xy[virtual_i][virtual_j]
													.setClusterLabel(point
															.getClusterLabel());
											relation_point_num++;
											relation_points[relation_point_num] = virtual_i
													* row + virtual_j;
										}
									}
								}
							}
							virtual_j++;
							compare_num_y++;
						}
						virtual_i--;
						compare_num_x++;
					}

					// 再次初始化虚拟比较的坐标
					virtual_i = v_i;
					compare_num_x = 0;
					compare_num_y = 0;

					// System.out.println(i +" and " + j +"    "+ virtual_j);
					while (virtual_i <= line && compare_num_x < range + 1) {
						virtual_j = v_j - 1;
						compare_num_y = 0;
						while (virtual_j > 1 && compare_num_y < range) {
							// System.out.println(i +" and " + j +"    "+
							// virtual_j);
							// if(virtual_j == 323) System.out.println(i
							// +" and " + j);
							if ((points_xy[virtual_i][virtual_j].getLine() - point
									.getLine()) < hs
									&& (points_xy[virtual_i][virtual_j]
											.getRow() - point.getRow()) < hs) {
								distance_xy = compute_distance_xy(point,
										points_xy[virtual_i][virtual_j]);
								if (distance_xy < hs * hs) {
									distance_luv = compute_distance_luv(point,
											points_xy[virtual_i][virtual_j]);
									if (distance_luv < hr * hr) {
										if (points_xy[virtual_i][virtual_j]
												.getClusterLabel() == -1) {
											points_xy[virtual_i][virtual_j]
													.setClusterLabel(point
															.getClusterLabel());
											relation_point_num++;
											relation_points[relation_point_num] = virtual_i
													* row + virtual_j;
										}
									}
								}
							}
							virtual_j--;
							compare_num_y++;
						}

						// 再次初始化坐标
						virtual_j = v_j;
						compare_num_y = 0;

						while (virtual_j <= row && compare_num_y < range + 1) {
							if ((points_xy[virtual_i][virtual_j].getLine() - point
									.getLine()) < hs
									&& (points_xy[virtual_i][virtual_j]
											.getRow() - point.getRow()) < hs) {
								distance_xy = compute_distance_xy(point,
										points_xy[virtual_i][virtual_j]);
								if (distance_xy < hs * hs) {
									distance_luv = compute_distance_luv(point,
											points_xy[virtual_i][virtual_j]);
									if (distance_luv < hr * hr) {
										if (points_xy[virtual_i][virtual_j]
												.getClusterLabel() == -1) {
											points_xy[virtual_i][virtual_j]
													.setClusterLabel(point
															.getClusterLabel());
											relation_point_num++;
											relation_points[relation_point_num] = virtual_i
													* row + virtual_j;
										}
									}
								}
							}
							virtual_j++;
							compare_num_y++;
						}
						virtual_i++;
						compare_num_x++;
					}
				} while (now_num <= relation_point_num);
			}
		}

		// 这里是第四步，取出每一类的点，计算平均值,并将平均值存回指定点
		System.out.println("4.将分类好的所有点放回原先的坐标中");
		// 将分类好的所有点放回原先的坐标中
		for (int k = 1; k <= line; k++) {
			for (int l = 1; l <= row; l++) {
				points_xy[k][l].setLine(k);
				points_xy[k][l].setRow(l);
			}
		}
		System.out.println("5.计算每类的平均luv值");
		// 5.计算每类的平均luv值
		double l_last = 0, u_last = 0, v_last = 0;
		int a_clusterLabel_num = 0;
		// System.out.println("总共有 "+ cluster_num +" 类");
		for (int now_count = 0; now_count < cluster_num; now_count++) {
			System.out.println("总共有" + cluster_num + "类，目前计算第" + now_count
					+ "类的平均luv值");
			l_last = 0;
			u_last = 0;
			v_last = 0;
			a_clusterLabel_num = 0;
			for (int k = 1; k <= line; k++) {
				for (int l = 1; l <= row; l++) {
					// 把相同类的luv值加在一起
					if (points_xy[k][l].getClusterLabel() == now_count) {
						l_last += points_xy[k][l].getCoord_X();
						u_last += points_xy[k][l].getCoord_Y();
						v_last += points_xy[k][l].getCoord_Z();
						a_clusterLabel_num++;
					}

				}
			}

			// 计算luv的平均值
			l_last /= a_clusterLabel_num;
			u_last /= a_clusterLabel_num;
			v_last /= a_clusterLabel_num;

			// 把平均的luv值赋值回去
			for (int k = 1; k <= line; k++) {
				for (int l = 1; l <= row; l++) {
					if (points_xy[k][l].getClusterLabel() == now_count) {
						points_xy[k][l].setCoord_X(l_last);
						points_xy[k][l].setCoord_Y(u_last);
						points_xy[k][l].setCoord_Z(v_last);
					}

				}
			}
		}

		return cluster_num;
	}
 
	/**
	 * 先对类进行分块，把相对比较少的块合并到附近的块中(跟那个颜色比较接近就合并到哪里)
	 * @author QiLiang Chen
	 * 
	 */
	public int combineCalsses(DataPoint[][] points_xy, int clusters_num,
			int line, int row, int combine_class_point) {
		DataPoint point = new DataPoint(0, 0, 0, 0, 0, 1, 1,
				new Color(0, 0, 0), 0, -1, -1, 0, 0, 0);
		DataPoint change_point = new DataPoint(0, 0, 0, 0, 0, 1, 1, new Color(
				0, 0, 0), 0, -1, -1, 0, 0, 0);

		int classes_num = 0; // 记录块的编号
		int combine_num = 0; // 记录被合并的块数

		int c_num = 0; // 存储四邻域的类别
		LinkedList<DataPoint> temp;

		int count; // 用于记录一类中点的个数
		double distance = 0.0, min_distance = Double.MAX_VALUE; // 记录两个点的最小距离

		System.out.println("请等待，现在进行分块中");
		// 对类进行分块1
		int cluster = 0;
		LinkedList<DataPoint> class_point;
		for (int i = 1; i <= line; i++)
			for (int j = 1; j <= row; j++) {
				if (points_xy[i][j].getSegmentChecked() == 0) {
					// System.out.println("分块过程中进行到i=" + i + "   和        j="+
					// j);
					cluster = points_xy[i][j].getClusterLabel();
					int k = 0;
					class_point = new LinkedList<DataPoint>();

					points_xy[i][j].setSegmentChecked(1);
					points_xy[i][j].setSegmentLabel(classes_num);
					class_point.add(points_xy[i][j]);

					do {
						temp = allFourNeighborFeilds(class_point.get(k));
						k++;

						for (DataPoint p : temp) {
							if (p.getClusterLabel() == cluster) {
								if (p.getSegmentChecked() == 0) {
									p.setSegmentChecked(1);
									p.setSegmentLabel(classes_num);
									class_point.add(p);
								}
							}
						}
					} while (k < class_point.size());
					classes_num++;
				}

			}


		int[] all_classes = new int[classes_num]; // 用于存储块的类别
		for (int classes = 0; classes < classes_num; classes++) {
			System.out.println("共" + classes_num + " 块，目前是分块处理到第" + classes
					+ " 块");
			distance = 0;
			c_num = 0;
			min_distance = Double.MAX_VALUE;
			// 初始化all_clusters的值，全为0
			for (int i = 0; i < classes_num; i++) {
				all_classes[i] = 0;
			}
			DataPoint[] class_datapoint = new DataPoint[classes_num];
			// 初始化datapoint的值，全为新的点
			for (int i = 0; i < classes_num; i++) {
				class_datapoint[i] = new DataPoint(0, 0, 0, 0, 0, 1, 1,
						new Color(0, 0, 0), 0, -1, -1, 0, 0, 0);
			}

			count = 0; // 计数器清0
			for (int i = 1; i <= line; i++) {
				for (int j = 1; j <= row; j++) {
					if (points_xy[i][j].getSegmentLabel() == classes) {
						count++;
						point = points_xy[i][j];
					}
				}
			}

			// 小于一定点数的块都被合并掉
			if (count <= combine_class_point) {
				combine_num++;
				for (int i = 1; i <= line; i++) {
					for (int j = 1; j <= row; j++) {
						// 取出当前点的上下左右四邻域，然后添加到相应的数组位置中，记录该类相邻的类有哪些
						if (points_xy[i][j].getSegmentLabel() == classes) {
							// 上邻域
							if (i - 1 >= 1) {
								c_num = points_xy[i - 1][j].getSegmentLabel();
								if (classes != c_num) {
									all_classes[c_num] += 1;
									if (class_datapoint[c_num]
											.getSegmentLabel() == -1)
										class_datapoint[c_num] = points_xy[i - 1][j];
								}

							}

							// 下邻域
							if (i + 1 <= line) {
								c_num = points_xy[i + 1][j].getSegmentLabel();
								if (c_num == -1)
									System.out.println(points_xy[i + 1][j]);
								if (classes != c_num) {
									all_classes[c_num] += 1;
									if (class_datapoint[c_num]
											.getSegmentLabel() == -1)
										class_datapoint[c_num] = points_xy[i + 1][j];
								}

							}

							// 左邻域
							if (j - 1 >= 1) {
								c_num = points_xy[i][j - 1].getSegmentLabel();
								if (classes != c_num) {
									all_classes[c_num] += 1;
									if (class_datapoint[c_num]
											.getSegmentLabel() == -1)
										class_datapoint[c_num] = points_xy[i][j - 1];
								}

							}

							// 右邻域
							if (j + 1 <= row) {
								c_num = points_xy[i][j + 1].getSegmentLabel();
								if (classes != c_num) {
									all_classes[c_num] += 1;
									if (class_datapoint[c_num]
											.getSegmentLabel() == -1)
										class_datapoint[c_num] = points_xy[i][j + 1];
								}

							}
						}

					}
				}

				// 找到距离要消除的点数量较少的类的旁边的那一个块的一个点
				for (int i = 0; i < classes_num; i++) {
					if (all_classes[i] != 0) {
						// System.out.println("这里进入了all_clusters[i]");
						distance = compute_distance_luv(point,
								class_datapoint[i]);
						if (min_distance > distance) {
							// System.out.println("这里进入了min_distance > distance");
							min_distance = distance;
							change_point = class_datapoint[i];
						}
					}
				}

				for (int i = 1; i <= line; i++) {
					for (int j = 1; j <= row; j++) {
						if (points_xy[i][j].getSegmentLabel() == classes) {
							points_xy[i][j].setSegmentLabel(change_point
									.getSegmentLabel());
						}
					}
				}

			}
		}

		int now_classes = classes_num - combine_num;
		System.out.println(now_classes);
		int now_cluster_num = 0; // 记录现在要比较的类的类型号码

		// 对新聚类后的点再次重新分好块，类数为now_classes
		for (int i = 0; i < now_classes; i++) {
			here1: for (int k = 1; k <= line; k++) {
				for (int l = 1; l <= row; l++) {
					if (points_xy[k][l].getSegmentLabel() >= i) {
						now_cluster_num = points_xy[k][l].getSegmentLabel();
						break here1;
					}
				}
			}

			for (int k = 1; k <= line; k++)
				for (int l = 1; l <= row; l++) {
					if (points_xy[k][l].getSegmentLabel() == now_cluster_num)
						points_xy[k][l].setSegmentLabel(i);
				}
		}

		// 对分类后的点从新计算luv平均值
		double points_l = 0, points_u = 0, points_v = 0;
		int count_same_classes = 0; // 用于计数
		for (int i = 0; i < now_classes; i++) {
			// 所有的记数都要清0
			points_l = 0;
			points_u = 0;
			points_v = 0;
			count_same_classes = 0;

			for (int k = 1; k <= line; k++) {
				for (int l = 1; l <= row; l++) {
					if (points_xy[k][l].getSegmentLabel() == i) {
						points_l += points_xy[k][l].getCoord_X();
						points_u += points_xy[k][l].getCoord_Y();
						points_v += points_xy[k][l].getCoord_Z();
						count_same_classes++;
					}
				}
			}

			points_l /= count_same_classes;
			points_u /= count_same_classes;
			points_v /= count_same_classes;

			for (int k = 1; k <= line; k++) {
				for (int l = 1; l <= row; l++) {
					if (points_xy[k][l].getSegmentLabel() == i) {
						points_xy[k][l].setCoord_X(points_l);
						points_xy[k][l].setCoord_Y(points_u);
						points_xy[k][l].setCoord_Z(points_v);
					}
				}
			}
		}

		return now_classes;
	}
	

	// 搜索指定数据点在原图像中的四邻域
	public LinkedList<DataPoint> allFourNeighborFeilds(DataPoint point) {
		LinkedList<DataPoint> allFourNeighboor = new LinkedList<DataPoint>();
		// 加入自身
		allFourNeighboor.add(point);
		// 添加上方的点
		if (point.getLine() > 1)
			allFourNeighboor
					.add(points_xy[(int) (point.getLine() - 1)][(int) point
							.getRow()]);
		// 添加左边的点
		if (point.getRow() > 1)
			allFourNeighboor
					.add(points_xy[(int) point.getLine()][(int) (point
							.getRow() - 1)]);
		// 添加右边的点
		if (point.getRow() <= row_num - 1)
			allFourNeighboor
					.add(points_xy[(int) point.getLine()][(int) (point
							.getRow()+1)]);
		// 添加下方的点
		if (point.getLine() <= line_num - 1)
			allFourNeighboor.add(points_xy[(int) (point.getLine()+1)][(int) point
					.getRow()]);

		return allFourNeighboor;

	}
}
