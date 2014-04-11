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
 * ��ֵƫ���㷨
 * @author QiLiang Chen
 * 
 */
public class MeanShift {

	// ��¼��ʼʱ��
	private long startTime;

	// ��¼�������ݵĵ�
	private DataPoint[] points = null;
	private DataPoint[] points_pre = null;
	private DataPoint[] points1 = null;
	DataPoint[][] points_xy;
	DataPoint[][] points_xy_pre;

	// ��¼���ݵ������
	// private int allPointsNum;
	int line_num, row_num;
	
	//��ֵƫ�Ƶ�hsֵ��hrֵ
	double hs = 7.0, hr = 8.0;
	
	//��¼��ƽ���Ŀ����ĵ�������ֵ
	int combine_class_point=60;

	// ���췽��,��ʼ��ͼ��RGBֵת������LUVֵ
	public MeanShift(int line_num, int row_num,double hs, double hr, int soomth_value) {

		startTime = System.currentTimeMillis();

		points1 = DataBase.getInstance().getPoints();

		// ��¼��̬����������
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
	 * process�����������Լ����Ż�
	 * process������mean shift�𽥾���(1.ȡ��luv3ά���ݣ�2.�����еĵ������� 3.��ʼ��j=1;Yi,1 = Xi ;
	 * 3.ͨ������ʽ�������ԱߵĲ�ֵ����ת����luvֵ4.��luvת����rgbֵ)
	 * @author QiLiang Chen
	 * 
	 */
	public void process() {

		startTime = System.currentTimeMillis();
		points_xy = new DataPoint[line_num+1][row_num+1];
		points_xy_pre = new DataPoint[line_num+1][row_num+1];
		

		// ��points���е�ŵ���ά����points_xy�У��Զ�ά��������Ϊ����

		for (int i = 1; i <= line_num; i++)
			for (int j = 1; j <= row_num; j++)
				if ((i-1) * row_num + j-1 < points.length) {
					points_xy[i][j] = points[(i-1) * row_num + j-1];
					points_xy_pre[i][j] = points[(i-1) * row_num + j-1];
					// System.out.println("����" + i + "and" + j);
				} else {
					points_xy[i][j] = points[points.length - 1];
					points_xy_pre[i][j] = points[points.length - 1];
				}
		
		smoothing(hs);
		

/*		
		// ���þ�ֵƫ���˲����������˲�����
		double move_range = 15;	
		if(hs<move_range)
		{
			smoothing((int)hs+1);			//���õ�������ľ�ֵƫ�ƣ��ȽϷ�ΧΪhsXhs
			
		}
		else
		{
			smoothing(move_range);			//���õ�������ľ�ֵƫ�ƣ��ȽϷ�ΧΪ15X15
		}
*/		
		// ����classify�������з��࣬Ȼ��combineCalsses����������зֿ�Ȼ��ȥ�루��ȥ���������ٵĿ飩
		int cluster_num = classify(points_xy, line_num, row_num, hs, hr); // ����������ķ���
		int class_num = combineCalsses(points_xy, cluster_num, line_num,
				row_num, combine_class_point);

		System.out.println("6.��points_xy_pre�����еĵ����points");
		// ��points_xy_pre�����еĵ����points
		for (int k = 1; k <= line_num; k++)
			for (int l = 1; l <= row_num; l++)
				if ((k-1) * row_num + l-1 < points.length) {
					points[(k-1) * row_num + l-1] = points_xy[k][l];
					// System.out.println("����" +k+" and "+l );
				} else {

					points[points.length - 1] = points_xy[k][l];
					// System.out.println("������" +k+" and "+l );
				}

		// ��luvװ����xyzȻ��ת����rgb
		System.out.println("��luvװ����xyzȻ��ת����rgb");
		LuvToRgb(points);

		// �����ݵ��ʼ�����ͼ������ص�
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

		// �����ͼAPI, ׼����ͼ
		System.out.println("�����ͼAPI, ׼����ͼ");
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
		g.drawString("��" + points.length + "���㣬�۳�" + class_num + "����", 0,
				DataBase.sHei - 25);

		g.drawString(time, 0, DataBase.sHei - 10);
		g.dispose();

		g.dispose();
		DataField.updateImagePanel(img);
	}

	// ��rgbװ����xyzȻ��ת����luv
	public void RgbToLuv(DataPoint[] points) {
		for (DataPoint p : points) {
			// ȡ�����е�RGB��ɫ��ʵ������
			double r = p.getCoord_X();
			double g = p.getCoord_Y();
			double b = p.getCoord_Z();
			// System.out.println(r+"  "+g+"   "+b);
			double x, y, z;
			double point_l, point_u, point_v;

			// RGBת����CIE XYZ��ת����CIE LUV���Ҵ��Datapoints��

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
			// ��luvֵ����֮ǰ��rgbֵ
			p.setCoord_X(point_l);
			p.setCoord_Y(point_u);
			p.setCoord_Z(point_v);

			// System.out.println(p.getCoord_X()+"   "+p.getCoord_Y()+"   "+
			// p.getCoord_Z());
		}
	}

	// ��luvװ����xyzȻ��ת����rgb
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

	// ��������XY�ľ���
	public double compute_distance_xy(DataPoint point1, DataPoint point2) {

		// ��������XY�ľ���
		double d_line = point1.getLine() - point2.getLine();
		double d_row = point1.getRow() - point2.getRow();
		double value = d_line * d_line + d_row * d_row;

		return value;

	}

	// ��������LUV�ľ���
	public double compute_distance_luv(DataPoint point1, DataPoint point2) {

		// ��������LUV�ľ���
		double d_point_l = point1.getCoord_X() - point2.getCoord_X();
		double d_point_u = point1.getCoord_Y() - point2.getCoord_Y();
		double d_point_v = point1.getCoord_Z() - point2.getCoord_Z();
		double value = d_point_l * d_point_l + d_point_u * d_point_u
				+ d_point_v * d_point_v;

		return value;

	}

	// ��������LUV�ľ���
	public double compute_distance(DataPoint point1, DataPoint point2) {

		// ��������LUV�ľ���
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
	 * �˲��ķ���
	 * @author QiLiang Chen
	 * 
	 */
	public void smoothing( double hs) 
	{
		// �������˲���ԭʼ����

		int test_count; // �����õı���
		int i, j = 1;
		DataPoint point = new DataPoint(0, 0, 0, 0, 0, 1, 1,
				new Color(0, 0, 0), 0, -1, -1, 0, 0, 0); // ��¼Ҫ�Ƚϵĵ㣬��Xi
		DataPoint compute_point = new DataPoint(0, 0, 0, 0, 0, 1, 1, new Color(
				0, 0, 0), 0, -1, -1, 0, 0, 0); // ����g(x)ʽ���м��ֵ��
		for (i = 1; i <= line_num; i++) {
			System.out.println("�����˲������������" + i);
			for (j = 1; j <= row_num; j++) {
				test_count = 0; // �����õı���
				double line = 0, row = 0, point_l = 0, point_u = 0, point_v = 0;
				double distance; // ��¼����ľ���
				Boolean key = true;
				point = points_xy[i][j];
				int count = 0; // ��¼�ж��ٸ�����ʵ����Ч�Ƚϵĵ�
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

						// �ٴγ�ʼ������
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

					// �ٴγ�ʼ������Ƚϵ�����
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

						// �ٴγ�ʼ������
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
						compute_point.setLine((line / count)); // ������double����
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
	 * ����������ķ���
	 * @author QiLiang Chen
	 * 
	 */
	public int classify(DataPoint[][] points_xy, int line, int row, double hs,
			double hr) {

		DataPoint point = new DataPoint(0, 0, 0, 0, 0, 1, 1,
				new Color(0, 0, 0), 0, -1, -1, 0, 0, 0); // ��¼Ҫ�Ƚϵĵ㣬��Xi
		// ����g(x)ʽ���м��ֵ��

		// �����ǵ�һ������ʼ����һ���㣬����Ϊ��һ�࣬��������Ϊ������
		int cluster_num = 0; // ���ڼ�¼�����
		double distance_xy = 0, distance_luv = 0;

		int range = 15;
		// ���ڴ���й�ϵ�ĵ����������
		int[] relation_points = new int[200000];
		for (int i = 0; i < 200000; i++) {
			relation_points[i] = -1;
		}
		int now_num = 0;
		int relation_point_num = 0;

		relation_points[0] = 0;

		for (int i = 1; i <= line; i++) {
			System.out.println("�����Ƿ��ദ���������" + i);
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
				relation_point_num = 0; // ��ʼ���й�ϵ�ĵ���
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

						// �ٴγ�ʼ������
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

					// �ٴγ�ʼ������Ƚϵ�����
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

						// �ٴγ�ʼ������
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

		// �����ǵ��Ĳ���ȡ��ÿһ��ĵ㣬����ƽ��ֵ,����ƽ��ֵ���ָ����
		System.out.println("4.������õ����е�Ż�ԭ�ȵ�������");
		// ������õ����е�Ż�ԭ�ȵ�������
		for (int k = 1; k <= line; k++) {
			for (int l = 1; l <= row; l++) {
				points_xy[k][l].setLine(k);
				points_xy[k][l].setRow(l);
			}
		}
		System.out.println("5.����ÿ���ƽ��luvֵ");
		// 5.����ÿ���ƽ��luvֵ
		double l_last = 0, u_last = 0, v_last = 0;
		int a_clusterLabel_num = 0;
		// System.out.println("�ܹ��� "+ cluster_num +" ��");
		for (int now_count = 0; now_count < cluster_num; now_count++) {
			System.out.println("�ܹ���" + cluster_num + "�࣬Ŀǰ�����" + now_count
					+ "���ƽ��luvֵ");
			l_last = 0;
			u_last = 0;
			v_last = 0;
			a_clusterLabel_num = 0;
			for (int k = 1; k <= line; k++) {
				for (int l = 1; l <= row; l++) {
					// ����ͬ���luvֵ����һ��
					if (points_xy[k][l].getClusterLabel() == now_count) {
						l_last += points_xy[k][l].getCoord_X();
						u_last += points_xy[k][l].getCoord_Y();
						v_last += points_xy[k][l].getCoord_Z();
						a_clusterLabel_num++;
					}

				}
			}

			// ����luv��ƽ��ֵ
			l_last /= a_clusterLabel_num;
			u_last /= a_clusterLabel_num;
			v_last /= a_clusterLabel_num;

			// ��ƽ����luvֵ��ֵ��ȥ
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
	 * �ȶ�����зֿ飬����ԱȽ��ٵĿ�ϲ��������Ŀ���(���Ǹ���ɫ�ȽϽӽ��ͺϲ�������)
	 * @author QiLiang Chen
	 * 
	 */
	public int combineCalsses(DataPoint[][] points_xy, int clusters_num,
			int line, int row, int combine_class_point) {
		DataPoint point = new DataPoint(0, 0, 0, 0, 0, 1, 1,
				new Color(0, 0, 0), 0, -1, -1, 0, 0, 0);
		DataPoint change_point = new DataPoint(0, 0, 0, 0, 0, 1, 1, new Color(
				0, 0, 0), 0, -1, -1, 0, 0, 0);

		int classes_num = 0; // ��¼��ı��
		int combine_num = 0; // ��¼���ϲ��Ŀ���

		int c_num = 0; // �洢����������
		LinkedList<DataPoint> temp;

		int count; // ���ڼ�¼һ���е�ĸ���
		double distance = 0.0, min_distance = Double.MAX_VALUE; // ��¼���������С����

		System.out.println("��ȴ������ڽ��зֿ���");
		// ������зֿ�1
		int cluster = 0;
		LinkedList<DataPoint> class_point;
		for (int i = 1; i <= line; i++)
			for (int j = 1; j <= row; j++) {
				if (points_xy[i][j].getSegmentChecked() == 0) {
					// System.out.println("�ֿ�����н��е�i=" + i + "   ��        j="+
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


		int[] all_classes = new int[classes_num]; // ���ڴ洢������
		for (int classes = 0; classes < classes_num; classes++) {
			System.out.println("��" + classes_num + " �飬Ŀǰ�Ƿֿ鴦����" + classes
					+ " ��");
			distance = 0;
			c_num = 0;
			min_distance = Double.MAX_VALUE;
			// ��ʼ��all_clusters��ֵ��ȫΪ0
			for (int i = 0; i < classes_num; i++) {
				all_classes[i] = 0;
			}
			DataPoint[] class_datapoint = new DataPoint[classes_num];
			// ��ʼ��datapoint��ֵ��ȫΪ�µĵ�
			for (int i = 0; i < classes_num; i++) {
				class_datapoint[i] = new DataPoint(0, 0, 0, 0, 0, 1, 1,
						new Color(0, 0, 0), 0, -1, -1, 0, 0, 0);
			}

			count = 0; // ��������0
			for (int i = 1; i <= line; i++) {
				for (int j = 1; j <= row; j++) {
					if (points_xy[i][j].getSegmentLabel() == classes) {
						count++;
						point = points_xy[i][j];
					}
				}
			}

			// С��һ�������Ŀ鶼���ϲ���
			if (count <= combine_class_point) {
				combine_num++;
				for (int i = 1; i <= line; i++) {
					for (int j = 1; j <= row; j++) {
						// ȡ����ǰ�����������������Ȼ����ӵ���Ӧ������λ���У���¼�������ڵ�������Щ
						if (points_xy[i][j].getSegmentLabel() == classes) {
							// ������
							if (i - 1 >= 1) {
								c_num = points_xy[i - 1][j].getSegmentLabel();
								if (classes != c_num) {
									all_classes[c_num] += 1;
									if (class_datapoint[c_num]
											.getSegmentLabel() == -1)
										class_datapoint[c_num] = points_xy[i - 1][j];
								}

							}

							// ������
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

							// ������
							if (j - 1 >= 1) {
								c_num = points_xy[i][j - 1].getSegmentLabel();
								if (classes != c_num) {
									all_classes[c_num] += 1;
									if (class_datapoint[c_num]
											.getSegmentLabel() == -1)
										class_datapoint[c_num] = points_xy[i][j - 1];
								}

							}

							// ������
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

				// �ҵ�����Ҫ�����ĵ��������ٵ�����Աߵ���һ�����һ����
				for (int i = 0; i < classes_num; i++) {
					if (all_classes[i] != 0) {
						// System.out.println("���������all_clusters[i]");
						distance = compute_distance_luv(point,
								class_datapoint[i]);
						if (min_distance > distance) {
							// System.out.println("���������min_distance > distance");
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
		int now_cluster_num = 0; // ��¼����Ҫ�Ƚϵ�������ͺ���

		// ���¾����ĵ��ٴ����·ֺÿ飬����Ϊnow_classes
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

		// �Է����ĵ���¼���luvƽ��ֵ
		double points_l = 0, points_u = 0, points_v = 0;
		int count_same_classes = 0; // ���ڼ���
		for (int i = 0; i < now_classes; i++) {
			// ���еļ�����Ҫ��0
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
	

	// ����ָ�����ݵ���ԭͼ���е�������
	public LinkedList<DataPoint> allFourNeighborFeilds(DataPoint point) {
		LinkedList<DataPoint> allFourNeighboor = new LinkedList<DataPoint>();
		// ��������
		allFourNeighboor.add(point);
		// ����Ϸ��ĵ�
		if (point.getLine() > 1)
			allFourNeighboor
					.add(points_xy[(int) (point.getLine() - 1)][(int) point
							.getRow()]);
		// �����ߵĵ�
		if (point.getRow() > 1)
			allFourNeighboor
					.add(points_xy[(int) point.getLine()][(int) (point
							.getRow() - 1)]);
		// ����ұߵĵ�
		if (point.getRow() <= row_num - 1)
			allFourNeighboor
					.add(points_xy[(int) point.getLine()][(int) (point
							.getRow()+1)]);
		// ����·��ĵ�
		if (point.getLine() <= line_num - 1)
			allFourNeighboor.add(points_xy[(int) (point.getLine()+1)][(int) point
					.getRow()]);

		return allFourNeighboor;

	}
}
