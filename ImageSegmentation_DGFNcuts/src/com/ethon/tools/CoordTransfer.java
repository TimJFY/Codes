package com.ethon.tools;

import com.ethon.DataBase;
import com.ethon.model.DataPoint;
import com.ethon.ui.DataField;

public class CoordTransfer {
	
	/**
	 * ���ݸ������ݵ㣬��ʼ���߳�ΪlengthOfImagePanel�Ĵ���ͼ������
	 * 
	 * @param points
	 *            ���������ݳ������ݵ�
	 * @param lengthOfImagePanel
	 *            ������ͼ������ı߳�
	 */
	public void initCoord(DataPoint[] points, int lengthOfImagePanel,int heightOfImagePanel ) {
		if (points == null) {
			DataBase.sLen = lengthOfImagePanel;
			DataBase.sHei = heightOfImagePanel;
			DataBase.xMin = 0.0;
			DataBase.yMax = 2.0;
			DataBase.yMin = 0.0;
			DataBase.yMax = 2.0;
//			DataBase.scale = lengthOfImagePanel / 4;
			return;
		}
		double xMin = Float.MAX_VALUE;
		double yMin = Float.MAX_VALUE;
		double xMax = 0;
		double yMax = 0;
		int num = points.length;

		for (int i = 0; i < num; i++) {
			DataPoint p = points[i];
			if (xMin > p.getCoord_X())
				xMin = p.getCoord_X();
			if (yMin > p.getCoord_Y())
				yMin = p.getCoord_Y();

			if (xMax < p.getCoord_X())
				xMax = p.getCoord_X();
			if (yMax < p.getCoord_Y())
				yMax = p.getCoord_Y();
		}

		double side = (xMax - xMin) > (yMax - yMin) ? (xMax - xMin)
				: (yMax - yMin);
/*
		if (side == 0)
			DataBase.scale = 20;
		else
			DataBase.scale = lengthOfImagePanel / (2 * side);
		// else DataBase.scale=lengthOfImagePanel/(side);
*/
		DataBase.xMin = xMin;
		DataBase.xMax = xMax;
		DataBase.yMin = yMin;
		DataBase.yMax = yMax;
		DataBase.sLen = lengthOfImagePanel;
		DataBase.sHei = heightOfImagePanel;
	}
	
	/**
	 * ����ʵ�����꣬������ڴ���ͼ���ϵ�����
	 * 
	 * @param coord_X
	 *            ��ʵ�ʵ�x����
	 * @param coord_Y
	 *            ��ʵ�ʵ�y����
	 * @return ����ͼ���ϵ����� int[0]=x��int[1]=y
	 */
	public int[] getPanelCoord(double coord_X, double coord_Y) {
//		DataBase.scale = 1.0;
		// double x = DataBase.sLen / 2.0 + DataBase.scale
		// * (coord_X - (DataBase.xMin + DataBase.xMax) / 2);
		// double y = DataBase.sLen / 2.0 + DataBase.scale
		// * (coord_Y - (DataBase.yMin + DataBase.yMax) / 2);
		double x = (DataBase.sLen - DataField.line) / 2.0 + coord_X;

//		System.out.println(DataBase.sLen +"   "+ DataField.line);
//		System.out.println(DataBase.sLen +"   "+ DataField.line+"   "+ x + "   "+ coord_X);
		double y = (DataBase.sHei - DataField.row) / 2.0-60 + coord_Y;
//		System.out.println(DataField.line+"  "+"  "+DataField.row+"  ");
//		System.out.println(y);
//		System.out.println(DataBase.sHei +"   "+ DataField.row);

		int[] array = new int[2];
		array[0] = (int) y;
		array[1] = (int) x;

		return array;
	}
	
	public static int getLenOfImagePanel(DataField field) {
		return field.getWidth()/2;
	}
	
	public static int getHeiOfImagePanel(DataField field) {
		return field.getHeight();
	}
}
