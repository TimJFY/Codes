package com.ethon;

import java.awt.Color;

import com.ethon.model.DataPoint;

/**
 * 公共的数据池
 * 
 * @author Alex
 * 
 */
public class DataBase {
	private static DataBase instance = null;
	private DataPoint[] points_public = null;
	private static DataPoint[] points;

	/*
	 * 数据场的相关常量START
	 */
	// 点的坐标的缩放程度=panel的长度/（2*点坐标的较大的边长）
//	public static double scale = 1.0;

	// 点的实际坐标的范围
	public static double xMin = Double.MAX_VALUE;
	public static double xMax = Double.MIN_VALUE;
	public static double yMin = Double.MAX_VALUE;
	public static double yMax = Double.MIN_VALUE;
	public static double zMin = Double.MAX_VALUE;
	public static double zMax = Double.MIN_VALUE;

	// 初始化坐标轴时坐标轴的长度
	public static int sLen;
	
	// 初始化坐标轴时坐标轴的高度
	public static int sHei;

	/*
	 * 数据场的相关常量END
	 */

	public DataBase() {
	}

	public DataPoint[] newPoints(int num) {
		points = new DataPoint[num];

		return points;
	}

	public DataPoint[] new_public_points(DataPoint[] points1) {
		DataPoint[] public_points = new DataPoint[points1.length];
		for (int points_length = 0; points_length < points1.length; points_length++) {
			public_points[points_length] = new DataPoint(0, 0, 0, 0, 0, 1, 1,
					new Color(0, 0, 0), 0, -1, -1, 0, 0, 0);
		}

		for (int points_length = 0; points_length < points1.length; points_length++) {
			public_points[points_length].setLine(points1[points_length]
					.getLine());
			public_points[points_length]
					.setRow(points1[points_length].getRow());
			public_points[points_length].setCoord_X(points1[points_length]
					.getCoord_X());
			public_points[points_length].setCoord_Y(points1[points_length]
					.getCoord_Y());
			public_points[points_length].setCoord_Z(points1[points_length]
					.getCoord_Z());
			public_points[points_length].setLineNum(points1[points_length]
					.getLineNum());

		}
		return public_points;
	}

	public static DataBase getInstance() {
		if (instance == null) {
			instance = new DataBase();
			points = null;
		}
		return instance;
	}

	public DataPoint[] getPoints() {
		return points;
	}

	public void setPoints(DataPoint[] points) {
		DataBase.points = points;
	}

	public boolean pointsIsNull() {
		if (points == null || points.length == 0)
			return true;
		else
			return false;
	}

	public void setPoint(int index, DataPoint point) {
		points[index] = point;
	}

	public DataPoint getPoint(int index) {
		return points[index];
	}

	public void removePoint(int index) {
		int len = points.length;
		if (index < 0 || index > len - 1)
			return;
		DataPoint[] result = new DataPoint[len - 1];
		for (int i = 0, j = 0; i < len; i++) {
			if (i == index)
				continue;
			result[j] = points[i];
			j++;
		}
		points = result;
	}

	public void clear() {
		points = new DataPoint[0];
	}

	public DataPoint[] insertPoint(DataPoint point) {
		if (points == null) {
			points = new DataPoint[1];
			points[0] = point;
		} else {
			int len = points.length;
			DataPoint[] newone = new DataPoint[len + 1];
			for (int i = 0; i < len; i++) {
				newone[i] = points[i];
			}
			newone[len] = point;
			points = newone;
		}

		return points;
	}

	public DataPoint[] insertPoints(int num) {
		points_public = new DataPoint[num];
		return points_public;
	}

	public void set_points_color(DataPoint[] points)
	{
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
	}
}
