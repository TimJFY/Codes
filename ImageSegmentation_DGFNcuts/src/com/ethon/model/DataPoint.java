package com.ethon.model;

import java.awt.Color;

public class DataPoint {

	// 数据点在图片中的坐标(x,y)=[line,row]
	private double line;
	private double row;

	private int lineNum;
	// 数据点的RGB值，[R,G,B]=[coord_X,coord_Y,coord_Z]
	private double coord_X;
	private double coord_Y;
	private double coord_Z;

	// 数据点的权值
	private double weight = 1;

	// 数据点在panel中的半径
	private int radium = 1;

	// 数据点的颜色
	private Color color = new Color(0, 0, 0);

	// 类别标签
	private int clusterLabel;

	// 块号标签
	private int segmentLabel;

	// 是否参与过平滑
	private int smoothed = 0;

	// 是否参与过块搜索
	private int segmentChecked = 0;

	// 是否参与过融合
	private int combined = 0;

	public double getLine() {
		return line;
	}

	public void setLine(double line) {
		this.line = line;
	}

	public double getRow() {
		return row;
	}

	public void setRow(double row) {
		this.row = row;
	}

	public int getSmoothed() {
		return smoothed;
	}

	public void setSmoothed(int smoothed) {
		this.smoothed = smoothed;
	}

	public int getSegmentChecked() {
		return segmentChecked;
	}

	public void setSegmentChecked(int segmentChecked) {
		this.segmentChecked = segmentChecked;
	}

	public int getCombined() {
		return combined;
	}

	public void setCombined(int combined) {
		this.combined = combined;
	}

	public int getClusterLabel() {
		return clusterLabel;
	}

	public void setClusterLabel(int clusterLabel) {
		this.clusterLabel = clusterLabel;
	}

	public int getSegmentLabel() {
		return segmentLabel;
	}

	public void setSegmentLabel(int segmentLabel) {
		this.segmentLabel = segmentLabel;
	}

	public DataPoint(double line, double row, double coord_X, double coord_Y,
			double coord_Z, double weight, int radium, Color color,
			int lineNum, int clusterLabel, int segmentLabel, int smoothed,
			int combined, int segmentChecked) {
		this.line = line;
		this.row = row;

		this.coord_X = coord_X;
		this.coord_Y = coord_Y;
		this.coord_Z = coord_Z;
		this.weight = weight > 0 ? weight : 100;
		this.radium = radium > 0 ? radium : 1;
		this.color = color == null ? new Color(0, 0, 0) : color;
		this.lineNum = lineNum;
		this.clusterLabel = clusterLabel;
		this.segmentLabel = segmentLabel;
		this.smoothed = smoothed;
		this.combined = combined;
		this.segmentChecked = segmentChecked;
	}

	public double getWeight() {
		return weight;
	}

	public void setWeight(double weight) {
		this.weight = weight;
	}

	public int getRadium() {
		return radium;
	}

	public void setRadium(int radium) {
		this.radium = radium;
	}

	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}

	public void setCoord_X(double coord_X) {
		this.coord_X = coord_X;
	}

	public void setCoord_Y(double coord_Y) {
		this.coord_Y = coord_Y;
	}

	public void setCoord_Z(double coord_Z) {
		this.coord_Z = coord_Z;
	}

	public double getCoord_X() {
		return coord_X;
	}

	public double getCoord_Y() {
		return coord_Y;
	}

	public double getCoord_Z() {
		return coord_Z;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null)
			return false;
		if (this.getClass() != obj.getClass())
			return false;

		DataPoint other = (DataPoint) obj;

		if (coord_X != other.coord_X)
			return false;
		if (coord_Y != other.coord_Y)
			return false;
		if (coord_Z != other.coord_Z)
			return false;

		if (lineNum != other.lineNum)
			return false;
		return true;
	}

	public int getLineNum() {
		return lineNum;
	}

	public void setLineNum(int lineNum) {
		this.lineNum = lineNum;
	}

	@Override
	public String toString() {
		return "line= " + line + "," + "row= " + row + "," + "coord X;Y;Z= "
				+ coord_X + "," + coord_Y + "," + coord_Z + "," + "weight= "
				+ weight + "," + "radium= " + radium + ",RGB[][][]= ["
				+ color.getRed() + "," + color.getGreen() + ","
				+ color.getBlue() + "]" + "," + "lineNum= " + lineNum + ","
				+ "clusterLabel= " + clusterLabel;
	}

	/**
	 * 根据读取的字符串更新相应的DataPoint数据
	 */
	public DataPoint update(String str, int lineNum, int line_num, int row_num) {
		DataPoint dp = DataPoint.getNewPoint(str, lineNum, line_num, row_num);

		double x = dp.getCoord_X();
		double y = dp.getCoord_Y();
		double z = dp.getCoord_Z();
		if (x >= 0)
			coord_X = x;
		if (y >= 0)
			coord_Y = y;
		if (z >= 0)
			coord_Z = z;
		weight = dp.getWeight();
		radium = dp.getRadium();
		color = dp.getColor();
		line = dp.getLine();
		row = dp.getRow();
		lineNum = dp.getLineNum();
		return this;
	}

	public static DataPoint getNewPoint(String str, int i, int line_num,
			int row_num) {
		str = str.trim();
		if (str.length() == 0)
			return null;
		double coord_X = -1;
		double coord_Y = -1;
		double coord_Z = -1;
		Color color = null;
		int line;
		int row;

		if ((int) i % row_num == 0) {
			line = (int) i / row_num;
			row = row_num;
		} else {
			line = (int) i / row_num + 1;
			row = i % row_num;
		}

		str = str.replaceAll("\\s+", " ");
		String[] array = str.split(" ");

		if (array.length == 3) {
			System.out.println("Yes! Reading...");
			if (array[0].trim().length() != 0)
				coord_X = Double.parseDouble(array[0].trim());
			if (array[1].trim().length() != 0)
				coord_Y = Double.parseDouble(array[1].trim());
			if (array[2].trim().length() != 0)
				coord_Z = Double.parseDouble(array[2].trim());
		}
		return new DataPoint(line, row, coord_X, coord_Y, coord_Z, 1, 1, color,
				i, -1, -1, 0, 0, 0);
	}

}
