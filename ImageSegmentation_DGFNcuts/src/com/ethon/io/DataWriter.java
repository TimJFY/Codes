package com.ethon.io;

import java.awt.Color;
import java.io.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;


import com.ethon.model.DataPoint;
import com.ethon.plugin.dataFeildLUV;
import com.ethon.plugin.LUVtest.Cluster;
import com.ethon.plugin.LUVtest.Coord;
import com.ethon.plugin.LUVtest.Grid;

public class DataWriter {

	public static DataPoint[] transformToImagePoints(DataPoint[] points,
			List<Cluster> clusters, Grid[][][] steps, Coord[][][] allCoords,
			int K) {

		DataPoint[] ImagePoints = new DataPoint[points.length];
		// 350*350
		// 所有边界小网格坐标对象
		LinkedList<Coord> allEdgeCoords = new LinkedList<Coord>();
				int clusterLabel;
		for (int i = 0; i < points.length; i++) {
			if (points[i].getClusterLabel() != -1
			) {
				clusterLabel = points[i].getClusterLabel();
			}
			// 依据欧式距离最近原则找回损失点,找回点需满足距离相近的理解条件
			else {
				// 各聚类的中心(centerX,centerY,centerZ)

				double clusterCenterX = 0;
				double clusterCenterY = 0;
				double clusterCenterZ = 0;
				// 损失点与聚类中心间的最短欧式距离
				double minDistance = 0;
				// 归并类号,及满足损失点与聚类中心间的欧式距离最短
				int minDistanceCluster = 0;

				for (int j = 0; j < clusters.size(); j++) {
					double dis = 0;
					Coord centerCoord = clusters.get(j).getClusterCenter();

					clusterCenterX = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getXloc();
					clusterCenterY = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getYloc();
					clusterCenterZ = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getZloc();
					// 计算距离
					dis = (clusterCenterX - points[i].getCoord_X())
							* (clusterCenterX - points[i].getCoord_X())
							+ (clusterCenterY - points[i].getCoord_Y())
							* (clusterCenterY - points[i].getCoord_Y())
							+ (clusterCenterZ - points[i].getCoord_Z())
							* (clusterCenterZ - points[i].getCoord_Z());

					if (j == 0 || dis < minDistance) {
						minDistance = dis;
						minDistanceCluster = j;
					}

				}

				clusterLabel = minDistanceCluster;
			}
			ImagePoints[i] = new DataPoint(points[i].getLine(),
					points[i].getRow(), points[i].getCoord_X(),
					points[i].getCoord_Y(), points[i].getCoord_Z(), 1, 1,
					points[i].getColor(), points[i].getLineNum(), clusterLabel,
					points[i].getSegmentLabel(), points[i].getSmoothed(),
					points[i].getCombined(), points[i].getSegmentChecked());
		}

		// 寻找所有边界网格,其中的点按最近距离重新分类
		for (int x = 0; x < K; x++) {
			for (int y = 0; y < K; y++) {
				for (int z = 0; z < K; z++) {
					int currentCoordClN = allCoords[x][y][z].getClusterNum();
					boolean notEdgeCoord = true;
					if (x > 0 && notEdgeCoord) {
						if (allCoords[x - 1][y][z].getClusterNum() != currentCoordClN)
							notEdgeCoord = false;
					}
					if (x < K - 1 && notEdgeCoord) {
						if (allCoords[x + 1][y][z].getClusterNum() != currentCoordClN)
							notEdgeCoord = false;
					}
					if (y > 0 && notEdgeCoord) {
						if (allCoords[x][y - 1][z].getClusterNum() != currentCoordClN)
							notEdgeCoord = false;
					}
					if (y < K - 1 && notEdgeCoord) {
						if (allCoords[x][y + 1][z].getClusterNum() != currentCoordClN)
							notEdgeCoord = false;
					}
					if (z > 0 && notEdgeCoord) {
						if (allCoords[x][y][z - 1].getClusterNum() != currentCoordClN)
							notEdgeCoord = false;
					}
					if (z < K - 1 && notEdgeCoord) {
						if (allCoords[x][y][z + 1].getClusterNum() != currentCoordClN)
							notEdgeCoord = false;
					}
					if (!notEdgeCoord)
						allEdgeCoords.add(allCoords[x][y][z]);
				}
			}
		}
		
		for (int m = 0; m < allEdgeCoords.size(); m++) {
			ArrayList<DataPoint> edgeCoordPoints = new ArrayList<DataPoint>();
			edgeCoordPoints = steps[allEdgeCoords.get(m).getX()][allEdgeCoords
					.get(m).getY()][allEdgeCoords.get(m).getZ()].getList();
			for (int n = 0; n < edgeCoordPoints.size(); n++) {
				int lineNum = edgeCoordPoints.get(n).getLineNum();
				lineNum--;
				// 依据欧式距离最近原重新为边界网格点分类
				// 各聚类的中心(centerX,centerY,centerZ)
				double clusterCenterX = 0;
				double clusterCenterY = 0;
				double clusterCenterZ = 0;
				// 损失点与聚类中心间的最短欧式距离
				double minDistance = 0;
				// 归并类号,及满足损失点与聚类中心间的欧式距离最短
				int minDistanceCluster = 0;

				for (int j = 0; j < clusters.size(); j++) {
					double dis = 0;
					Coord centerCoord = clusters.get(j).getClusterCenter();

					clusterCenterX = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getXloc();
					clusterCenterY = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getYloc();
					clusterCenterZ = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getZloc();
					// 计算距离
					dis = (clusterCenterX - points[lineNum].getCoord_X())
							* (clusterCenterX - points[lineNum].getCoord_X())
							+ (clusterCenterY - points[lineNum].getCoord_Y())
							* (clusterCenterY - points[lineNum].getCoord_Y())
							+ (clusterCenterZ - points[lineNum].getCoord_Z())
							* (clusterCenterZ - points[lineNum].getCoord_Z());

					if (j == 0 || dis < minDistance) {
						minDistance = dis;
						minDistanceCluster = j;
					}

				}

				clusterLabel = minDistanceCluster;
				ImagePoints[lineNum] = new DataPoint(points[lineNum].getLine(),
						points[lineNum].getRow(), points[lineNum].getCoord_X(),
						points[lineNum].getCoord_Y(),
						points[lineNum].getCoord_Z(), 1, 1,
						points[lineNum].getColor(),
						points[lineNum].getLineNum(), clusterLabel,
						points[lineNum].getSegmentLabel(),
						points[lineNum].getSmoothed(),
						points[lineNum].getCombined(),
						points[lineNum].getSegmentChecked());
			}
		}
		
		return ImagePoints;
	}
	
	public static DataPoint[] transformToImagePoints2(DataPoint[] points,
			List<dataFeildLUV.Cluster> clusters, dataFeildLUV.Grid[][][] steps,
			dataFeildLUV.Coord[][][] allCoords, int K) {

		DataPoint[] ImagePoints = new DataPoint[points.length];
		// 200*330

		// 所有边界小网格坐标对象
		LinkedList<dataFeildLUV.Coord> allEdgeCoords = new LinkedList<dataFeildLUV.Coord>();
		Color pointColor;
		int clusterLabel;
		// 记录按欧式距离找回的损失点的个数
		int lostRestore = 0;
		// 记录按欧式距离重新分类的边界点的个数
		int edgeRestore = 0;

		for (int i = 0; i < points.length; i++) {
			if (points[i].getClusterLabel() != -1
			// || points[i].getClusterLabel() == -1
			) {
				// pointColor = new Color((int) points[i].getCoord_X(),
				// (int) points[i].getCoord_Y(),
				// (int) points[i].getCoord_Z());
				clusterLabel = points[i].getClusterLabel();
			}
			// 依据欧式距离最近原则找回损失点,找回点需满足距离相近的理解条件
			else {
				// 各聚类的中心(centerX,centerY,centerZ)
				double distanceValve = 25 * 25 * 3;

				double clusterCenterX = 0;
				double clusterCenterY = 0;
				double clusterCenterZ = 0;
				// 损失点与聚类中心间的最短欧式距离
				double minDistance = 0;
				// 归并类号,及满足损失点与聚类中心间的欧式距离最短
				int minDistanceCluster = 0;
				// 损失点颜色值(centerX,centerY,centerZ)
				int centerX;
				int centerY;
				int centerZ;

				for (int j = 0; j < clusters.size(); j++) {
					double dis = 0;
					dataFeildLUV.Coord centerCoord = clusters.get(j)
							.getClusterCenter();

					clusterCenterX = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getXloc();
					clusterCenterY = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getYloc();
					clusterCenterZ = steps[centerCoord.getX()][centerCoord
							.getY()][centerCoord.getZ()].getZloc();
					// 计算距离
					dis = (clusterCenterX - points[i].getCoord_X())
							* (clusterCenterX - points[i].getCoord_X())
							+ (clusterCenterY - points[i].getCoord_Y())
							* (clusterCenterY - points[i].getCoord_Y())
							+ (clusterCenterZ - points[i].getCoord_Z())
							* (clusterCenterZ - points[i].getCoord_Z());

					if (j == 0 || dis < minDistance) {
						minDistance = dis;
						minDistanceCluster = j;
					}

				}

				// if (minDistance <= distanceValve) {
				lostRestore++;
				dataFeildLUV.Cluster restoreToCluster = clusters
						.get(minDistanceCluster);
				dataFeildLUV.Coord centerCoord = restoreToCluster
						.getClusterCenter();
				centerX = (int) steps[centerCoord.getX()][centerCoord.getY()][centerCoord
						.getZ()].getXloc();
				centerY = (int) steps[centerCoord.getX()][centerCoord.getY()][centerCoord
						.getZ()].getYloc();
				centerZ = (int) steps[centerCoord.getX()][centerCoord.getY()][centerCoord
						.getZ()].getZloc();

				//pointColor = new Color(centerX, centerY, centerZ);
				clusterLabel = minDistanceCluster;

				points[i].setCoord_X(centerX);
				points[i].setCoord_Y(centerY);
				points[i].setCoord_Z(centerZ);

				// else {
				// pointColor = new Color((int) points[i].getCoord_X(),
				// (int) points[i].getCoord_Y(),
				// (int) points[i].getCoord_Z());
				// clusterLabel = points[i].getClusterLabel();
				// }

			}
			ImagePoints[i] = new DataPoint(points[i].getLine(), points[i]
					.getRow(), points[i].getCoord_X(), points[i].getCoord_Y(),
					points[i].getCoord_Z(), 1, 1, points[i].getColor(),
					points[i].getLineNum(), clusterLabel, points[i]
							.getSegmentLabel(), points[i].getSmoothed(),
					points[i].getCombined(), points[i].getSegmentChecked());
		}
		return ImagePoints;
	}
	
}
