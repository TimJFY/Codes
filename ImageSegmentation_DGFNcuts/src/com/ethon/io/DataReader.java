package com.ethon.io;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import com.ethon.DataBase;
import com.ethon.model.DataPoint;
import com.ethon.tools.ErrorHandler;

/**
 * 
 * @author ethonCHAN �������ݵĶ�ȡ����
 */
public class DataReader {

	// ���ļ��ж�ȡ��groupNum�����ݣ�groupNum��0��ʼ
	public DataPoint[] getDataPoints(File file, int groupNum, int line, int row) {
		if (groupNum < 0) {
			ErrorHandler.log("groupNum ����Ϊ��");
			return null;
		}
		DataPoint[] points = DataBase.getInstance().getPoints();
		int i = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String str = reader.readLine();
			int points_num = line * row; // ��¼�ܵ���

			points = DataBase.getInstance().newPoints(points_num);

			for (int p = 0; p < points_num; p++) {
				points[p] = new DataPoint(0, 0, 0, 0, 0, 1, 1, new Color(0, 0,
						0), 0, -1, -1, 0, 0, 0);
			}

			while (str != null) {
				if (str.length() != 0 && !str.startsWith("//")) {
					if (str.startsWith(".new"))
						break;
					i++;
					points[i - 1] = DataPoint.getNewPoint(str, i, line, row);
					System.out.println("����  -- " + i + " -- ��");
					// points = db.insertPoint(dp);
				}
				str = reader.readLine();
			}
			int length = 0;
			while (groupNum > -1 && str != null) {
				if (str.length() != 0 && !str.startsWith("//")) {
					if (str.startsWith(".new"))
						groupNum--;
					else {
						i++;
						points[length % points.length]
								.update(str, i, line, row);
						length++;
					}
					str = reader.readLine();
				}
			}
			reader.close();
			return points;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
}
