package com.ethon.tools;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;

import javax.swing.JOptionPane;

import com.ethon.DataBase;
import com.ethon.model.DataPoint;
import com.ethon.tools.CoordTransfer;

public class ImageGenerator {
	public static BufferedImage drawImage(BufferedImage bg, int slen, int shei) {
		BufferedImage img = new BufferedImage(slen, shei,
				BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < slen; i++) {
			for (int j = 0; j < shei; j++) {
				img.setRGB(i, j, new Color(255, 255, 255, 0).getRGB());
			}
		}
		if (bg != null) {
			Graphics g = img.getGraphics();
			g.drawImage(bg, 0, 0, slen-20, shei-60, null);
			g.dispose();
		}
		return img;
	}
	
	public static BufferedImage show_DataPoints_on_Image(DataPoint[] points,
			BufferedImage img) {
		if (points == null || points.length == 0) {
			Graphics g = img.getGraphics();
			g.setColor(new Color(0, 0, 0));
			g.drawString("ÔÝÎÞÊý¾Ý", img.getHeight(null) / 2 - 5,
					img.getHeight(null) / 2);
			g.dispose();
		} else {
			Graphics g = img.getGraphics();
			CoordTransfer tfr = new CoordTransfer();
			for (int i = 0; i < points.length; i++) {
				g.setColor(points[i].getColor());
				int[] coord = tfr.getPanelCoord(points[i].getCoord_X(),
						points[i].getCoord_Y());
				g.fillOval(coord[0] - points[i].getRadium(), coord[1]
						- points[i].getRadium(), points[i].getRadium() * 2,
						points[i].getRadium() * 2);
			}
			g.dispose();
		}
		return img;
	}
	
	public static BufferedImage show_DataPoint_on_Image(BufferedImage image,
			DataPoint point, Color color) {
		Graphics g = image.getGraphics();
		
		int[] coord = new CoordTransfer().getPanelCoord(point.getLine(),
				point.getRow());

//		System.out.println(coord[0]+"  "+coord[1]+"");
		g.setColor(color);
		int radium = point.getRadium();
		g.fillOval(coord[0] - radium, coord[1] - radium, radium*2, radium*2);
		g.dispose();
		return image;
	}
}
