package com.ethon.ui;

import java.awt.Graphics;
import java.awt.image.BufferedImage;

import javax.swing.JPanel;

import com.ethon.tools.CoordTransfer;

public class ImagePanel extends JPanel{
	private static final long serialVersionUID = 1L;
	DataField parent;
	private BufferedImage image;
	private BufferedImage image2;
	
	ImagePanel(DataField parent){
		this.parent=parent;
		image=DataField.tempImg;
		image2=DataField.tempImg;
	}
	
	/**
	 * 重载JPanel的方法，刷新panel的时候自动调用
	 */
	protected void paintComponent(Graphics g){
		super.paintComponent(g);
		int slen = CoordTransfer.getLenOfImagePanel(parent);
//		int slen = parent.getWidth()/2;
		int shei = CoordTransfer.getHeiOfImagePanel(parent);
		
//		int xs=(parent.getWidth()-slen)/2;
		g.drawImage(image2, 0, 0, slen-20, shei-60, this);
		g.drawImage(image, slen+10, 0, slen-20, shei-60, this);
		
		g.dispose();
	}
	
	void update(BufferedImage image){
		this.image=image;		
		repaint();
	}
	void original_picture(BufferedImage image){
		this.image2=image;		
		repaint();
	}
}
