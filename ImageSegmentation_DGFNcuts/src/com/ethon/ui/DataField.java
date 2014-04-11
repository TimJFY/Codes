package com.ethon.ui;

import java.awt.Button;
import java.awt.Color;
import java.awt.Dialog;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.TextArea;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.UnsupportedLookAndFeelException;

import com.ethon.DataBase;
import com.ethon.io.DataReader;
import com.ethon.model.DataPoint;
import com.ethon.plugin.FHS;
import com.ethon.plugin.LUVtest;
import com.ethon.plugin.MeanShift;
import com.ethon.plugin.MeanShiftAndNcuts;

import com.ethon.tools.CoordTransfer;
import com.ethon.tools.ImageGenerator;

public class DataField extends JFrame {
	private static final long serialVersionUID = 1L;

	private static ImagePanel ipanel;
	static BufferedImage tempImg;
	public static int line = 0, row = 0; 			// 记录行数和列数
	File file;
	static Frame f;
	
	public DataField() {
		// 设置UI--START
		String ui = UIManager.getSystemLookAndFeelClassName();
		LookAndFeelInfo[] info = UIManager.getInstalledLookAndFeels();
		for (LookAndFeelInfo lk : info) {
			if (lk.getName().equals("Nimbus")) {
				ui = lk.getClassName();
				break;
			}
		}
		try {
			UIManager.setLookAndFeel(ui);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		} catch (UnsupportedLookAndFeelException e) {
			e.printStackTrace();
		}
		// 设置UI--END

		// 设置菜单栏--START
		MenuBar mb = new MenuBar();
		setJMenuBar(mb);
		// 设置菜单栏--END

		// 设置窗口大小和位置--START
		Dimension dim = getToolkit().getScreenSize();
		int sw = (int) dim.getWidth() / 9;
		int sh = (int) dim.getHeight() / 9;
		setLocation(sw, sh);
		setSize(sw * 8, sh * 8);
		// 设置窗口大小和位置--END

		int len = CoordTransfer.getLenOfImagePanel(this);
		int hei = CoordTransfer.getHeiOfImagePanel(this);
		tempImg = ImageGenerator.drawImage(null, len , hei);

		// 设置显示面板--START
		ipanel = new ImagePanel(DataField.this);
		add(ipanel);
		// 设置显示面板--END

		setTitle("图像分割处理系统");
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setLocationRelativeTo(null);
	}
	
	
	
	/**
	 * 窗口的菜单栏
	 * 
	 */
	class MenuBar extends JMenuBar {
		private static final long serialVersionUID = 1L;

		MenuBar() {
			JMenu menu1 = new JMenu("文件");
			JMenuItem item11 = new JMenuItem("打开");
			item11.addActionListener(new FileOpenListener(""));
//			JMenuItem item12 = new JMenuItem("保存");
//			item12.addActionListener(new FileSaveListener());
//			JMenuItem item13 = new JMenuItem("另存为");
//			item13.addActionListener(new FileSaveToOtherListener());
			JMenuItem item14 = new JMenuItem("退出");
			item14.addActionListener(new ExitListener());
//
			JMenu menu2 = new JMenu("GDF");
			JMenuItem item21 = new JMenuItem("GDF分割");
			item21.addActionListener(new LUVtestListener("TEST"));
			
			JMenu menu3 = new JMenu("GDF & Ncut");
			JMenuItem item31 = new JMenuItem("GDF&Ncut分割");
			item31.addActionListener(new LUVtestListener("TEST"));

			JMenu menu4 = new JMenu("MS");
			JMenuItem item41 = new JMenuItem("MS分割");
			item41.addActionListener(new MeanShiftListener());
			
			JMenu menu5 = new JMenu("MS & Ncut");
			JMenuItem item51 = new JMenuItem("MS & Ncut分割");
			item51.addActionListener(new MeanShiftAndNcutsListener());
			
			//FHS算法
			JMenu menu6 = new JMenu("FHS");
			JMenuItem item61 = new JMenuItem("FHS分割");
			item61.addActionListener(new FHSListener());
		
			JMenu menu7 = new JMenu("帮助");
			JMenuItem item71 = new JMenuItem("使用说明");
			item71.addActionListener(new HelpListener());
			JMenuItem item72 = new JMenuItem("版权所有");
			item72.addActionListener(new AuthorListener());
			
			menu1.add(item11);
//			menu1.add(item12);
//			menu1.add(item13);
			menu1.add(item14);

			menu2.add(item21);

			menu3.add(item31);
//			
			menu4.add(item41);
			
			menu5.add(item51);
			
			menu6.add(item61);
			
			menu7.add(item71);
			menu7.add(item72);
//			
			add(menu1);
			add(menu2);
			add(menu3);
			add(menu4);
			add(menu5);
			add(menu6);
			add(menu7);
		}
	}
	
	/**
	 * 重置数据场
	 */
	public void resetDataField() {
		
	}
	
	/**
	 * 刷新窗口图像
	 *
	 */
	public static void updateImagePanel(BufferedImage image) {
		tempImg = image;
		ipanel.update(image);
	}
	
	/**
	 * 退出程序
	 * 
	 * 
	 */
	private class ExitListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			System.exit(0);
		}
	}
	
	
	/**
	 * 读取数据文件
	 * 
	 */
	class FileOpenListener extends Frame implements ActionListener {

		myDialog0 myd0;

		FileOpenListener(String s) {
			super(s);
			setLayout(new GridLayout(1, 1));
			setSize(160, 170);
			setBackground(Color.white);
			setVisible(false);
			myd0 = new myDialog0(this, "图像参数设置", false);
			myd0.setLocation(
					(int) (getToolkit().getScreenSize().getWidth() / 2.5),
					(int) (getToolkit().getScreenSize().getHeight() / 2.5));

		}

		public void actionPerformed(ActionEvent e) {

			JFileChooser jfc = new JFileChooser();
			int value = jfc.showOpenDialog(DataField.this);

			if (value == JFileChooser.APPROVE_OPTION) {
				file = jfc.getSelectedFile();
				DataField.this.setTitle("图像分割处理系统--" + file.getName());
			}
			if(file != null)
			{
				myd0.setVisible(true);
			}
		}
	}
	
	class myDialog0 extends Dialog implements ActionListener {
		Button but, but1;
		TextField text1, text2;
		Label lab1, lab2;
		GridLayout gl;

		myDialog0(Frame f, String s, boolean b) {
			super(f, s, b);
			but = new Button("确定");
			but1 = new Button("重置");
			text1 = new TextField(10);
			text2 = new TextField(10);

			text1.setText("350");
			text2.setText("350");

			lab1 = new Label("请输入图像行宽");
			lab2 = new Label("请输入图像列宽");

			gl = new GridLayout(3, 2);
			setLayout(gl);
			setSize(240, 120);
			setVisible(false);
			setModal(false);

			add(lab1);
			add(text1);
			add(lab2);
			add(text2);

			add(but1);
			add(but);

			but.addActionListener(this);
			but1.addActionListener(this);
			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});
		}

		public void actionPerformed(ActionEvent e) {

			if (e.getSource() == but) {

				line = Integer.valueOf(text1.getText());
				row = Integer.valueOf(text2.getText());
				this.setVisible(false);

				resetDataField();
				DataReader reader = new DataReader();

				DataBase.getInstance().clear();

				DataPoint[] points = reader.getDataPoints(file, 0, line, row);// 读取数据

				CoordTransfer ctf = new CoordTransfer();
				int len = CoordTransfer.getLenOfImagePanel(DataField.this);
				int hei = CoordTransfer.getHeiOfImagePanel(DataField.this);
				ctf.initCoord(points, len, hei);

//				BufferedImage image = ImageGenerator.drawImage(null, len, hei);
//				image = ImageGenerator.show_DataPoints_on_Image(points, image);
//				updateImagePanel(image);
												
				BufferedImage img = ImageGenerator.drawImage(null, len , hei);
				Graphics g = img.getGraphics();
				
				DataBase db = new DataBase();
				db.set_points_color(points);
				
				for (DataPoint p : points) {
					ImageGenerator.show_DataPoint_on_Image(img, p, p.getColor());
				}
				g.setColor(Color.black);				
				g.dispose();
				
				ipanel.original_picture(img);

			} else if (e.getSource() == but1) {
				text1.setText("");
				text2.setText("");
			}
		}
	}
	
	//GDF分割算法
	class LUVtestListener extends Frame implements ActionListener {
		myDialog myd;

		LUVtestListener(String s) {
			super(s);
			setLayout(new GridLayout(1, 1));
			setSize(160, 170);
			setBackground(Color.white);
			setVisible(false);
			myd = new myDialog(this, "网格参数及影响因子设置", false);
			myd.setLocation(
					(int) (getToolkit().getScreenSize().getWidth() / 4.2),
					(int) (getToolkit().getScreenSize().getHeight() / 3.5));
		}

		public void actionPerformed(ActionEvent e) {
			myd.text1.setText(line + "");
			myd.text2.setText(row + "");
			myd.text1.setEditable(false);
			myd.text2.setEditable(false);
//			myd.text12.setEditable(false);

			myd.setVisible(true);

		}
	}
	
	class myDialog extends Dialog implements ActionListener {
		Button but, but1, but2;
		TextField text1, text2, text3, text4, text5, text6, text7, text8,
				text9, text10, text11, text12, text13, text14, text15, text16;
		Label lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8, lab9, lab10,
				lab11, lab12, lab13, lab14, lab15, lab16;
		GridLayout gl, gl2;

		myDialog(Frame f, String s, boolean b) {
			super(f, s, b);
			but = new Button("确定");
			but1 = new Button("重置");
			but2 = new Button("查看推荐影响因子");

			text1 = new TextField(5);
			text2 = new TextField(5);
			text3 = new TextField(5);
			text4 = new TextField(5);
			text5 = new TextField(5);
			text6 = new TextField(5);
			text7 = new TextField(5);
			text8 = new TextField(5);
			text9 = new TextField(5);
			text10 = new TextField(5);
			text11 = new TextField(5);
			text12 = new TextField(5);
			text13 = new TextField(5);
			text14 = new TextField(5);
			text15 = new TextField(5);
			text16 = new TextField(5);

			text3.setText("12");
			text4.setText("2");
			text5.setText("4.8054");
			text6.setText("2.6593");
			text7.setText("4.8318");
			text8.setText("70");
			text9.setText("70");
			text10.setText("40");
			text11.setText("15");
			text12.setText("3");
			text13.setText("0.25");
			text14.setText("2.1");
			text15.setText("2.1");
			text16.setText("2.1");

			lab1 = new Label("图像行宽");
			lab2 = new Label("图像列宽");
			lab3 = new Label("请输入单边小网格数");
			lab4 = new Label("大网格维度数");
			lab5 = new Label("特征空间Simga-L");
			lab6 = new Label("特征空间Simga-U");
			lab7 = new Label("特征空间Simga-V");
			lab8 = new Label("坐标空间Simga-X");
			lab9 = new Label("坐标空间Simga-Y");
			lab10 = new Label("聚类平滑阈值");
			lab11 = new Label("分块Simga-Matrix");
			lab12 = new Label("辅助结点个数");
			lab13 = new Label("Ncut分割阈值");
			lab14 = new Label("  * (times of)");
			lab15 = new Label("  * (times of)");
			lab16 = new Label("  * (times of)");

			gl = new GridLayout(9, 4);

			setLayout(gl);
			setSize(500, 300);
			setVisible(false);
			setModal(false);

			add(lab1);
			add(text1);
			add(lab3);
			add(text3);
			add(lab2);
			add(text2);
			add(lab4);
			add(text4);

			add(lab5);
			add(text14);
			add(lab14);
			add(text5);

			add(lab6);
			add(text15);
			add(lab15);
			add(text6);

			add(lab7);
			add(text16);
			add(lab16);
			add(text7);

			add(lab8);
			add(text8);
			add(lab10);
			add(text10);

			add(lab9);
			add(text9);
			add(lab12);
			add(text12);

			add(lab11);
			add(text11);

			add(lab13);
			add(text13);

			add(but2);
			add(but1);
			add(but);

			but.addActionListener(this);
			but1.addActionListener(this);
			but2.addActionListener(this);

			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});
		}

		myDialog(Frame f, String s, boolean b, int i) {
			super(f, s, b);

			Label[] allLabels = new Label[100];

			allLabels[0] = new Label("文件名");
			allLabels[1] = new Label("特征空间Simga-L");
			allLabels[2] = new Label("特征空间Simga-U");
			allLabels[3] = new Label("特征空间Simga-V");
			allLabels[4] = new Label("Boeing777");
			allLabels[5] = new Label("6.2605");
			allLabels[6] = new Label("2.9413");
			allLabels[7] = new Label("6.0131");
			allLabels[8] = new Label("Bryant");
			allLabels[9] = new Label("5.0156");
			allLabels[10] = new Label("4.5163");
			allLabels[11] = new Label("5.8458");
			allLabels[12] = new Label("Bund");
			allLabels[13] = new Label("3.5453");
			allLabels[14] = new Label("2.6086");
			allLabels[15] = new Label("4.3622");
			allLabels[16] = new Label("dog");
			allLabels[17] = new Label("4.2862");
			allLabels[18] = new Label("1.9572");
			allLabels[19] = new Label("1.8298");
			allLabels[20] = new Label("fruit");
			allLabels[21] = new Label("2.1893");
			allLabels[22] = new Label("4.6951");
			allLabels[23] = new Label("4.3831");
			allLabels[24] = new Label("goose");
			allLabels[25] = new Label("3.3087");
			allLabels[26] = new Label("4.4261");
			allLabels[27] = new Label("3.7204");
			allLabels[28] = new Label("Hawii");
			allLabels[29] = new Label("4.8181");
			allLabels[30] = new Label("4.6572");
			allLabels[31] = new Label("3.9113");
			allLabels[32] = new Label("hofn");
			allLabels[33] = new Label("3.6616");
			allLabels[34] = new Label("6.7378");
			allLabels[35] = new Label("1.7220");
			allLabels[36] = new Label("HuangShan");
			allLabels[37] = new Label("4.7032");
			allLabels[38] = new Label("1.8004");
			allLabels[39] = new Label("3.8797");
			allLabels[40] = new Label("ImperialPalace");
			allLabels[41] = new Label("5.2110");
			allLabels[42] = new Label("4.6848");
			allLabels[43] = new Label("3.0038");
			allLabels[44] = new Label("ImperialPalace_2");
			allLabels[45] = new Label("4.6108");
			allLabels[46] = new Label("6.4276");
			allLabels[47] = new Label("5.1219");
			allLabels[48] = new Label("kingfisher");
			allLabels[49] = new Label("4.2029");
			allLabels[50] = new Label("7.1579");
			allLabels[51] = new Label("7.9442");
			allLabels[52] = new Label("Lamborghini");
			allLabels[53] = new Label("5.1733");
			allLabels[54] = new Label("2.7564");
			allLabels[55] = new Label("5.1010");
			allLabels[56] = new Label("litchi");
			allLabels[57] = new Label("3.0391");
			allLabels[58] = new Label("6.3032");
			allLabels[59] = new Label("2.9623");
			allLabels[60] = new Label("Messi");
			allLabels[61] = new Label("4.7611");
			allLabels[62] = new Label("4.8905");
			allLabels[63] = new Label("6.1772");
			allLabels[64] = new Label("peony");
			allLabels[65] = new Label("3.5819");
			allLabels[66] = new Label("8.4991");
			allLabels[67] = new Label("6.2137");
			allLabels[68] = new Label("Pete");
			allLabels[69] = new Label("4.8863");
			allLabels[70] = new Label("2.2866");
			allLabels[71] = new Label("2.5271");
			allLabels[72] = new Label("QinghaiLake");
			allLabels[73] = new Label("5.1536");
			allLabels[74] = new Label("2.1954");
			allLabels[75] = new Label("2.3420");
			allLabels[76] = new Label("QinghaiLake2");
			allLabels[77] = new Label("3.9078");
			allLabels[78] = new Label("2.8950");
			allLabels[79] = new Label("3.0903");
			allLabels[80] = new Label("tulip");
			allLabels[81] = new Label("3.3528");
			allLabels[82] = new Label("10.3355");
			allLabels[83] = new Label("3.3740");
			allLabels[84] = new Label("white house");
			allLabels[85] = new Label("4.8054");
			allLabels[86] = new Label("2.6593");
			allLabels[87] = new Label("4.8318");

			gl = new GridLayout(22, 4);

			setLayout(gl);
			setSize(480, 420);
			setVisible(false);
			setModal(false);

			for (int num = 0; num < 88; num++) {
				add(allLabels[num]);
			}
			// add(txtArea);
			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});

		}

		public void actionPerformed(ActionEvent e) {

			if (e.getSource() == but) {
				int K, IFP, auxiliaryNum;
				double sig1, sig2, sig3, sigX, sigY, sigMatrix, smoothValve, nCutValve;

				K = Integer.valueOf(text3.getText());
				IFP = Integer.valueOf(text4.getText());
				auxiliaryNum = Integer.valueOf(text12.getText());

				sig1 = Double.valueOf(text5.getText()).doubleValue()
						* Double.valueOf(text14.getText()).doubleValue();
				sig2 = Double.valueOf(text6.getText()).doubleValue()
						* Double.valueOf(text15.getText()).doubleValue();
				sig3 = Double.valueOf(text7.getText()).doubleValue()
						* Double.valueOf(text16.getText()).doubleValue();

				sigX = Double.valueOf(text8.getText()).doubleValue();
				sigY = Double.valueOf(text9.getText()).doubleValue();
				sigMatrix = Double.valueOf(text11.getText()).doubleValue();

				smoothValve = Double.valueOf(text10.getText()).doubleValue();
				nCutValve = Double.valueOf(text13.getText()).doubleValue();

				this.setVisible(false);
				LUVtest luvTest = new LUVtest(K, IFP, (int) line, (int) row,
						sig1, sig2, sig3, sigX, sigY, sigMatrix, smoothValve,
						nCutValve, auxiliaryNum);
				luvTest.process();

			} else if (e.getSource() == but1) {
				text3.setText("");
				text4.setText("");
				text5.setText("");
				text6.setText("");
				text7.setText("");
				text8.setText("");
				text8.setText("");
				text9.setText("");
				text10.setText("");
				text11.setText("");
				text13.setText("");
				text14.setText("");
				text15.setText("");
				text16.setText("");

			} else if (e.getSource() == but2) {
				myDialog myDialog1 = new myDialog(DataField.f, "网格参数及影响因子设置",
						false, 1);
				myDialog1.setLocation((int) (getToolkit().getScreenSize()
						.getWidth() / 1.6), (int) (getToolkit().getScreenSize()
						.getHeight() / 3.5));
				myDialog1.setVisible(true);
			}
		}
	}
	
	//MeanShift算法
	class MeanShiftListener extends Frame implements ActionListener {
		MeanShiftDialog ms_dialog;
		
		MeanShiftListener() {
			super();
			setLayout(new GridLayout(1, 1));
			setSize(160, 170);
			setBackground(Color.white);
			setVisible(false);
			ms_dialog = new MeanShiftDialog(this, "均值偏移参数设置", false);
			ms_dialog.setLocation(
					(int) (getToolkit().getScreenSize().getWidth() / 2.5),
					(int) (getToolkit().getScreenSize().getHeight() / 3.5));
		}
		
		public void actionPerformed(ActionEvent arg0) {

			ms_dialog.text1.setText(line + "");
			ms_dialog.text2.setText(row + "");
			ms_dialog.text1.setEditable(false);
			ms_dialog.text2.setEditable(false);

			ms_dialog.setVisible(true);

		}
	}
	
	class MeanShiftDialog extends Dialog implements ActionListener {
		Button but, but1, but2;
		TextField text1, text2, text3, text4, text5;
		Label lab1, lab2, lab3, lab4, lab5;
		GridLayout gl, gl2;

		MeanShiftDialog(Frame f, String s, boolean b) {
			super(f, s, b);
			but = new Button("确定");
			but1 = new Button("重置");
			but2 = new Button("参数填写说明");

			text1 = new TextField(5);
			text2 = new TextField(5);
			text3 = new TextField(5);
			text4 = new TextField(5);
			text5 = new TextField(5);

			text3.setText("7.0");
			text4.setText("8.0");
			text5.setText("60");

			lab1 = new Label("图像行宽");
			lab2 = new Label("图像列宽");
			lab3 = new Label("均值偏移hs值");
			lab4 = new Label("均值偏移hr值");
			lab5 = new Label("聚类平滑阈值");

			gl = new GridLayout(7, 2);

			setLayout(gl);
			setSize(300, 300);
			setVisible(false);
			setModal(false);

			add(lab1);
			add(text1);
			
			add(lab2);
			add(text2);
			
			add(lab3);
			add(text3);
			
			add(lab4);
			add(text4);

			add(lab5);
			add(text5);

			
			add(but1);
			add(but);
			add(but2);

			but2.addActionListener(this);
			but.addActionListener(this);
			but1.addActionListener(this);

			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});
		}

		MeanShiftDialog(Frame f, String s, boolean b,int value) {
			super(f, s, b);
			lab1 = new Label("参数填写说明",Label.CENTER);
			lab2 = new Label("参数hs：计算均值偏移时xy坐标的欧式距离比较值",Label.CENTER);
			lab3 = new Label("参数hr：计算均值偏移时luv的欧式距离比较值",Label.CENTER);
			lab4 = new Label("hs和hr值越大，分割块数越少，图像也就越模糊",Label.CENTER);
			lab5 = new Label("聚类平滑阈值：被平滑的点数阈值，某块点数小于该值则被临块合并",Label.CENTER);
			gl = new GridLayout(5,1);

			setLayout(gl);
			setSize(400, 300);
			setVisible(false);
			setModal(false);

			add(lab1);
			add(lab2);
			add(lab3);
			add(lab4);
			add(lab5);

			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});
		}
		
		public void actionPerformed(ActionEvent e) {

			if (e.getSource() == but) {
				double hs,hr;
				int smoothValve;

				hs = Double.valueOf(text3.getText()).doubleValue();
				hr = Double.valueOf(text4.getText()).doubleValue();
				smoothValve = Integer.valueOf(text5.getText());

				this.setVisible(false);
				MeanShift meanshift = new MeanShift(line, row,hs,hr,smoothValve);
				meanshift.process();

			} else if (e.getSource() == but1) {
				text3.setText("");
				text4.setText("");
				text5.setText("");
			} else if (e.getSource() == but2) {
				MeanShiftDialog ms_p = new MeanShiftDialog(DataField.f, "参数填写说明",
						false, 1);
				ms_p.setLocation((int) (getToolkit().getScreenSize()
						.getWidth() / 2.5), (int) (getToolkit().getScreenSize()
						.getHeight() / 3.5));
				ms_p.setVisible(true);
			}
		}
	}
	
	//MeanShiftAndNcuts算法
	class MeanShiftAndNcutsListener extends Frame implements ActionListener {

			MeanShiftAndNcutsDialog ms_ncut_dialog;
			
			MeanShiftAndNcutsListener() {
				super();
				setLayout(new GridLayout(1, 1));
				setSize(160, 170);
				setBackground(Color.white);
				setVisible(false);
				ms_ncut_dialog = new MeanShiftAndNcutsDialog(this, "均值偏移参数设置", false);
				ms_ncut_dialog.setLocation(
						(int) (getToolkit().getScreenSize().getWidth() / 2.5),
						(int) (getToolkit().getScreenSize().getHeight() / 3.5));
			}
			
			public void actionPerformed(ActionEvent arg0) {
				ms_ncut_dialog.text1.setText(line + "");
				ms_ncut_dialog.text2.setText(row + "");
				ms_ncut_dialog.text1.setEditable(false);
				ms_ncut_dialog.text2.setEditable(false);

				ms_ncut_dialog.setVisible(true);

		}
	}
	
	class MeanShiftAndNcutsDialog extends Dialog implements ActionListener {
		Button but, but1, but2;
		TextField text1, text2, text3, text4, text5, text6;
		Label lab1, lab2, lab3, lab4, lab5, lab6;
		GridLayout gl, gl2;

		MeanShiftAndNcutsDialog(Frame f, String s, boolean b) {
			super(f, s, b);
			but = new Button("确定");
			but1 = new Button("重置");
			but2 = new Button("参数填写说明");

			text1 = new TextField(5);
			text2 = new TextField(5);
			text3 = new TextField(5);
			text4 = new TextField(5);
			text5 = new TextField(5);
			text6 = new TextField(5);

			text3.setText("7.0");
			text4.setText("8.0");
			text5.setText("60");			
			text6.setText("0.25");

			lab1 = new Label("图像行宽");
			lab2 = new Label("图像列宽");
			lab3 = new Label("均值偏移hs值");
			lab4 = new Label("均值偏移hr值");
			lab5 = new Label("聚类平滑阈值");
			lab6 = new Label("Ncut分割阈值");

			gl = new GridLayout(8, 2);

			setLayout(gl);
			setSize(300, 300);
			setVisible(false);
			setModal(false);

			add(lab1);
			add(text1);
			
			add(lab2);
			add(text2);
			
			add(lab3);
			add(text3);
			
			add(lab4);
			add(text4);

			add(lab5);
			add(text5);
			
			add(lab6);
			add(text6);
			
			add(but1);
			add(but);
			add(but2);

			but2.addActionListener(this);
			but.addActionListener(this);
			but1.addActionListener(this);

			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});
		}

		MeanShiftAndNcutsDialog(Frame f, String s, boolean b,int value) {
			super(f, s, b);
			lab1 = new Label("参数填写说明",Label.CENTER);
			lab2 = new Label("参数hs：计算均值偏移时xy坐标的欧式距离比较值",Label.CENTER);
			lab3 = new Label("参数hr：计算均值偏移时luv的欧式距离比较值",Label.CENTER);
			lab4 = new Label("hs和hr值越大，分割块数越少，图像也就越模糊",Label.CENTER);
			lab5 = new Label("聚类平滑阈值：被平滑的点数阈值，某块点数小于该值则被临块合并",Label.CENTER);
			lab6 = new Label("Ncut分割阈值：nuct算法中的参数，该值越大，分割的块数越多",Label.CENTER);
			gl = new GridLayout(5,1);

			setLayout(gl);
			setSize(400, 300);
			setVisible(false);
			setModal(false);

			add(lab1);
			add(lab2);
			add(lab3);
			add(lab4);
			add(lab5);
			add(lab6);

			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});
		}
		
		public void actionPerformed(ActionEvent e) {

			if (e.getSource() == but) {
				double hs,hr;
				int smoothValve;
				double nCutValve;

				hs = Double.valueOf(text3.getText()).doubleValue();
				hr = Double.valueOf(text4.getText()).doubleValue();
				smoothValve = Integer.valueOf(text5.getText());
				nCutValve = Double.valueOf(text6.getText()).doubleValue();

				this.setVisible(false);
				MeanShiftAndNcuts meanShiftAndNcuts = new MeanShiftAndNcuts(line, row,hs,hr,smoothValve, nCutValve);
				meanShiftAndNcuts.process();

			} else if (e.getSource() == but1) {
				text3.setText("");
				text4.setText("");
				text5.setText("");
			} else if (e.getSource() == but2) {
				MeanShiftAndNcutsDialog ms_ncut_p = new MeanShiftAndNcutsDialog(DataField.f, "参数填写说明",
						false, 1);
				ms_ncut_p.setLocation((int) (getToolkit().getScreenSize()
						.getWidth() / 2.5), (int) (getToolkit().getScreenSize()
						.getHeight() / 3.5));
				ms_ncut_p.setVisible(true);
			}
		}
	}
	
	
	
	//使用说明菜单
	class HelpListener extends Frame implements ActionListener{
		
		HelpDialog help;

		HelpListener() {
			super();
			setLayout(new GridLayout(1, 1));
			setSize(160, 170);
			setBackground(Color.white);
			setVisible(false);
			help = new HelpDialog(this, "使用说明", false);
			help.setLocation(
					(int) (getToolkit().getScreenSize().getWidth() / 2.5),
					(int) (getToolkit().getScreenSize().getHeight() / 2.5));

		}
		
		public void actionPerformed(ActionEvent e) {
			help.setVisible(true);
		}
	}
	
	class HelpDialog extends Dialog {
		Label lab1,lab2,lab3,lab4,lab5;
		GridLayout gl;

		HelpDialog(Frame f, String s, boolean b) {
			super(f, s, b);
			lab1 = new Label("使用帮助",Label.CENTER);
			lab2 = new Label("点击文件->打开，然后选择要读入的文件",Label.CENTER);
			lab3 = new Label("必须是txt格式，每行为一个像素点的RGB值",Label.CENTER);
			lab4 = new Label("每个RGB值的中间应该以空格符号隔开",Label.CENTER);
			lab5 = new Label("然后选择你所需要的进行分割处理的算法",Label.CENTER);
			gl = new GridLayout(6,1);

			setLayout(gl);

			setSize(400, 300);
			setVisible(false);
			setModal(false);

			add(lab1);
			add(lab2);
			add(lab3);
			add(lab4);
			add(lab5);

			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});
		}

	}
	
	//版权所有菜单
	class AuthorListener extends Frame implements ActionListener{
		
		AuthorDialog help;

		AuthorListener() {
			super();
			setLayout(new GridLayout(1, 1));
			setSize(160, 170);
			setBackground(Color.white);
			setVisible(false);
			help = new AuthorDialog(this, "版权所有", false);
			help.setLocation(
					(int) (getToolkit().getScreenSize().getWidth() / 2.5),
					(int) (getToolkit().getScreenSize().getHeight() / 2.5));
		}
		
		public void actionPerformed(ActionEvent e) {
			help.setVisible(true);
		}
	}
	
	class AuthorDialog extends Dialog {
		Label lab1,lab2,lab3,lab4;
		GridLayout gl;

		AuthorDialog(Frame f, String s, boolean b) {
			super(f, s, b);
			lab1 = new Label("版权所有",Label.CENTER);
			lab2 = new Label("作者：王树良，李英，尹劲飞，陈其良，李伟",Label.CENTER);
			lab3 = new Label("最终解释权归作者所有",Label.CENTER);
			lab4 = new Label("email：slwang2005@whu.edu.cn",Label.CENTER);
			gl = new GridLayout(5,1);

			setLayout(gl);
			setSize(500, 300);
			setVisible(false);
			setModal(false);

			add(lab1);
			add(lab2);
			add(lab3);
			add(lab4);

			this.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					setVisible(false);
				}
			});
		}

	}
	
	//FHS算法
	class FHSListener implements ActionListener {
		public void actionPerformed(ActionEvent arg0) {

			int h0 = 16;
			String h0_str = JOptionPane.showInputDialog("输入参数h0值", "16");
			if (h0_str != null)
				h0 = Integer.parseInt(h0_str);
			//先调用构造函数
			FHS fhs = new FHS(h0, line, row);
			fhs.process();

		}
	}

}
