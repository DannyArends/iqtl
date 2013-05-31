/**
 * \file View.java
 * \brief Code file containing the view class
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * Copyright (c) 2010 Danny Arends
 * 
 **/
 
import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;

import javax.swing.Timer;

public class View extends Component implements ActionListener {
  private static final long serialVersionUID = 1L;
  private String plottype;
  public int zoomlevel, trait, marker;
  public int plotby = 0;
  public BufferedImage disp;

  View() {
    zoomlevel = 1;
    trait = 0;
    marker = 0;
    plotby = 0;
    this.plottype = "Overview";
    setBackground(Color.white);
    disp = new BufferedImage(10, 10, BufferedImage.TYPE_INT_RGB);
    Timer t = new Timer(25, this);
    t.start();
  }

  public void paint(Graphics g) {
    Graphics2D g2d = (Graphics2D) g;
    synchronized (disp) {
      g2d.drawImage(disp, 0, 0, this);
    }
  }

  public void update_display(Image nextimage) {
    synchronized (disp) {
      disp = (BufferedImage) nextimage;
      System.out.println("Updated display to nextimage");
    }
  }

  public void setPlottype(String plottype) {
    this.plottype = plottype;
  }

  public String getPlottype() {
    return plottype;
  }

  @Override
  public void actionPerformed(ActionEvent e) {
    this.repaint();
  }
}
