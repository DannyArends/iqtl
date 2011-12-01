/**
 * \file ImageControler.java
 * \brief Code file containing the image controler class
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * Copyright (c) 2010 Danny Arends
 * 
 **/
 
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;

public class ImageControler extends Thread {
  private View myview;
  private Model mydata;
  private PlotFunctions myplots;

  ImageControler(PlotFunctions functions, View view, Model data) {
    myview = view;
    mydata = data;
    myplots = functions;
  }

  public void run() {
    System.out.println("Painting to nextimage");
    Graphics2D g2D = null;
    Image nextimage = null;
    try {
      if (myview.getWidth() > 0 && myview.getHeight() > 0) {
        nextimage = new BufferedImage(myview.getWidth(), myview.getHeight(),
            BufferedImage.TYPE_INT_RGB);
      }
      g2D = (Graphics2D) nextimage.getGraphics();
      g2D.setColor(Color.white);
      g2D.fillRect(0, 0, myview.getWidth(), myview.getHeight());
      myplots.header(g2D, myview.getPlottype());
      if (myview.getPlottype().equals("Overview")) {
        myplots.overview(g2D, mydata, myview);
      }
      if (myview.getPlottype().contains("Genotype")) {
        myplots.genotype(g2D, mydata, myview);
      }
      if (myview.getPlottype().contains("Heatmap")) {
        myplots.heatmap(g2D, mydata, myview);
      }
      if (myview.getPlottype().contains("Cis/Trans")) {
        myplots.cistrans(g2D, mydata, myview);
      }
      if (myview.getPlottype().contains("Circle")) {
        myplots.circle(g2D, mydata, myview);
      }
      if (myview.getPlottype().contains("Profile")) {
        myplots.profile(g2D, mydata, myview);
      }
      System.out.println("Thread done updating nettimage");
    } catch (Exception e) {
      System.out.println("Strange race condition " + e.getStackTrace());
    }
    myview.update_display(nextimage);
    myview.repaint();
  }

}
