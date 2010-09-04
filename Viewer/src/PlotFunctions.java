import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;

public class PlotFunctions {

  String myTruncate(String toTruncate) {
    int nchar = toTruncate.length();
    if (nchar > 20) {
      return (toTruncate.substring(0, 9) + ".." + toTruncate.substring(
          nchar - 9, nchar));
    } else {
      return (toTruncate);
    }
  }

  Color colorFromLod(double lod, double maxlod) {
    float cR = (float) (0.5 + (lod / (2 * maxlod)));
    float cG = (float) (1 - (lod / maxlod));
    float cB = (float) ((lod / (2 * maxlod)));
    Color c = new Color(cR, cG, cB);
    return c;
  }

  void header(Graphics2D g2D, String Header) {
    g2D.setColor(Color.black);
    g2D.setFont(new Font("sansserif", Font.BOLD, 12));
    g2D.drawString(Header, 10, 12);
    g2D.setFont(new Font("sansserif", Font.PLAIN, 10));
  }

  void overview(Graphics2D g2D, Model mydata, View myview) {
    g2D.setColor(Color.black);
    g2D.drawString("Traits: " + mydata.ntraits, 10, 36);
    g2D.drawString("Chromosomes: " + mydata.nchromosomes, 10, 48);
    String distances = "Lengths: ";
    String s;
    for (int c = 0; c < mydata.nchromosomes; c++) {
      distances += mydata.chrlengths[c] + " ";
    }
    g2D.drawString(distances, 10, 60);
    g2D.drawString("Markers: " + mydata.nmarkers, 10, 72);
    for (int t = 0; t < mydata.ntraits; t++) {
      s = myTruncate(mydata.qtlmatrix[t].name) + " with significant markers: ";
      for (int m = 0; m < mydata.nmarkers; m++) {
        if (mydata.modelmatrix != null) {
          if (mydata.modelmatrix[t].scores[m] == 1) {
            s = s + ", " + mydata.markers[m].name;
          }
        }
      }
      // g2D.drawString(s, 10, 84+12*t);
    }
  }

  void genotype(Graphics2D g2D, Model mydata, View myview) {
    g2D.setColor(Color.black);
    g2D.drawString("Traits: " + mydata.ntraits, 10, 36);
    g2D.drawString("Chromosomes: " + mydata.nchromosomes, 10, 48);
    String distances = "Lengths: ";
    for (int c = 0; c < mydata.nchromosomes; c++) {
      distances += mydata.chrlengths[c] + " ";
    }
    g2D.drawString(distances, 10, 60);
    g2D.drawString("Markers: " + mydata.nmarkers, 10, 72);
  }

  void heatmap(Graphics2D g2D, Model mydata, View myview) {
    double lod;
    String name;
    for (int t = myview.trait; t < mydata.ntraits; t++) {
      name = mydata.qtlmatrix[t].name;
      g2D.setColor(Color.black);
      g2D.drawString(myTruncate(name), 10, 10 * (t - myview.trait) + 36);
      for (int m = myview.marker; m < mydata.nmarkers; m++) {
        lod = mydata.qtlmatrix[t].scores[m];
        double y = 0;
        double cm = 0;
        double pcm = 0;

        if (myview.plotby == 0) {
          y = 150
              + myview.zoomlevel
              * ((m - myview.marker) + 25 * (mydata.markers[m].chromosome - mydata.markers[myview.marker].chromosome));
        }
        if (myview.plotby == 1) {
          cm = mydata.markers[m].location
              + mydata.chrlengths[mydata.markers[m].chromosome];
          pcm = mydata.markers[myview.marker].location
              + mydata.chrlengths[mydata.markers[myview.marker].chromosome];
          y = 150
              + myview.zoomlevel
              * ((cm - pcm) + 25 * ((mydata.markers[m].chromosome) - (mydata.markers[myview.marker].chromosome)));
        }
        g2D.setColor(colorFromLod(lod, mydata.maxqtl));
        if (mydata.modelmatrix != null) {
          if (mydata.modelmatrix[t].scores[m] == 1) {
            g2D.setColor(Color.blue);
          }
        }
        g2D.fillRect((int) y, 10 * (t - myview.trait) + 30,
            myview.zoomlevel * 1, 6);
      }
    }
  }

  void cistrans(Graphics2D g2D, Model mydata, View myview) {
    double lod;
    String name;
    for (int t = myview.trait; t < mydata.ntraits; t++) {
      name = mydata.qtlmatrix[t].name;
      g2D.setColor(Color.black);
      g2D.drawString(myTruncate(name), 10, 10 * (t - myview.trait) + 36);
      for (int m = myview.marker; m < mydata.nmarkers; m++) {
        lod = mydata.qtlmatrix[t].scores[m];
        g2D.setColor(colorFromLod(lod, mydata.maxqtl));
        if (mydata.modelmatrix != null) {
          if (mydata.modelmatrix[t].scores[m] == 1) {
            g2D.setColor(Color.blue);
          }
        }
        g2D.fillRect(150 + myview.zoomlevel * (m - myview.marker) + 10
            * mydata.markers[m].chromosome, 10 * (t - myview.trait) + 30,
            myview.zoomlevel * 1, 6);
      }
    }
  }

  void profile(Graphics2D g2D, Model mydata, View myview) {
    double lod;
    String name;
    int t = myview.trait;
    name = mydata.qtlmatrix[t].name;
    g2D.setColor(Color.black);
    g2D.setFont(new Font("sansserif", Font.BOLD, 14));
    g2D.drawString(name, 10, 10 * (t - myview.trait) + 36);
    for (int m = myview.marker; m < mydata.nmarkers; m++) {
      lod = mydata.qtlmatrix[t].scores[m];
      double y = 0;
      double cm = 0;
      double pcm = 0;
      if (myview.plotby == 0) {
        y = 15
            + myview.zoomlevel
            * (2 * (m - myview.marker) + 10 * (mydata.markers[m].chromosome - mydata.markers[myview.marker].chromosome));
      }
      if (myview.plotby == 1) {
        cm = mydata.markers[m].location
            + mydata.chrlengths[mydata.markers[m].chromosome];
        pcm = mydata.markers[myview.marker].location
            + mydata.chrlengths[mydata.markers[myview.marker].chromosome];
        y = 25
            + myview.zoomlevel
            * ((cm - pcm) + 10 * ((mydata.markers[m].chromosome) - (mydata.markers[myview.marker].chromosome)));
      }
      g2D.setColor(colorFromLod(lod, mydata.maxqtl));
      if (mydata.modelmatrix != null) {
        if (mydata.modelmatrix[t].scores[m] == 1) {
          g2D.setColor(Color.blue);
          g2D
              .drawString(
                  mydata.markers[m].name,
                  (int) y,
                  (int) (myview.getHeight() - ((lod * (myview.getHeight() / (1.1 * mydata.maxqtl))) + 45)));
        }
      }
      g2D.fillOval((int) y, (int) (myview.getHeight() - ((lod * (myview
          .getHeight() / (1.1 * mydata.maxqtl))) + 30)), myview.zoomlevel + 2,
          +myview.zoomlevel + 2);
    }
  }

  void circle(Graphics2D g2D, Model mydata, View myview) {
    double lod;
    String name;
    int t = myview.trait;
    name = mydata.qtlmatrix[t].name;
    g2D.setColor(Color.black);
    g2D.setFont(new Font("sansserif", Font.BOLD, 14));
    g2D.drawString(name, 10, 10 * (t - myview.trait) + 36);
    double maxy = 0;
    double y = 0;
    double cm = 0;
    double pcm = 0;
    if (myview.plotby == 0) {
      maxy = ((mydata.nmarkers - 0) + 25 * (1 + mydata.markers[mydata.nmarkers - 1].chromosome - mydata.markers[0].chromosome));
    }
    if (myview.plotby == 1) {
      cm = mydata.markers[mydata.nmarkers - 1].location
          + mydata.chrlengths[mydata.markers[mydata.nmarkers - 1].chromosome];
      pcm = mydata.markers[0].location
          + mydata.chrlengths[mydata.markers[0].chromosome];
      maxy = ((cm - pcm) + 25 * ((1 + mydata.markers[mydata.nmarkers - 1].chromosome) - (mydata.markers[0].chromosome)));
    }
    for (int m = 0; m < mydata.nmarkers; m++) {
      lod = mydata.qtlmatrix[t].scores[m];
      g2D.setColor(colorFromLod(lod, mydata.maxqtl));
      if (myview.plotby == 0) {
        y = ((m - myview.marker) + 25 * (mydata.markers[m].chromosome - mydata.markers[myview.marker].chromosome));
      }
      if (myview.plotby == 1) {
        cm = mydata.markers[m].location
            + mydata.chrlengths[mydata.markers[m].chromosome];
        pcm = mydata.markers[myview.marker].location
            + mydata.chrlengths[mydata.markers[myview.marker].chromosome];
        y = ((cm - pcm) + 25 * ((mydata.markers[m].chromosome) - (mydata.markers[myview.marker].chromosome)));
      }
      g2D.setColor(colorFromLod(lod, mydata.maxqtl));
      if (mydata.modelmatrix != null) {
        if (mydata.modelmatrix[t].scores[m] == 1) {
          g2D.setColor(Color.blue);
          g2D
              .drawString(mydata.markers[m].name, (int) (int) (myview
                  .getWidth() / 2 + ((myview.getWidth() - 100) / 3)
                  * Math.sin((2 * Math.PI) * (y / maxy))), (int) (myview
                  .getHeight() / 2 - ((myview.getHeight() - 100) / 3)
                  * Math.cos((2 * Math.PI) * (y / maxy))));
        }
      }
      g2D.fillOval((int) (myview.getWidth() / 2 + (myview.getWidth() / 3)
          * Math.sin((2 * Math.PI) * (y / maxy))),
          (int) (myview.getHeight() / 2 - (myview.getHeight() / 3)
              * Math.cos((2 * Math.PI) * (y / maxy))), (int) (5 + lod),
          (int) (5 + lod));
    }
  }
}
