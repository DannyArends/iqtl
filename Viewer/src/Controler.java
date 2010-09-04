import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;

import javax.swing.JFrame;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class Controler extends ComponentAdapter implements ActionListener, ChangeListener {
  private ImageControler imageControl;
  private JFrame myapplication;
  private View myview;
  private Model mydata;
  private PlotFunctions myplots;

  public Controler(JFrame parent, View window, Model data) {
    System.out.println("Controler started");
    myapplication = parent;
    myview = window;
    mydata = data;
    myplots = new PlotFunctions();
  }

  public void componentResized(ComponentEvent e) {
    System.out.println("JFrame was resized");
    imageControl = new ImageControler(myplots, myview, mydata);
    new Thread(imageControl).start();

  }

  public void componentHidden(ComponentEvent e) {
    System.out.println("JFrame was closed");
    System.exit(0);
  }

  @Override
  public void actionPerformed(ActionEvent e) {
    System.out.println("Command: " + e.getActionCommand());
    if (e.getActionCommand().equals("Open")) {
      if (!mydata.loadnewModel(myapplication)) {
        mydata = mydata.previous;
        mydata.previous = null;
        System.out.println("Previous model saved, Load valid files ;)");
      }
    }
    if (e.getActionCommand().equals("Save")) {

    }
    if (e.getActionCommand().equals("Marker")) {
      myview.plotby = 0;
    }
    if (e.getActionCommand().equals("CentiMorgan")) {
      myview.plotby = 1;
    }
    if (e.getActionCommand().equals("Basepair")) {
      myview.plotby = 2;
    }
    if (e.getActionCommand().equals("Up")) {
      if (myview.trait > 0) {
        myview.trait--;
      }
    }
    if (e.getActionCommand().equals("Left")) {
      if (myview.marker > 0) {
        myview.marker--;
      }
    }
    if (e.getActionCommand().equals("Right")) {
      if (myview.marker < (mydata.nmarkers - 1)) {
        myview.marker++;
      }
    }
    if (e.getActionCommand().equals("Down")) {
      if (myview.trait < (mydata.ntraits - 1)) {
        myview.trait++;
      }
    }
    if (e.getActionCommand().equals("Overview")) {
      myview.setPlottype("Overview");
    }
    if (e.getActionCommand().contains("plot")) {
      myview.setPlottype(e.getActionCommand());
    }
    imageControl = new ImageControler(myplots, myview, mydata);
    new Thread(imageControl).start();
  }

  @Override
  public void stateChanged(ChangeEvent e) {
    JSlider slider = (JSlider) e.getSource();
    int value;
    if (!slider.getValueIsAdjusting()) {
      value = (int) slider.getValue();
      if (slider.getName().equals("rSlider")) {
        myview.trait = mydata.ntraits - value;
      }
      if (slider.getName().equals("bSlider")) {
        myview.marker = value;
      }
      if (slider.getName().equals("zSlider")) {
        myview.zoomlevel = value;
      }
      imageControl = new ImageControler(myplots, myview, mydata);
      new Thread(imageControl).start();
    }
  }
}