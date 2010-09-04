import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class QTLviewer {

  public JFrame window;
  public View myview;
  public Model mydata;
  private Controler mycontroler;
  private Menu mymenu;
  private WindowControls usercontrols;

  public QTLviewer() {
    window = new JFrame("QTLviewer - A Java(TM) Technology");
    mydata = new Model("./data.dat");
    myview = new View();
    myview.setSize(1024, 668);
    mycontroler = new Controler(window, myview, mydata);
    usercontrols = new WindowControls(mycontroler, window, myview, mydata);
    usercontrols.setSize(100, 100);
    window.addComponentListener(mycontroler);

    mymenu = new Menu(mycontroler);
    window.setJMenuBar(mymenu);
    window.getContentPane().setLayout(new BorderLayout());
    window.getContentPane().add(myview, BorderLayout.CENTER);
    window.getContentPane().add(usercontrols, BorderLayout.NORTH);
    window.getContentPane().add(new JPanel(), BorderLayout.WEST);
    window.getContentPane().add(usercontrols.leftslider(), BorderLayout.EAST);
    window.getContentPane().add(usercontrols.bottomslider(), BorderLayout.SOUTH);
    window.getContentPane().doLayout();
    window.pack();
    window.setSize(1024, 768);
    window.setBackground(Color.white);
    window.setVisible(true);

  }

  public static void main(String[] args) {
    new QTLviewer();
  }
}
