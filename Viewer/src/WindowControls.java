/**
 * \file WindowControls.java
 * \brief Code file containing the controls class
 *
 * last modified Sep, 2010
 * first written Apr, 2010
 * Copyright (c) 2010 Danny Arends
 * 
 **/
 
import java.awt.BorderLayout;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;

public class WindowControls extends JPanel {
  private static final long serialVersionUID = 1L;
  private Model mydata;
  private Controler mycontroler;
  private JPanel mypanel;

  WindowControls(Controler controler, JFrame parent, View view, Model data) {
    super();
    mydata = data;
    JButton n;
    mycontroler = controler;
    mypanel = new JPanel();
    mypanel.setLayout(new BorderLayout());
    n = new JButton("Up");
    n.addActionListener(mycontroler);
    mypanel.add(n, BorderLayout.NORTH);
    n = new JButton("Left");
    n.addActionListener(mycontroler);
    mypanel.add(n, BorderLayout.WEST);
    n = new JButton("Right");
    n.addActionListener(mycontroler);
    mypanel.add(n, BorderLayout.EAST);
    n = new JButton("Down");
    n.addActionListener(mycontroler);
    mypanel.add(n, BorderLayout.SOUTH);
    add(new JLabel("Move"), BorderLayout.CENTER);
    mypanel.doLayout();

    JSlider nslide = new JSlider(JSlider.HORIZONTAL, 1, 25, 1);
    nslide.setName("zSlider");
    nslide.addChangeListener(mycontroler);

    setLayout(new BorderLayout());
    add(mypanel, BorderLayout.EAST);

    add(nslide, BorderLayout.SOUTH);
    doLayout();
  }

  JSlider leftslider() {
    JSlider nslide = new JSlider(JSlider.VERTICAL, 0, mydata.ntraits,
        mydata.ntraits);
    nslide.setName("rSlider");
    nslide.addChangeListener(mycontroler);
    return nslide;
  }

  JSlider bottomslider() {
    JSlider nslide = new JSlider(JSlider.HORIZONTAL, mydata.nmarkers, 0);
    nslide.setName("bSlider");
    nslide.addChangeListener(mycontroler);
    return nslide;
  }
}
