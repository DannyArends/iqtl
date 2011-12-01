import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

public class Menu extends JMenuBar {
  private static final long serialVersionUID = 1L;

  Menu(Controler mycntlr) {
    super();
    JMenu menu;
    JMenuItem item;

    menu = new JMenu("File");

    item = new JMenuItem("Open");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("Save");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("Exit");
    item.addActionListener(mycntlr);
    menu.add(item);

    add(menu);

    menu = new JMenu("View");

    item = new JMenuItem("Overview");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("Genotype plot");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("Heatmap plot");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("Cis/Trans plot");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("Circle plot");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("Profile plot");
    item.addActionListener(mycntlr);
    menu.add(item);

    add(menu);

    menu = new JMenu("Scale");

    item = new JMenuItem("Marker");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("CentiMorgan");
    item.addActionListener(mycntlr);
    menu.add(item);

    item = new JMenuItem("Basepair");
    item.addActionListener(mycntlr);
    menu.add(item);

    add(menu);
  }
}
