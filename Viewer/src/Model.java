import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.swing.JFileChooser;
import javax.swing.JFrame;

class marker {
  public String name;
  public int chromosome;
  public double location;
}

class qtl {
  public String name;
  public double[] scores;
}

class qtlmodel {
  public String name;
  public int[] scores;
}

public class Model {
  public int nmarkers, nindividuals, nchromosomes, ntraits, maxqtl, minqtl;
  public int[] chrlengths;
  public qtl[] qtlmatrix;
  public qtlmodel[] modelmatrix;
  public marker[] markers;
  public marker[] traits;
  public Model previous;

  public Model(String path) {
    System.out.println("MODEL: started parsing from string: " + path);
    File myfile = new File(path);
    readmodeldata(myfile);
    System.out.println("MODEL: Parsing done");
  }

  public boolean loadnewModel(JFrame frame) {
    previous = this;
    System.out.println("MODEL: opening new model from file");
    JFileChooser chooser = new JFileChooser();
    chooser.showOpenDialog(frame);
    File file = chooser.getSelectedFile();
    if (file != null && file.exists() && file.isFile()) {
      readmodeldata(file);
      System.out.println("MODEL: Parsing done");
      return true;
    } else {
      return false;
    }
  }

  public void readmodeldata(File myfile) {
    modelmatrix = null;
    traits = null;
    try {
      BufferedReader input = new BufferedReader(new FileReader(myfile));
      try {
        String line = null; // entire line
        String param = null; // param ( Before = )
        String value = null; // value ( After = )
        String[] values; // split of values with ,
        while ((line = input.readLine()) != null) {
          try {
            param = line.split("=")[0];
            value = line.split("=")[1];
          } catch (Exception e) {
            System.out.println("Unparsable line: " + line);
          }
          if (param.equals("markers")) {
            nmarkers = Integer.parseInt(value);
            markers = new marker[nmarkers];
          }
          if (param.equals("individuals")) {
            nindividuals = Integer.parseInt(value);
            markers = new marker[nmarkers];
          }
          if (param.equals("traits")) {
            ntraits = Integer.parseInt(value);
            qtlmatrix = new qtl[ntraits];
          }
          if (param.equals("locdata")) {
            if (Integer.parseInt(value) == 1) {
              System.out.println("Traits with genetic locations");
              traits = new marker[ntraits];
            }
          }
          if (param.equals("modeldata")) {
            if (Integer.parseInt(value) == 1) {
              System.out.println("MQM Model present");
              modelmatrix = new qtlmodel[ntraits];
            }
          }
          if (param.equals("chromosomes")) {
            nchromosomes = Integer.parseInt(value);
            chrlengths = new int[nchromosomes + 1];
          }
          if (param.equals("maxQTL")) {
            maxqtl = Integer.parseInt(value);
          }
          if (param.equals("minQTL")) {
            minqtl = Integer.parseInt(value);
          }
          if (param.equals("length")) {
            values = value.split(",");
            for (int x = 0; x < values.length; x++) {
              chrlengths[x] = Integer.parseInt(values[x]);
            }
          }
          try {
            String[] split = param.split("\\.");
            if (split[0].equals("qtldata")) {
              qtl n = new qtl();
              values = value.split(",");
              n.name = values[0];
              n.scores = new double[nmarkers];
              for (int x = 1; x < values.length; x++) {
                n.scores[x - 1] = Double.parseDouble(values[x]);
              }
              qtlmatrix[Integer.parseInt(split[1])] = n;
            }
            if (split[0].equals("model")) {
              qtlmodel n = new qtlmodel();
              values = value.split(",");
              n.name = qtlmatrix[(Integer.parseInt(split[1]))].name;
              n.scores = new int[nmarkers];
              for (int x = 0; x < values.length; x++) {
                n.scores[x] = Integer.parseInt(values[x]);
              }
              modelmatrix[Integer.parseInt(split[1])] = n;
            }
            if (split[0].equals("loc")) {
              marker n = new marker();
              values = value.split(",");
              n.name = values[0];
              n.chromosome = (Integer.parseInt(values[1]) - 1);
              n.location = Double.parseDouble(values[2]);
              traits[Integer.parseInt(split[1])] = n;
            }
            if (split[0].equals("mapdata")) {
              marker n = new marker();
              values = value.split(",");
              n.name = values[0];
              n.chromosome = (Integer.parseInt(values[1]) - 1);
              n.location = Double.parseDouble(values[2]);
              markers[Integer.parseInt(split[1])] = n;
            }
          } catch (Exception e) {
            System.out.println("Line: " + line);
          }
        }
      } finally {
        input.close();
      }
    } catch (IOException ex) {
      ex.printStackTrace();
    }

  }
}
