import java.util.ArrayList;
import java.util.Collections;

public class RandSepTest {
  public static void main(String args[]) {
    double r;
    ArrayList<Double> rs = new ArrayList<Double>();
    RandSep randSep = new RandSep(0.1);
    for(int i = 0; i < 8; i++) {
      r = randSep.next();
      System.out.format("%.3g\n",r);
      rs.add(new Double(r));
    }
    while(rs.remove(new Double(-1)));
    System.out.println("\nSorted List:");
    Collections.sort(rs);
    for(Double value:rs) {
      System.out.format("%.3g\n",value);
    }
    
    double min_deviation = 1;
    double previous = rs.get(rs.size()-1)-1;
    for(Double value:rs) {
      min_deviation = Math.min(min_deviation,value-previous);
      previous = value;
    }
    
    System.out.format("\nThe minimum deviation is %.3g\n",min_deviation);
  }
}
