import java.util.PriorityQueue;

public class RandSep {
/*
 * Generate random numbers in [0,1] with guaranteed separation
 */
  private class Interval implements Comparable<Interval> {
    public double start;
    public double length;
    public Interval(double start,double length) {
      this.start = start;
      this.length = length;
    }
    double end() {return start + length; }
    public int compareTo(Interval o) {
      return (new Double(start).compareTo(new Double(o.start)));
    }
    void setStart(double new_start) {
      length = end() - new_start;
      start = new_start;
    }
    void setEnd(double new_end) {
      length = new_end - start;
    }
  }

  private PriorityQueue<Interval> intervals;
  private double sep;

  public RandSep(double sep) {
   this.sep = sep;
   intervals = new PriorityQueue<Interval>();
   intervals.add(new Interval(0,1));
  }
  
  private double total_length() {
    double sum = 0;
    for(Interval I: intervals)
      sum += I.length;
    return sum;
  }

  private Interval pickRandomInterval() {
    double spinner = total_length()*Math.random();
    double sum = 0;
    for(Interval I: intervals) {
      if((spinner > sum) && (spinner <= sum + I.length))
       return I;
      sum+=I.length;
     }
    return null;
  }

  private Interval getLastInterval() {
    /* Return, if any, an interval containing 1 */
    if(intervals.size()==0)
      return null;
    for(Interval I: intervals) {
      if(I.start + I.length==1)
        return I;
    }
    return null;
  }

  private Interval getFirstInterval() {
    /* Return, if any, an interval containing 0 */
    if(intervals.size()==0)
      return null;
    if(intervals.peek().start==0)
      return intervals.peek();
    else
      return null;

  }

  private void updateIntervals(double v,Interval I) {
    // Wrap beyond 0
    if(v<=sep) {
      // Remove or shrink left interval
      if(v+sep>I.start + I.length) intervals.remove(I);
      else I.setStart(v+sep);
      // Remove or shrink right interval
      Interval last = getLastInterval();
      if(last!=null) {
        if(last.start < 1+v-sep) last.setEnd(1+v-sep);
        else intervals.remove(last);
      }
      return;
    }
    // Wrap beyond 1
    if(v+sep>1) {
      // Remove or shrink right interval
      if(v-sep < I.start) intervals.remove(I);
      else I.setEnd(v-sep);
      Interval first = getFirstInterval();
      if(first!=null) {
        if(v+sep-1 < first.end()) {
          first.setStart(v+sep-1);
        } else intervals.remove(first);
      }
      return;
    }
    if(v-sep > I.start) {
      if(v+sep <= I.start + I.length) {
        // Delta_v in I
        intervals.remove(I);
        intervals.add(new Interval(I.start,v-sep-I.start));
        intervals.add(new Interval(v+sep,I.end()-v-sep));
      } 
      else I.setEnd(v - sep);
    }
    else {
      if(v+sep > I.start + I.length) intervals.remove(I);
      else I.setStart(v+sep);
    }
  }
  
  private void printIntervals() {
    
    for(Interval I: intervals) {
      System.out.format("\t\t(%.3g,%.3g)\n",I.start,I.end());
    }
  }
  
  public double next() {
    double v;
    if(intervals.size()==0) return -1;
    Interval I = pickRandomInterval();
    v = I.start + I.length*Math.random();
    updateIntervals(v,I);
    /* System.out.println("\tPrinting " + intervals.size() + " intervals");
    printIntervals(); */
    return v;
  }
}
