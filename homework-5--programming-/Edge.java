public class Edge {

  public double distance;
  public Vertex source;
  public Vertex target;

  public Edge(Vertex vertex1, Vertex vertex2, double weight) {
    source = vertex1;
    target = vertex2;
    this.distance = weight;
  }

  public String toString() {
    return source + " - " + target;
  }
  
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null) {
      return false;
    }
    if (!(o instanceof Edge)) {
      return false;
    }
    Edge oEdge = (Edge) o;

    return (source.equals(oEdge.source) && target.equals(oEdge.target) 
        && distance == oEdge.distance) 
      || (source.equals(oEdge.target) && target.equals(oEdge.source) 
        && distance == oEdge.distance);
  }
}