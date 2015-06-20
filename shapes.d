module poly2trid.shapes;

import std.math;

class Point 
{
	double x, y;
	
	/// The edges this point constitutes an upper ending point
	Edge[] edge_list;
	
	/// Construct using coordinates.
	this(double x, double y)
	{
		this.x = x;
		this.y = y;
	}
	
	/// Set this point to all zeros.
	void set_zero()
	{
		x = 0.0;
		y = 0.0;
	}
	
	/// Set this point to some specified coordinates.
	void set(double x_, double y_)
	{
		x = x_;
		y = y_;
	}

	/// Negate this point.
	Point opUnary(string s)() if (s == "-")
	{
		return new Point(-x, -y);
	}

	/// Add a point to this point.
	void opAddAssign(Point v)
	{
		x += v.x;
		y += v.y;
	}

	/// Subtract a point from this point.
	void opSubAssign(Point v)
	{
		x -= v.x;
		y -= v.y;
	}
	
	/// Multiply this point by a scalar.
	void opMulAssign(double a)
	{
		x *= a;
		y *= a;
	}

	/// Add two points_ component-wise.
	Point opAdd(Point a)
	{
		return new Point(x + a.x, y + a.y);
	}

	/// Subtract two points_ component-wise.
	Point opSub(Point a)
	{
		return new Point(x - a.x, y - a.y);
	}

	/// Multiply point by scalar
	Point opMul(Point a)
	{
		return new Point(x * a.x, y * a.y);
	}

	bool opEquals(U)(U other) pure const nothrow
		if (is(U : Point))
	{
		return a.x == x && a.y == y;
	}

	/// Peform the dot product on two vectors.
	double Dot(const Point a, const Point b)
	{
	 	return a.x * b.x + a.y * b.y;
	}
	
	/// Get the length of this point (the norm).
	double Length() const
	{
		return sqrt(x * x + y * y);
	}
	
	/// Convert this point into a unit point. Returns the Length.
	double Normalize()
	{
		double len = Length();
		x /= len;
		y /= len;
		return len;
	}
}

// Represents a simple polygon's edge
struct Edge {
	
	Point p, q;
	
	/// Constructor
	this(Point p1, Point p2)
	{
		p = p1;
		q = p2;
		if (p1.y > p2.y) {
			q = p1;
			p = p2;
		} else if (p1.y == p2.y) {
			if (p1.x > p2.x) {
				q = p1;
				p = p2;
			} else if (p1.x == p2.x) {
				// Repeat points
				assert(false, "Repeat points!");
			}
		}
		
		q.edge_list ~= this;
	}
}

/// Triangle-based data structures are know to have better performance than quad-edge structures
/// See: J. Shewchuk, "Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator"
///      "Triangulations in CGAL"
class Triangle 
{
	
	/// Constructor
	this(Point a, Point b, Point c)
	{
		points_ = [a, b, c];
		neighbors_ = [0: null, 1: null, 2: null];
		constrained_edge = [false, false, false];
		delaunay_edge = [false, false, false];
		interior_ = false;
	}
	
	/// Flags to determine if an edge is a Constrained edge
	bool[] constrained_edge;
	/// Flags to determine if an edge is a Delauney edge
	bool[] delaunay_edge;
	
	Point GetPoint(const int index)
	{
		return points_[index];
	}
	Point PointCW(Point point)
	{
		if (point == points_[0]) {
			return points_[2];
		} else if (point == points_[1]) {
			return points_[0];
		} else if (point == points_[2]) {
			return points_[1];
		}
		assert(0);
	}
	Point PointCCW(Point point)
	{
		if (point == points_[0]) {
			return points_[1];
		} else if (point == points_[1]) {
			return points_[2];
		} else if (point == points_[2]) {
			return points_[0];
		}
		assert(0);
	}
	Point OppositePoint(Triangle t, Point p)
	{
		Point cw = t.PointCW(p);
		double x = cw.x;
		double y = cw.y;
		x = p.x;
		y = p.y;
		return PointCW(cw);
	}
	
	Triangle GetNeighbor(const int index)
	{
		return neighbors_.get(index, null);
		//assert(0);
		//return null;
	}



	void MarkNeighbor(Point p1, Point p2, Triangle t)
	{
		if ((p1 == points_[2] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[2]))
    		neighbors_[0] = t;
		else if ((p1 == points_[0] && p2 == points_[2]) || (p1 == points_[2] && p2 == points_[0]))
			neighbors_[1] = t;
		else if ((p1 == points_[0] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[0]))
			neighbors_[2] = t;
		else
			assert(0);
	}
	void MarkNeighbor(Triangle t)
	{
		if (t.Contains(points_[1], points_[2])) {
			neighbors_[0] = t;
			t.MarkNeighbor(points_[1], points_[2], this);
		} else if (t.Contains(points_[0], points_[2])) {
			neighbors_[1] = t;
			t.MarkNeighbor(points_[0], points_[2], this);
		} else if (t.Contains(points_[0], points_[1])) {
			neighbors_[2] = t;
			t.MarkNeighbor(points_[0], points_[1], this);
		}
	}
	
	void MarkConstrainedEdge(const int index)
	{
		constrained_edge[index] = true;
	}
	void MarkConstrainedEdge(Edge edge)
	{
		MarkConstrainedEdge(edge.p, edge.q);
	}
	void MarkConstrainedEdge(Point p, Point q)
	{
		if ((q == points_[0] && p == points_[1]) || (q == points_[1] && p == points_[0])) {
			constrained_edge[2] = true;
		} else if ((q == points_[0] && p == points_[2]) || (q == points_[2] && p == points_[0])) {
			constrained_edge[1] = true;
		} else if ((q == points_[1] && p == points_[2]) || (q == points_[2] && p == points_[1])) {
			constrained_edge[0] = true;
		}
	}
	
	int Index(const Point p)
	{
		if (p == points_[0]) {
			return 0;
		} else if (p == points_[1]) {
			return 1;
		} else if (p == points_[2]) {
			return 2;
		}
		assert(0);
	}
	int EdgeIndex(const Point p1, const Point p2)
	{
		if (points_[0] == p1) {
			if (points_[1] == p2) {
				return 2;
			} else if (points_[2] == p2) {
				return 1;
			}
		} else if (points_[1] == p1) {
			if (points_[2] == p2) {
				return 0;
			} else if (points_[0] == p2) {
				return 2;
			}
		} else if (points_[2] == p1) {
			if (points_[0] == p2) {
				return 1;
			} else if (points_[1] == p2) {
				return 0;
			}
		}
		return -1;
	}
	
	Triangle NeighborCW(Point point)
	{
		if (point == points_[0]) {
			return neighbors_.get(1, null);
		} else if (point == points_[1]) {
			return neighbors_.get(2, null);
		}
		return neighbors_.get(0, null);
	}
	Triangle NeighborCCW(Point point)
	{
		if (point == points_[0]) {
			return neighbors_.get(2, null);
		} else if (point == points_[1]) {
			return neighbors_.get(0, null);
		}
		return neighbors_.get(1, null);
	}
	bool GetConstrainedEdgeCCW(Point p)
	{
		if (p == points_[0]) {
			return constrained_edge[2];
		} else if (p == points_[1]) {
			return constrained_edge[0];
		}
		return constrained_edge[1];
	}
	bool GetConstrainedEdgeCW(Point p)
	{
		if (p == points_[0]) {
			return constrained_edge[1];
		} else if (p == points_[1]) {
			return constrained_edge[2];
		}
		return constrained_edge[0];
	}
	void SetConstrainedEdgeCCW(Point p, bool ce)
	{
		if (p == points_[0]) {
			constrained_edge[2] = ce;
		} else if (p == points_[1]) {
			constrained_edge[0] = ce;
		} else {
			constrained_edge[1] = ce;
		}
	}
	void SetConstrainedEdgeCW(Point p, bool ce)
	{
		if (p == points_[0]) {
			constrained_edge[1] = ce;
		} else if (p == points_[1]) {
			constrained_edge[2] = ce;
		} else {
			constrained_edge[0] = ce;
		}
	}
	bool GetDelunayEdgeCCW(Point p)
	{
		if (p == points_[0]) {
			return delaunay_edge[2];
		} else if (p == points_[1]) {
			return delaunay_edge[0];
		}
		return delaunay_edge[1];
	}
	bool GetDelunayEdgeCW(Point p)
	{
		if (p == points_[0]) {
			return delaunay_edge[1];
		} else if (p == points_[1]) {
			return delaunay_edge[2];
		}
		return delaunay_edge[0];
	}
	void SetDelunayEdgeCCW(Point p, bool e)
	{
		if (p == points_[0]) {
			delaunay_edge[2] = e;
		} else if (p == points_[1]) {
			delaunay_edge[0] = e;
		} else {
			delaunay_edge[1] = e;
		}
	}
	void SetDelunayEdgeCW(Point p, bool e)
	{
		if (p == points_[0]) {
			delaunay_edge[1] = e;
		} else if (p == points_[1]) {
			delaunay_edge[2] = e;
		} else {
			delaunay_edge[0] = e;
		}
	}
	
	bool Contains(Point p)
	{
		return p == points_[0] || p == points_[1] || p == points_[2];
	}
	bool Contains(Edge e)
	{
		return Contains(e.p) && Contains(e.q);
	}
	bool Contains(Point p, Point q)
	{
		return Contains(p) && Contains(q);
	}
	void Legalize(Point point)
	{
		points_[1] = points_[0];
		points_[0] = points_[2];
		points_[2] = point;
	}
	void Legalize(Point opoint, Point npoint)
	{
		if (opoint == points_[0]) {
			points_[1] = points_[0];
			points_[0] = points_[2];
			points_[2] = npoint;
		} else if (opoint == points_[1]) {
			points_[2] = points_[1];
			points_[1] = points_[0];
			points_[0] = npoint;
		} else if (opoint == points_[2]) {
			points_[0] = points_[2];
			points_[2] = points_[1];
			points_[1] = npoint;
		} else {
			assert(0);
		}
	}
	/**
   * Clears all references to all other triangles and points
   */
	void Clear()
	{
		Triangle t;
	    for( int i=0; i<3; i++ )
	    {
	        t = neighbors_[i];
	        if( t !is null )
	        {
	            t.ClearNeighbor( this );
	        }
	    }
	    ClearNeighbors();
	    points_ = [];
	}
	void ClearNeighbor(Triangle triangle )
	{
		if( neighbors_[0] == triangle )
	    {
	        neighbors_[0] = null;
	    }
	    else if( neighbors_[1] == triangle )
	    {
	        neighbors_[1] = null;            
	    }
	    else
	    {
	        neighbors_[2] = null;
	    }
	}
	void ClearNeighbors()
	{
		Triangle[int] n;
		neighbors_ = n;
	}
	void ClearDelunayEdges()
	{
		delaunay_edge = [false, false, false];
	}
	
	bool IsInterior()
	{
		return interior_;
	}
	void IsInterior(bool b)
	{
		interior_ = b;
	}
	
	Triangle NeighborAcross(Point opoint)
	{
		if (opoint == points_[0]) {
			return neighbors_[0];
		} else if (opoint == points_[1]) {
			return neighbors_[1];
		}
		return neighbors_[2];
	}
	
	void DebugPrint()
	{
		import std.stdio, std.conv;
		write(text(points_[0].x, ",", points_[0].y, " "));
		write(text(points_[1].x, ",", points_[1].y, " "));
		write(text(points_[2].x, ",", points_[2].y));
		writeln();
	}
	
private:
	
	/// Triangle points
	Point[] points_;
	/// Neighbor list
	Triangle[int] neighbors_;
	
	/// Has this triangle been marked as an interior triangle?
	bool interior_;
}

public
{
	bool pt_cmp(const Point a, const Point b)
	{
		if (a.y < b.y) {
			return true;
		} else if (a.y == b.y) {
			// Make sure q is point with greater x value
			if (a.x < b.x) {
			  return true;
			}
		}
		return false;
	}
	
	/// Perform the cross product on two vectors. In 2D this produces a scalar.
	double Cross(const Point a, const Point b)
	{
	 	return a.x * b.y - a.y * b.x;
	}

	/// Perform the cross product on a point and a scalar. In 2D this produces
	/// a point.
	Point Cross(const Point a, double s)
	{
		return new Point(s * a.y, -s * a.x);
	}

	/// Perform the cross product on a scalar and a point. In 2D this produces
	/// a point.
	Point Cross(const double s, const Point a)
	{
	 	return new Point(-s * a.y, s * a.x);
	}
}