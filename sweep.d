module poly2trid.sweep;

import poly2trid;
import std.math;

public
{
	const double kAlpha = 0.3;
}

class SweepContext
{
	/// Constructor
	this(Point[] polyline)
	{
		basin = Basin();
		edge_event = EdgeEvent();
		points_ = polyline;
		InitEdges(points_);
	}

	void set_head(Point p1) 
	{ 
		head_ = p1; 
	}

	Point head() 
	{ 
		return head_; 
	}

	void set_tail(Point p1) 
	{ 
		tail_ = p1; 
	}

	Point tail() 
	{ 
		return tail_; 
	}

	int point_count() 
	{ 
		return points_.length; 
	}

	Node LocateNode(Point point)
	{
		return front_.LocateNode(point.x);
	}

	void CreateAdvancingFront(Node[] nodes)
	{
		// Initial triangle
		Triangle triangle = new Triangle(points_[0], tail_, head_);

		map_ ~= triangle;

		af_head_ = new Node(triangle.GetPoint(1), triangle);
		af_middle_ = new Node(triangle.GetPoint(0), triangle);
		af_tail_ = new Node(triangle.GetPoint(2));
		front_ = new AdvancingFront(af_head_, af_tail_);

		// TODO: More intuitive if head is middles next and not previous?
		//       so swap head and tail
		af_head_.next = af_middle_;
		af_middle_.next = af_tail_;
		af_middle_.prev = af_head_;
		af_tail_.prev = af_middle_;
	}

	/// Try to map a node to all sides of this triangle that don't have a neighbor
	void MapTriangleToNodes(Triangle t)
	{
		for (int i = 0; i < 3; i++) {
			if (!t.GetNeighbor(i)) {
				Node n = front_.LocatePoint(t.PointCW(t.GetPoint(i)));
				if (n)
					n.triangle = t;
			}
		}
	}

	void AddToMap(Triangle triangle)
	{
		map_ ~= triangle;
	}

	Point GetPoint(const int index)
	{
		return points_[index];
	}

	void RemoveFromMap(Triangle triangle)
	{
		for(int i = 0; i < map_.length; i++)
			if(map_[i] == triangle)
			{
				map_ = map_[0 .. i] ~ map_[i + 1 .. $];
				return;
			}
	}

	void AddHole(Point[] polyline)
	{
		InitEdges(polyline);
		points_ ~= polyline;
	}

	void AddPoint(Point point)
	{
		points_ ~= point;
	}

	AdvancingFront front() { return front_; }

	void MeshClean(Triangle triangle)
	{
		Triangle[] triangles = [triangle];

		while(triangles.length > 0){
			Triangle t = triangles[$ - 1];
			triangles = triangles[0 .. $ - 1];

			if (t !is null && !t.IsInterior()) {
				t.IsInterior(true);
				triangles_ ~= t;
				for (int i = 0; i < 3; i++) {
					if (!t.constrained_edge[i])
						triangles ~= t.GetNeighbor(i);
				}
			}
		}
	}

	Triangle[] GetTriangles()
	{
		return triangles_;
	}
	Triangle[] GetMap()
	{
		return map_;
	}

	Edge[] edge_list;

	struct Basin {
		Node left_node;
		Node bottom_node;
		Node right_node;
		double width;
		bool left_highest;

		void Clear()
		{
			left_node = null;
			bottom_node = null;
			right_node = null;
			width = 0.0;
			left_highest = false;
		}
	}

	struct EdgeEvent {
		Edge constrained_edge;
		bool right;
	}

	Basin basin;
	EdgeEvent edge_event;

private:
	Triangle[] triangles_;
	Triangle[] map_;
	Point[] points_;

	// Advancing front
	AdvancingFront front_;
	// head point used with advancing front
	Point head_;
	// tail point used with advancing front
	Point tail_;

	Node af_head_, af_middle_, af_tail_;

	void InitTriangulation()
	{
		double xmax = points_[0].x;
		double xmin = points_[0].x;
		double ymax = points_[0].y;
		double ymin = points_[0].y;

		// Calculate bounds.
		for (int i = 0; i < points_.length; i++) {
			Point p = points_[i];
			if (p.x > xmax)
				xmax = p.x;
			if (p.x < xmin)
				xmin = p.x;
			if (p.y > ymax)
				ymax = p.y;
			if (p.y < ymin)
				ymin = p.y;
		}

		double dx = kAlpha * (xmax - xmin);
		double dy = kAlpha * (ymax - ymin);
		head_ = new Point(xmax + dx, ymin - dy);
		tail_ = new Point(xmin - dx, ymin - dy);

		// Sort points along y-axis
		import std.algorithm;
		sort!(pt_cmp)(points_);
	}
	void InitEdges(Point[] polyline)
	{
		int num_points = polyline.length;
		for (int i = 0; i < num_points; i++) {
			int j = i < num_points - 1 ? i + 1 : 0;
			edge_list ~= Edge(polyline[i], polyline[j]);
		}
	}
}

/**
 * Sweep-line, Constrained Delauney Triangulation (CDT) See: Domiter, V. and
 * Zalik, B.(2008)'Sweep-line algorithm for constrained Delaunay triangulation',
 * International Journal of Geographical Information Science
 *
 * "FlipScan" Constrained Edge Algorithm invented by Thomas Åhlén, thahlen@gmail.com
 */

class Sweep
{
public:

	/**
	* Triangulate
	* 
	* @param tcx
	*/
	void Triangulate(SweepContext tcx)
	{
		tcx.InitTriangulation();
		tcx.CreateAdvancingFront(nodes_);
		// Sweep points; build mesh
		SweepPoints(tcx);
		// Clean up
		FinalizationPolygon(tcx);
	}

private:
	/**
	* Start sweeping the Y-sorted point set from bottom to top
	*/
	void SweepPoints(SweepContext tcx)
	{
		for (int i = 1; i < tcx.point_count(); i++) {
			Point point = tcx.GetPoint(i);
			Node node = PointEvent(tcx, point);
			for (int j = 0; j < point.edge_list.length; j++) {
				EdgeEvent(tcx, point.edge_list[j], node);
			}
		}
	}

	/**
	* Find closes node to the left of the new point and
	* create a new triangle. If needed new holes and basins
	* will be filled to.
	*/
	Node PointEvent(SweepContext tcx, Point point)
	{
		Node node = tcx.LocateNode(point);
		Node new_node = NewFrontTriangle(tcx, point, node);

		// Only need to check +epsilon since point never have smaller
		// x value than node due to how we fetch nodes from the front
		if (point.x <= node.point.x + EPSILON) {
			Fill(tcx, node);
		}
		FillAdvancingFront(tcx, new_node);
		return new_node;
	}

	void EdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		tcx.edge_event.constrained_edge = edge;
		tcx.edge_event.right = (edge.p.x > edge.q.x);

		if (IsEdgeSideOfTriangle(node.triangle, edge.p, edge.q)) {
			return;
		}

		// For now we will do all needed filling
		// TODO: integrate with flip process might give some better performance
		//       but for now this avoid the issue with cases that needs both flips and fills
		FillEdgeEvent(tcx, edge, node);
		EdgeEvent(tcx, edge.p, edge.q, node.triangle, edge.q);
	}

	void EdgeEvent(SweepContext tcx, Point ep, Point eq, Triangle triangle, Point point)
	{
		if (IsEdgeSideOfTriangle(triangle, ep, eq)) {
			return;
		}

		Point p1 = triangle.PointCCW(point);
		Orientation o1 = Orient2d(eq, p1, ep);
		if (o1 == Orientation.COLLINEAR) {
			if( triangle.Contains(eq, p1)) {
				triangle.MarkConstrainedEdge(eq, p1 );
				// We are modifying the constraint maybe it would be better to 
				// not change the given constraint and just keep a variable for the new constraint
				tcx.edge_event.constrained_edge.q = p1;
				triangle = triangle.NeighborAcross(point);
				EdgeEvent( tcx, ep, p1, triangle, p1 );
			} else {
				assert(0, "EdgeEvent - collinear points not supported");
			}
			return;
		}

		Point p2 = triangle.PointCW(point);
		Orientation o2 = Orient2d(eq, p2, ep);
		if (o2 == Orientation.COLLINEAR) {
			if( triangle.Contains(eq, p2)) {
				triangle.MarkConstrainedEdge(eq, p2 );
				// We are modifying the constraint maybe it would be better to 
				// not change the given constraint and just keep a variable for the new constraint
				tcx.edge_event.constrained_edge.q = p2;
				triangle = triangle.NeighborAcross(point);
				EdgeEvent( tcx, ep, p2, triangle, p2 );
			} else {
				assert(0, "EdgeEvent - collinear points not supported");
			}
			return;
		}

		if (o1 == o2) {
			// Need to decide if we are rotating CW or CCW to get to a triangle
			// that will cross edge
			if (o1 == Orientation.CW) {
				triangle = triangle.NeighborCCW(point);
			} else {
				triangle = triangle.NeighborCW(point);
			}
			EdgeEvent(tcx, ep, eq, triangle, point);
		} else {
			// This triangle crosses constraint so lets flippin start!
			FlipEdgeEvent(tcx, ep, eq, triangle, point);
		}
	}

	/**
	* Creates a new front triangle and legalize it
	*/
	Node NewFrontTriangle(SweepContext tcx, Point point, Node node)
	{
		Triangle triangle = new Triangle(point, node.point, node.next.point);

		triangle.MarkNeighbor(node.triangle);
		tcx.AddToMap(triangle);

		Node new_node = new Node(point);
		nodes_ ~= new_node;

		new_node.next = node.next;
		new_node.prev = node;
		node.next.prev = new_node;
		node.next = new_node;

		if (!Legalize(tcx, triangle)) {
			tcx.MapTriangleToNodes(triangle);
		}

		return new_node;
	}

	/**
	* Adds a triangle to the advancing front to fill a hole.
	* @param tcx
	* @param node - middle node, that is the bottom of the hole
	*/
	void Fill(SweepContext tcx, Node node)
	{
		Triangle triangle = new Triangle(node.prev.point, node.point, node.next.point);

		// TODO: should copy the constrained_edge value from neighbor triangles
		//       for now constrained_edge values are copied during the legalize
		triangle.MarkNeighbor(node.prev.triangle);
		triangle.MarkNeighbor(node.triangle);

		tcx.AddToMap(triangle);

		// Update the advancing front
		node.prev.next = node.next;
		node.next.prev = node.prev;

		// If it was legalized the triangle has already been mapped
		if (!Legalize(tcx, triangle)) {
			tcx.MapTriangleToNodes(triangle);
		}
	}

	/**
	* Returns true if triangle was legalized
	*/
	bool Legalize(SweepContext tcx, Triangle t)
	{
		// To legalize a triangle we start by finding if any of the three edges
		// violate the Delaunay condition
		for (int i = 0; i < 3; i++) {
			if (t.delaunay_edge[i])
				continue;

			Triangle ot = t.GetNeighbor(i);

			if (ot) {
				Point p = t.GetPoint(i);
				Point op = ot.OppositePoint(t, p);
				int oi = ot.Index(op);

				// If this is a Constrained Edge or a Delaunay Edge(only during recursive legalization)
				// then we should not try to legalize
				if (ot.constrained_edge[oi] || ot.delaunay_edge[oi]) {
					t.constrained_edge[i] = ot.constrained_edge[oi];
					continue;
				}

				bool inside = Incircle(p, t.PointCCW(p), t.PointCW(p), op);

				if (inside) {
					// Lets mark this shared edge as Delaunay
					t.delaunay_edge[i] = true;
					ot.delaunay_edge[oi] = true;

					// Lets rotate shared edge one vertex CW to legalize it
					RotateTrianglePair(t, p, ot, op);

					// We now got one valid Delaunay Edge shared by two triangles
					// This gives us 4 new edges to check for Delaunay

					// Make sure that triangle to node mapping is done only one time for a specific triangle
					bool not_legalized = !Legalize(tcx, t);
					if (not_legalized) {
						tcx.MapTriangleToNodes(t);
					}

					not_legalized = !Legalize(tcx, ot);
					if (not_legalized)
					tcx.MapTriangleToNodes(ot);

					// Reset the Delaunay edges, since they only are valid Delaunay edges
					// until we add a new triangle or point.
					// XXX: need to think about this. Can these edges be tried after we
					//      return to previous recursive level?
					t.delaunay_edge[i] = false;
					ot.delaunay_edge[oi] = false;

					// If triangle have been legalized no need to check the other edges since
					// the recursive legalization will handles those so we can end here.
					return true;
				}
			}
		}
		return false;
	}

	/**
	* Requirement:
	* 1. a,b and c form a triangle.
	* 2. a and d is know to be on opposite side of bc
	* 
	*                a
	*                +
	*               / \
	*              /   \
	*            b/     \c
	*            +-------+
	*           /    d    \
	*          /           \
	* 
	* Fact: d has to be in area B to have a chance to be inside the circle formed by
	*  a,b and c
	*  d is outside B if orient2d(a,b,d) or orient2d(c,a,d) is CW
	*  This preknowledge gives us a way to optimize the incircle test
	* @param a - triangle point, opposite d
	* @param b - triangle point
	* @param c - triangle point
	* @param d - point opposite a
	* @return true if d is inside circle, false if on circle edge
	*/
	bool Incircle(Point pa, Point pb, Point pc, Point pd)
	{
		double adx = pa.x - pd.x;
		double ady = pa.y - pd.y;
		double bdx = pb.x - pd.x;
		double bdy = pb.y - pd.y;

		double adxbdy = adx * bdy;
		double bdxady = bdx * ady;
		double oabd = adxbdy - bdxady;

		if (oabd <= 0)
			return false;

		double cdx = pc.x - pd.x;
		double cdy = pc.y - pd.y;

		double cdxady = cdx * ady;
		double adxcdy = adx * cdy;
		double ocad = cdxady - adxcdy;

		if (ocad <= 0)
			return false;

		double bdxcdy = bdx * cdy;
		double cdxbdy = cdx * bdy;

		double alift = adx * adx + ady * ady;
		double blift = bdx * bdx + bdy * bdy;
		double clift = cdx * cdx + cdy * cdy;

		double det = alift * (bdxcdy - cdxbdy) + blift * ocad + clift * oabd;

		return det > 0;
	}

	/**
	* Rotates a triangle pair one vertex CW
	*
	*       n2                    n2
	*  P +-----+             P +-----+
	*    | t  /|               |\  t |
	*    |   / |               | \   |
	*  n1|  /  |n3           n1|  \  |n3
	*    | /   |    after CW   |   \ |
	*    |/ oT |               | oT \|
	*    +-----+ oP            +-----+
	*       n4                    n4
	* 
	*/
	void RotateTrianglePair(Triangle t, Point p, Triangle ot, Point op)
	{
		Triangle n1, n2, n3, n4;
		n1 = t.NeighborCCW(p);
		n2 = t.NeighborCW(p);
		n3 = ot.NeighborCCW(op);
		n4 = ot.NeighborCW(op);

		bool ce1, ce2, ce3, ce4;
		ce1 = t.GetConstrainedEdgeCCW(p);
		ce2 = t.GetConstrainedEdgeCW(p);
		ce3 = ot.GetConstrainedEdgeCCW(op);
		ce4 = ot.GetConstrainedEdgeCW(op);

		bool de1, de2, de3, de4;
		de1 = t.GetDelunayEdgeCCW(p);
		de2 = t.GetDelunayEdgeCW(p);
		de3 = ot.GetDelunayEdgeCCW(op);
		de4 = ot.GetDelunayEdgeCW(op);

		t.Legalize(p, op);
		ot.Legalize(op, p);

		// Remap delaunay_edge
		ot.SetDelunayEdgeCCW(p, de1);
		t.SetDelunayEdgeCW(p, de2);
		t.SetDelunayEdgeCCW(op, de3);
		ot.SetDelunayEdgeCW(op, de4);

		// Remap constrained_edge
		ot.SetConstrainedEdgeCCW(p, ce1);
		t.SetConstrainedEdgeCW(p, ce2);
		t.SetConstrainedEdgeCCW(op, ce3);
		ot.SetConstrainedEdgeCW(op, ce4);

		// Remap neighbors
		// XXX: might optimize the markNeighbor by keeping track of
		//      what side should be assigned to what neighbor after the
		//      rotation. Now mark neighbor does lots of testing to find
		//      the right side.
		t.ClearNeighbors();
		ot.ClearNeighbors();
		if (n1) ot.MarkNeighbor(n1);
		if (n2) t.MarkNeighbor(n2);
		if (n3) t.MarkNeighbor(n3);
		if (n4) ot.MarkNeighbor(n4);
		t.MarkNeighbor(ot);
	}

	/**
	* Fills holes in the Advancing Front
	*/
	void FillAdvancingFront(SweepContext tcx, Node n)
	{
		// Fill right holes
		Node node = n.next;

		while (node.next) {
			// if HoleAngle exceeds 90 degrees then break.
			if (LargeHole_DontFill(node)) 
				break;
			Fill(tcx, node);
			node = node.next;
		}

		// Fill left holes
		node = n.prev;

		while (node.prev) {
			// if HoleAngle exceeds 90 degrees then break.
			if (LargeHole_DontFill(node)) 
				break;
			Fill(tcx, node);
			node = node.prev;
		}

		// Fill right basins
		if (n.next && n.next.next) {
			double angle = BasinAngle(n);
			if (angle < PI_3div4) {
				FillBasin(tcx, n);
			}
		}
	}

	// Decision-making about when to Fill hole. 
	// Contributed by ToolmakerSteve2
	bool LargeHole_DontFill(Node node)
	{
		Node nextNode = node.next;
		Node prevNode = node.prev;
		if (!AngleExceeds90Degrees(node.point, nextNode.point, prevNode.point))
			return false;

		// Check additional points on front.
		Node next2Node = nextNode.next;
		// "..Plus.." because only want angles on same side as point being added.
		if ((next2Node !is null) && !AngleExceedsPlus90DegreesOrIsNegative(node.point, next2Node.point, prevNode.point))
			return false;

		Node prev2Node = prevNode.prev;
		// "..Plus.." because only want angles on same side as point being added.
		if ((prev2Node !is null) && !AngleExceedsPlus90DegreesOrIsNegative(node.point, nextNode.point, prev2Node.point))
			return false;

		return true;
	}
	bool AngleExceeds90Degrees(Point origin, Point pa, Point pb)
	{
		double angle = Angle(origin, pa, pb);
		bool exceeds90Degrees = ((angle > PI_div2) || (angle < -PI_div2));
		return exceeds90Degrees;
	}
	bool AngleExceedsPlus90DegreesOrIsNegative(Point origin, Point pa, Point pb)
	{
		double angle = Angle(origin, pa, pb);
		bool exceedsPlus90DegreesOrIsNegative = (angle > PI_div2) || (angle < 0);
		return exceedsPlus90DegreesOrIsNegative;
	}
	double Angle(Point origin, Point pa, Point pb)
	{
		/* Complex plane
		* ab = cosA +i*sinA
		* ab = (ax + ay*i)(bx + by*i) = (ax*bx + ay*by) + i(ax*by-ay*bx)
		* atan2(y,x) computes the principal value of the argument function
		* applied to the complex number x+iy
		* Where x = ax*bx + ay*by
		*       y = ax*by - ay*bx
		*/
		double px = origin.x;
		double py = origin.y;
		double ax = pa.x- px;
		double ay = pa.y - py;
		double bx = pb.x - px;
		double by = pb.y - py;
		double x = ax * by - ay * bx;
		double y = ax * bx + ay * by;
		double angle = atan2(x, y);
		return angle;
	}

	/**
	*
	* @param node - middle node
	* @return the angle between 3 front nodes
	*/
	double HoleAngle(Node node)
	{
		/* Complex plane
		* ab = cosA +i*sinA
		* ab = (ax + ay*i)(bx + by*i) = (ax*bx + ay*by) + i(ax*by-ay*bx)
		* atan2(y,x) computes the principal value of the argument function
		* applied to the complex number x+iy
		* Where x = ax*bx + ay*by
		*       y = ax*by - ay*bx
		*/
		double ax = node.next.point.x - node.point.x;
		double ay = node.next.point.y - node.point.y;
		double bx = node.prev.point.x - node.point.x;
		double by = node.prev.point.y - node.point.y;
		return atan2(ax * by - ay * bx, ax * bx + ay * by);
	}

	/**
	* The basin angle is decided against the horizontal line [1,0]
	*/
	double BasinAngle(Node node)
	{
		double ax = node.point.x - node.next.next.point.x;
		double ay = node.point.y - node.next.next.point.y;
		return atan2(ay, ax);
	}

	/**
	* Fills a basin that has formed on the Advancing Front to the right
	* of given node.
	* First we decide a left,bottom and right node that forms the
	* boundaries of the basin. Then we do a reqursive fill.
	*
	* @param tcx
	* @param node - starting node, this or next node will be left node
	*/
	void FillBasin(SweepContext tcx, Node node)
	{
		if (Orient2d(node.point, node.next.point, node.next.next.point) == Orientation.CCW) {
			tcx.basin.left_node = node.next.next;
		} else {
			tcx.basin.left_node = node.next;
		}

		// Find the bottom and right node
		tcx.basin.bottom_node = tcx.basin.left_node;
		while (tcx.basin.bottom_node.next && tcx.basin.bottom_node.point.y >= tcx.basin.bottom_node.next.point.y) {
			tcx.basin.bottom_node = tcx.basin.bottom_node.next;
		}
		if (tcx.basin.bottom_node == tcx.basin.left_node) {
			// No valid basin
			return;
		}

		tcx.basin.right_node = tcx.basin.bottom_node;
		while (tcx.basin.right_node.next && tcx.basin.right_node.point.y < tcx.basin.right_node.next.point.y) {
			tcx.basin.right_node = tcx.basin.right_node.next;
		}
		if (tcx.basin.right_node == tcx.basin.bottom_node) {
			// No valid basins
			return;
		}

		tcx.basin.width = tcx.basin.right_node.point.x - tcx.basin.left_node.point.x;
		tcx.basin.left_highest = tcx.basin.left_node.point.y > tcx.basin.right_node.point.y;

		FillBasinReq(tcx, tcx.basin.bottom_node);
	}

	/**
	* Recursive algorithm to fill a Basin with triangles
	*
	* @param tcx
	* @param node - bottom_node
	* @param cnt - counter used to alternate on even and odd numbers
	*/
	void FillBasinReq(SweepContext tcx, Node node)
	{
		// if shallow stop filling
		if (IsShallow(tcx, node)) {
			return;
		}

		Fill(tcx, node);

		if (node.prev == tcx.basin.left_node && node.next == tcx.basin.right_node) {
			return;
		} else if (node.prev == tcx.basin.left_node) {
			Orientation o = Orient2d(node.point, node.next.point, node.next.next.point);
			if (o == Orientation.CW) {
				return;
			}
			node = node.next;
		} else if (node.next == tcx.basin.right_node) {
			Orientation o = Orient2d(node.point, node.prev.point, node.prev.prev.point);
			if (o == Orientation.CCW) {
				return;
			}
			node = node.prev;
		} else {
			// Continue with the neighbor node with lowest Y value
			if (node.prev.point.y < node.next.point.y) {
				node = node.prev;
			} else {
				node = node.next;
			}
		}

		FillBasinReq(tcx, node);
	}

	bool IsShallow(SweepContext tcx, Node node)
	{
		double height;

		if (tcx.basin.left_highest) {
			height = tcx.basin.left_node.point.y - node.point.y;
		} else {
			height = tcx.basin.right_node.point.y - node.point.y;
		}

		// if shallow stop filling
		if (tcx.basin.width > height) {
			return true;
		}
		return false;
	}

	bool IsEdgeSideOfTriangle(Triangle triangle, Point ep, Point eq)
	{
		int index = triangle.EdgeIndex(ep, eq);

		if (index != -1) {
			triangle.MarkConstrainedEdge(index);
			Triangle t = triangle.GetNeighbor(index);
			if (t) {
				t.MarkConstrainedEdge(ep, eq);
			}
			return true;
		}
		return false;
	}

	void FillEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		if (tcx.edge_event.right) {
			FillRightAboveEdgeEvent(tcx, edge, node);
		} else {
			FillLeftAboveEdgeEvent(tcx, edge, node);
		}
	}

	void FillRightAboveEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		while (node.next.point.x < edge.p.x) {
			// Check if next node is below the edge
			if (Orient2d(edge.q, node.next.point, edge.p) == Orientation.CCW) {
				FillRightBelowEdgeEvent(tcx, edge, node);
			} else {
				node = node.next;
			}
		}
	}

	void FillRightBelowEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		if (node.point.x < edge.p.x) {
			if (Orient2d(node.point, node.next.point, node.next.next.point) == Orientation.CCW) {
			// Concave
			FillRightConcaveEdgeEvent(tcx, edge, node);
			} else{
				// Convex
				FillRightConvexEdgeEvent(tcx, edge, node);
				// Retry this one
				FillRightBelowEdgeEvent(tcx, edge, node);
			}
		}
	}

	void FillRightConcaveEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		Fill(tcx, node.next);
		if (node.next.point != edge.p) {
			// Next above or below edge?
			if (Orient2d(edge.q, node.next.point, edge.p) == Orientation.CCW) {
				// Below
				if (Orient2d(node.point, node.next.point, node.next.next.point) == Orientation.CCW) {
					// Next is concave
					FillRightConcaveEdgeEvent(tcx, edge, node);
				} else {
					// Next is convex
				}
			}
		}
	}

	void FillRightConvexEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		// Next concave or convex?
		if (Orient2d(node.next.point, node.next.next.point, node.next.next.next.point) == Orientation.CCW) {
			// Concave
			FillRightConcaveEdgeEvent(tcx, edge, node.next);
		} else{
			// Convex
			// Next above or below edge?
			if (Orient2d(edge.q, node.next.next.point, edge.p) == Orientation.CCW) {
				// Below
				FillRightConvexEdgeEvent(tcx, edge, node.next);
			} else{
				// Above
			}
		}
	}

	void FillLeftAboveEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		while (node.prev.point.x > edge.p.x) {
			// Check if next node is below the edge
			if (Orient2d(edge.q, node.prev.point, edge.p) == Orientation.CW) {
				FillLeftBelowEdgeEvent(tcx, edge, node);
			} else {
				node = node.prev;
			}
		}
	}

	void FillLeftBelowEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		if (node.point.x > edge.p.x) {
			if (Orient2d(node.point, node.prev.point, node.prev.prev.point) == Orientation.CW) {
				// Concave
				FillLeftConcaveEdgeEvent(tcx, edge, node);
			} else {
				// Convex
				FillLeftConvexEdgeEvent(tcx, edge, node);
				// Retry this one
				FillLeftBelowEdgeEvent(tcx, edge, node);
			}
		}
	}

	void FillLeftConcaveEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		Fill(tcx, node.prev);
		if (node.prev.point != edge.p) {
			// Next above or below edge?
			if (Orient2d(edge.q, node.prev.point, edge.p) == Orientation.CW) {
				// Below
				if (Orient2d(node.point, node.prev.point, node.prev.prev.point) == Orientation.CW) {
					// Next is concave
					FillLeftConcaveEdgeEvent(tcx, edge, node);
				} else{
					// Next is convex
				}
			}
		}
	}

	void FillLeftConvexEdgeEvent(SweepContext tcx, Edge edge, Node node)
	{
		// Next concave or convex?
		if (Orient2d(node.prev.point, node.prev.prev.point, node.prev.prev.prev.point) == Orientation.CW) {
			// Concave
			FillLeftConcaveEdgeEvent(tcx, edge, node.prev);
		} else{
			// Convex
			// Next above or below edge?
			if (Orient2d(edge.q, node.prev.prev.point, edge.p) == Orientation.CW) {
				// Below
				FillLeftConvexEdgeEvent(tcx, edge, node.prev);
			} else{
				// Above
			}
		}
	}

	void FlipEdgeEvent(SweepContext tcx, Point ep, Point eq, Triangle t, Point p)
	{
		Triangle ot = t.NeighborAcross(p);
		Point op = ot.OppositePoint(t, p);

		if (ot is null) {
			// If we want to integrate the fillEdgeEvent do it here
			// With current implementation we should never get here
			//throw new RuntimeException( "[BUG:FIXME] FLIP failed due to missing triangle");
			assert(0);
		}

		if (InScanArea(p, t.PointCCW(p), t.PointCW(p), op)) {
			// Lets rotate shared edge one vertex CW
			RotateTrianglePair(t, p, ot, op);
			tcx.MapTriangleToNodes(t);
			tcx.MapTriangleToNodes(ot);

			if (p == eq && op == ep) {
				if (eq == tcx.edge_event.constrained_edge.q && ep == tcx.edge_event.constrained_edge.p) {
					t.MarkConstrainedEdge(ep, eq);
					ot.MarkConstrainedEdge(ep, eq);
					Legalize(tcx, t);
					Legalize(tcx, ot);
				} else {
					// XXX: I think one of the triangles should be legalized here?
				}
			} else {
				Orientation o = Orient2d(eq, op, ep);
				t = NextFlipTriangle(tcx, cast(int)o, t, ot, p, op);
				FlipEdgeEvent(tcx, ep, eq, t, p);
			}
		} else {
			Point newP = NextFlipPoint(ep, eq, ot, op);
			FlipScanEdgeEvent(tcx, ep, eq, t, ot, newP);
			EdgeEvent(tcx, ep, eq, t, p);
		}
	}

	/**
	* After a flip we have two triangles and know that only one will still be
	* intersecting the edge. So decide which to contiune with and legalize the other
	* 
	* @param tcx
	* @param o - should be the result of an orient2d( eq, op, ep )
	* @param t - triangle 1
	* @param ot - triangle 2
	* @param p - a point shared by both triangles 
	* @param op - another point shared by both triangles
	* @return returns the triangle still intersecting the edge
	*/
	Triangle NextFlipTriangle(SweepContext tcx, int o, Triangle  t, Triangle ot, Point p, Point op)
	{
		if (o == Orientation.CCW) {
			// ot is not crossing edge after flip
			int edge_index = ot.EdgeIndex(p, op);
			ot.delaunay_edge[edge_index] = true;
			Legalize(tcx, ot);
			ot.ClearDelunayEdges();
			return t;
		}

		// t is not crossing edge after flip
		int edge_index = t.EdgeIndex(p, op);

		t.delaunay_edge[edge_index] = true;
		Legalize(tcx, t);
		t.ClearDelunayEdges();
		return ot;
	}

	/**
	 * When we need to traverse from one triangle to the next we need 
	 * the point in current triangle that is the opposite point to the next
	 * triangle. 
	 * 
	 * @param ep
	 * @param eq
	 * @param ot
	 * @param op
	 * @return
	 */
	Point NextFlipPoint(Point ep, Point eq, Triangle ot, Point op)
	{
		Orientation o2d = Orient2d(eq, op, ep);
		if (o2d == Orientation.CW) {
			// Right
			return ot.PointCCW(op);
		} else if (o2d == Orientation.CCW) {
			// Left
			return ot.PointCW(op);
		} else{
			//throw new RuntimeException("[Unsupported] Opposing point on constrained edge");
			assert(0);
		}
	}

	/**
	 * Scan part of the FlipScan algorithm<br>
	 * When a triangle pair isn't flippable we will scan for the next 
	 * point that is inside the flip triangle scan area. When found 
	 * we generate a new flipEdgeEvent
	 * 
	 * @param tcx
	 * @param ep - last point on the edge we are traversing
	 * @param eq - first point on the edge we are traversing
	 * @param flipTriangle - the current triangle sharing the point eq with edge
	 * @param t
	 * @param p
	 */
	void FlipScanEdgeEvent(SweepContext tcx, Point ep, Point eq, Triangle flip_triangle, Triangle t, Point p)
	{
		Triangle ot = t.NeighborAcross(p);
		Point op = ot.OppositePoint(t, p);

		if (t.NeighborAcross(p) is null) {
			// If we want to integrate the fillEdgeEvent do it here
			// With current implementation we should never get here
			//throw new RuntimeException( "[BUG:FIXME] FLIP failed due to missing triangle");
			assert(0);
		}

		if (InScanArea(eq, flip_triangle.PointCCW(eq), flip_triangle.PointCW(eq), op)) {
			// flip with new edge op->eq
			FlipEdgeEvent(tcx, eq, op, ot, op);
			// TODO: Actually I just figured out that it should be possible to
			//       improve this by getting the next ot and op before the the above
			//       flip and continue the flipScanEdgeEvent here
			// set new ot and op here and loop back to inScanArea test
			// also need to set a new flip_triangle first
			// Turns out at first glance that this is somewhat complicated
			// so it will have to wait.
		} else{
			Point newP = NextFlipPoint(ep, eq, ot, op);
			FlipScanEdgeEvent(tcx, ep, eq, flip_triangle, ot, newP);
		}
	}

	void FinalizationPolygon(SweepContext tcx)
	{
		// Get an Internal triangle to start with
		Triangle t = tcx.front().head().next.triangle;
		Point p = tcx.front().head().next.point;
		while (t !is null && !t.GetConstrainedEdgeCW(p)) {
			t = t.NeighborCCW(p);
		}

		// Collect interior triangles constrained by edges
		tcx.MeshClean(t);
	}

	Node[] nodes_;
}