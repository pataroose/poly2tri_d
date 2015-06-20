module poly2trid.cdt;

import poly2trid;

class CDT
{
public:
	/**
	* Constructor - add polyline with non repeating points
	* 
	* @param polyline
	*/
	this(Point[] polyline)
	{
		sweep_context_ = new SweepContext(polyline);
		sweep_ = new Sweep();
	}

	/**
	* Add a hole
	* 
	* @param polyline
	*/
	void AddHole(Point[] polyline)
	{
		sweep_context_.AddHole(polyline);
	}

	/**
	* Add a steiner point
	* 
	* @param point
	*/
	void AddPoint(Point point)
	{
		sweep_context_.AddPoint(point);
	}

	/**
	* Triangulate - do this AFTER you've added the polyline, holes, and Steiner points
	*/
	void Triangulate()
	{
		sweep_.Triangulate(sweep_context_);
	}

	/**
	* Get CDT triangles
	*/
	Triangle[] GetTriangles()
	{
		return sweep_context_.GetTriangles();
	}

	/**
	* Get triangle map
	*/
	Triangle[] GetMap()
	{
		return sweep_context_.GetMap();
	}

private:
	SweepContext sweep_context_;
	Sweep sweep_;
}