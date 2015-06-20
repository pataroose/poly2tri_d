module poly2trid.utils;

import poly2trid.shapes;

import std.math;

public
{
	const double PI_3div4 = 3 * PI_4;
	const double PI_div2 = 1.57079632679489661923;
	const double EPSILON = 1e-12;

	enum Orientation { CW, CCW, COLLINEAR };

	/**
	* Forumla to calculate signed area
	* Positive if CCW
	* Negative if CW
	* 0 if collinear
	*/
	Orientation Orient2d(Point pa, Point pb, Point pc)
	{
		double detleft = (pa.x - pc.x) * (pb.y - pc.y);
		double detright = (pa.y - pc.y) * (pb.x - pc.x);
		double val = detleft - detright;
		if (val > -EPSILON && val < EPSILON) {
			return Orientation.COLLINEAR;
		} else if (val > 0) {
			return Orientation.CCW;
		}
		return Orientation.CW;
	}

	bool InScanArea(Point pa, Point pb, Point pc, Point pd)
	{
		double oadb = (pa.x - pb.x)*(pd.y - pb.y) - (pd.x - pb.x)*(pa.y - pb.y);
		if (oadb >= -EPSILON) {
			return false;
		}

		double oadc = (pa.x - pc.x)*(pd.y - pc.y) - (pd.x - pc.x)*(pa.y - pc.y);
		if (oadc <= EPSILON) {
			return false;
		}
		return true;
	}
}