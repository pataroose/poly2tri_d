module poly2trid.advancing_front;

import poly2trid;

class Node
{
	Point point;
	Triangle triangle;

	Node next;
	Node prev;

	double value;

	this(Point p)
	{
		point = p;
		value = p.x;
	}

	this(Point p, Triangle t)
	{
		point = p;
		triangle = t;
		value = p.x;
	}
}

class AdvancingFront 
{	
	this(Node head, Node tail)
	{
		head_ = head;
		tail_ = tail;
		search_node_ = head;
	}

	Node head()
	{
		return head_;
	}
	void set_head(Node node)
	{
		head_ = node;
	}
	Node tail()
	{
		return tail_;
	}
	void set_tail(Node node)
	{
		tail_ = node;
	}
	Node search()
	{
		return search_node_;
	}
	void set_search(Node node)
	{
		search_node_ = node;
	}

	/// Locate insertion point along advancing front
	Node LocateNode(const double x)
	{
		Node node = search_node_;

		if (x < node.value) {
			while ((node = node.prev) !is null) {
				if (x >= node.value) {
					search_node_ = node;
					return node;
				}
			}
		} else {
			while ((node = node.next) !is null) {
				if (x < node.value) {
					search_node_ = node.prev;
					return node.prev;
				}
			}
		}
		return null;
	}

	Node LocatePoint(const Point point)
	{
		const double px = point.x;
		Node node = FindSearchNode(px);
		const double nx = node.point.x;

		if (px == nx) {
			if (point != node.point) {
				// We might have two nodes with same x value for a short time
				if (point == node.prev.point) {
					node = node.prev;
				} else if (point == node.next.point) {
					node = node.next;
				} else {
					assert(0);
				}
			}
		} else if (px < nx) {
			while ((node = node.prev) !is null) {
				if (point == node.point) {
					break;
				}
			}
		} else {
			while ((node = node.next) !is null) {
				if (point == node.point)
					break;
			}
		}
		if(node) 
			search_node_ = node;
		return node;
	}

private:
	Node head_, tail_, search_node_;
	Node FindSearchNode(const double x)
	{
		return search_node_;
	}
}