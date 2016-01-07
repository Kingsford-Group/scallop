#include "sgraph.h"

sgraph::sgraph(const bbase &b)
	:bundle(b)
{}

sgraph::~sgraph()
{}

int sgraph::solve()
{
	build();
	return 0;
}

int sgraph::build()
{
	// vertices: each region, start, end
	for(int i = 0; i < regions.size() + 2; i++) add_vertex(gr);

	// edges: each bridge
	for(int i = 0; i < bridges.size(); i++)
	{
		bridge &b = bridges[i];
		
		/*
		printf("add edge for bridge %d: %d -> %d\n", i, b.lrgn, b.rrgn);
		b.print();
		printf("\n");
		*/

		pair<edge_descriptor, bool> p = add_edge(b.lrgn, b.rrgn, gr);
		assert(p.second == true);
		MEI::iterator it = e2b.find(p.first);
		assert(it == e2b.end());
		e2b.insert(PEI(p.first, i));
	}

	// edges: connecting start/end and regions
	MPI mb;
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundary &b = boundaries[i];
		assert(mb.find(b.pos) == mb.end());
		mb.insert(PPI(b.pos, b.type));
	}

	int ss = regions.size();
	int tt = ss + 1;
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		if(mb.find(r.lpos) != mb.end())
		{
			int lt = mb[r.lpos];
			if(lt == LEFT_BOUNDARY) add_edge(ss, i, gr); 
			if(lt == START_BOUNDARY) add_edge(ss, i, gr);
			if(r.lpos < r.asc_pos) add_edge(ss, i, gr);
		}

		if(mb.find(r.rpos) != mb.end())
		{
			int rt = mb[r.rpos];
			if(rt == RIGHT_BOUNDARY) add_edge(i, tt, gr);
			if(rt == END_BOUNDARY) add_edge(i, tt, gr);
			if(r.rpos > r.desc_pos) add_edge(i, tt, gr);
		}
	}

	return 0;
}

int sgraph::print()
{
	printf("Bundle: ");
	printf("tid = %d, #hits = %lu, range = %s:%d-%d\n", tid, hits.size(), chrm.c_str(), lpos, rpos);
	// print hits
	/*
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].print();
	}
	*/

	// print bridges 
	for(int i = 0; i < bridges.size(); i++)
	{
		bridges[i].print();
	}

	// print boundaries
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundaries[i].print();
	}

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print();
	}

	printf("\n");
	return 0;
}
