#include "sgraph.h"
#include "draw.h"

sgraph::sgraph(const bbase &b)
	:bundle(b)
{}

sgraph::~sgraph()
{}

int sgraph::solve()
{
	build();
	print();
	draw("sgraph.tex");
	//check();
	return 0;
}

int sgraph::build()
{
	// vertices: start, each region, end
	for(int i = 0; i < regions.size() + 2; i++) add_vertex(gr);

	// edges: each bridge
	for(int i = 0; i < bridges.size(); i++)
	{
		bridge &b = bridges[i];
		pair<edge_descriptor, bool> p = add_edge(b.lrgn + 1, b.rrgn + 1, gr);
		assert(p.second == true);
		MEI::iterator it = e2b.find(p.first);
		assert(it == e2b.end());
		e2b.insert(PEI(p.first, i));
	}

	// edges: connecting start/end and regions
	int ss = 0;
	int tt = regions.size() + 1;
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];

		if(r.lpos < r.asc_pos) add_edge(ss, i + 1, gr);
		else if(r.ltype == LEFT_BOUNDARY) add_edge(ss, i + 1, gr);
		else if(r.ltype == START_BOUNDARY) add_edge(ss, i + 1, gr);

		if(r.rpos > r.desc_pos) add_edge(i + 1, tt, gr);
		else if(r.rtype == RIGHT_BOUNDARY) add_edge(i + 1, tt, gr);
		else if(r.rtype == END_BOUNDARY) add_edge(i + 1, tt, gr);
	}

	return 0;
}

int sgraph::check()
{
	// check
	if(regions.size() == 1) return 0;
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		if(r.empty == false) continue;
		assert(r.empty == true);

		if(in_degree(i + 1, gr) != 0) printf("problem at region %d\n", i);
		assert(in_degree(i + 1, gr) == 0);

		if(out_degree(i + 1, gr) != 0) printf("problem at region %d\n", i);
		assert(out_degree(i + 1, gr) == 0);
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
		bridges[i].print(i);
	}

	// print boundaries
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundaries[i].print(i);
	}

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	printf("\n");
	return 0;
}

int sgraph::draw(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail())
	{
		printf("open file %s error.\n", file.c_str());
		return 0;
	}

	draw_header(fout);

	double len = 1.5;
	fout<<"\\def\\len{"<<len<<"cm}\n";

	// draw vertices
	char sx[1024];
	char sy[1024];
	for(int i = 0; i < num_vertices(gr); i++)
	{
		sprintf(sx, "s%d", i);
		double px = i * len;
		double py = 0.0;
		fout<<"\\node[mycircle, \\colx, draw] ("<<sx<<") at ("<<px<<", "<<py<<") {"<<i<<"};\n";
	}

	// draw edges
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		int s = source(*it1, gr);
		int t = target(*it1, gr);
		sprintf(sx, "s%d", s);
		sprintf(sy, "s%d", t);
		assert(s < t);
		
		double bend = 30;
		if(s + 1 == t) bend = 0;
		else if( (s + t) % 2 == 0 ) bend = -30;

		fout<<"\\draw[thick, ->, \\colx, bend right = "<< bend <<"] ("<<sx<<") to ("<<sy<<");\n";
	}

	draw_footer(fout);

	fout.close();
	return 0;
}

