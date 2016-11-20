#include "hyper_set.h"
#include "config.h"
#include <algorithm>

int hyper_set::clear()
{
	nodes.clear();
	edges.clear();
	e2s.clear();
	ecnts.clear();
	return 0;
}

int hyper_set::add_node_list(const set<int> &s)
{
	return add_node_list(s, 1);
}

int hyper_set::add_node_list(const set<int> &s, int c)
{
	vector<int> v(s.begin(), s.end());
	return add_node_list(v, c);
}

int hyper_set::add_node_list(const vector<int> &s, int c)
{
	vector<int> v = s;
	sort(v.begin(), v.end());
	for(int i = 0; i < v.size(); i++) v[i]++;
	if(nodes.find(v) == nodes.end()) nodes.insert(PVII(v, c));
	else nodes[v] += c;
	return 0;
}

int hyper_set::build(directed_graph &gr, MEI& e2i)
{
	build_edges(gr, e2i);
	build_index();
	return 0;
}

int hyper_set::build_edges(directed_graph &gr, MEI& e2i)
{
	edges.clear();
	for(MVII::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		int c = it->second;
		if(c < min_router_count) continue;

		const vector<int> &vv = it->first;
		vector<int> ve;
		for(int k = 0; k < vv.size() - 1; k++)
		{
			PEB p = gr.edge(vv[k], vv[k + 1]);
			if(p.second == false) ve.push_back(-1);
			else ve.push_back(e2i[p.first]);
		}
		edges.push_back(ve);
		ecnts.push_back(c);
	}

	return 0;
}

int hyper_set::build_index()
{
	e2s.clear();
	for(int i = 0; i < edges.size(); i++)
	{
		vector<int> &v = edges[i];
		for(int j = 0; j < v.size(); j++)
		{
			int e = v[j];
			if(e == -1) continue;
			if(e2s.find(e) == e2s.end())
			{
				set<int> s;
				s.insert(i);
				e2s.insert(PISI(e, s));
			}
			else
			{
				e2s[e].insert(i);
			}
		}
	}
	return 0;
}

set<int> hyper_set::get_intersection(const vector<int> &v)
{
	set<int> ss;
	if(v.size() == 0) return ss;
	assert(v[0] >= 0);
	if(e2s.find(v[0]) == e2s.end()) return ss;
	ss = e2s[v[0]];
	vector<int> vv(ss.size());
	for(int i = 1; i < v.size(); i++)
	{
		assert(v[i] >= 0);
		set<int> s;
		if(e2s.find(v[i]) == e2s.end()) return s;
		s = e2s[v[i]];
		vector<int>::iterator it = set_intersection(ss.begin(), ss.end(), s.begin(), s.end(), vv.begin());
		ss = set<int>(vv.begin(), it);
	}
	return ss;
}

set<int> hyper_set::get_successors(int e)
{
	set<int> s;
	if(e2s.find(e) == e2s.end()) return s;
	set<int> &ss = e2s[e];
	for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
	{
		vector<int> &v = edges[*it];
		for(int i = 0; i < v.size(); i++)
		{
			if(v[i] != e) continue;
			if(i >= v.size() - 1) continue;
			if(v[i + 1] == -1) continue;
			s.insert(v[i + 1]);
		}
	}
	return s;
}

set<int> hyper_set::get_predecessors(int e)
{
	set<int> s;
	if(e2s.find(e) == e2s.end()) return s;
	set<int> &ss = e2s[e];
	for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
	{
		vector<int> &v = edges[*it];
		for(int i = 0; i < v.size(); i++)
		{
			if(v[i] != e) continue;
			if(i == 0) continue;
			if(v[i - 1] == -1) continue;
			s.insert(v[i - 1]);
		}
	}
	return s;
}

vector<PI> hyper_set::get_routes(int x, directed_graph &gr, MEI &e2i)
{
	edge_iterator it1, it2;
	vector<PI> v;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		assert(e2i.find(*it1) != e2i.end());
		int e = e2i[*it1];
		set<int> s = get_successors(e);
		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			PI p(e, *it);
			v.push_back(p);
		}
	}
	return v;
}

int hyper_set::get_routes(int x, directed_graph &gr, MEI &e2i, MPII &mpi)
{
	edge_iterator it1, it2;
	mpi.clear();
	int total = 0;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		assert(e2i.find(*it1) != e2i.end());
		int e = e2i[*it1];

		if(e2s.find(e) == e2s.end()) continue;
		set<int> &ss = e2s[e];
		for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
		{
			int k = *it;
			assert(k >= 0 && k < edges.size());
			assert(k >= 0 && k < ecnts.size());
			vector<int> &v = edges[k];
			int cnt = ecnts[k];
			for(int i = 0; i < v.size(); i++)
			{
				if(v[i] != e) continue;
				if(i == v.size() - 1) continue;
				if(v[i + 1] == -1) continue;
				PI p(e, v[i + 1]);
				total += cnt;
				if(mpi.find(p) != mpi.end()) mpi[p] += cnt;
				else mpi.insert(PPII(p, cnt));
			}
		}
	}
	return total;
}

int hyper_set::replace(int x, int e)
{
	vector<int> v;
	v.push_back(x);
	replace(v, e);
	return 0;
}

int hyper_set::replace(int x, int y, int e)
{
	vector<int> v;
	v.push_back(x);
	v.push_back(y);
	replace(v, e);
	return 0;
}

int hyper_set::replace(const vector<int> &v, int e)
{
	if(v.size() == 0) return 0;
	set<int> s = get_intersection(v);
	
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		int b = consecutive_subset(vv, v);
		assert(b >= 0);
		vv[b] = e;
		vv.erase(vv.begin() + b + 1, vv.begin() + b + v.size());

		for(int i = 0; i < v.size(); i++) 
		{
			int x = v[i];
			assert(e2s.find(x) != e2s.end());
			assert(e2s[x].find(k) != e2s[x].end());
			e2s[x].erase(k);
		}

		if(e2s.find(e) == e2s.end())
		{
			set<int> ss;
			ss.insert(k);
			e2s.insert(PISI(e, ss));
		}
		else
		{
			e2s[e].insert(k);
		}
	}

	return 0;
}

int hyper_set::remove(const set<int> &s)
{
	return remove(vector<int>(s.begin(), s.end()));
}

int hyper_set::remove(const vector<int> &v)
{
	for(int i = 0; i < v.size(); i++) remove(v[i]);
	return 0;
}

int hyper_set::remove(int e)
{
	if(e2s.find(e) == e2s.end()) return 0;
	set<int> s = e2s[e];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 0; i < vv.size(); i++)
		{
			if(vv[i] == e) vv[i] = -1;
		}

		/*
		if(vv[0] == e) vv.erase(vv.begin());
		else if(vv[vv.size() - 1] == e) vv.pop_back();
		else assert(false);
		*/
	}
	e2s.erase(e);
	return 0;
}

int hyper_set::remove_pair(int x, int y)
{
	if(e2s.find(x) == e2s.end()) return 0;
	set<int> s = e2s[x];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 0; i < vv.size(); i++)
		{
			if(i == vv.size() - 1) continue;
			if(vv[i] != x) continue;
			if(vv[i + 1] != y) continue;
			vv.insert(vv.begin() + i + 1, -1);
			break;
		}
	}
	return 0;
}

bool hyper_set::extend(int e)
{
	return (left_extend(e) || right_extend(e));
}

bool hyper_set::left_extend(int e)
{
	if(e2s.find(e) == e2s.end()) return false;
	set<int> s = e2s[e];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 1; i < vv.size(); i++)
		{
			if(vv[i] == e && vv[i - 1] != -1) return true; 
		}
	}
	return false;
}

bool hyper_set::right_extend(int e)
{
	if(e2s.find(e) == e2s.end()) return false;
	set<int> s = e2s[e];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 0; i < vv.size() - 1; i++)
		{
			if(vv[i] == e && vv[i + 1] != -1) return true; 
		}
	}
	return false;
}

bool hyper_set::left_extend(const set<int> &s)
{
	for(set<int>::const_iterator it = s.begin(); it != s.end(); it++)
	{
		if(left_extend(*it) == true) return true;
	}
	return false;
}

bool hyper_set::right_extend(const set<int> &s)
{
	for(set<int>::const_iterator it = s.begin(); it != s.end(); it++)
	{
		if(right_extend(*it) == true) return true;
	}
	return false;
}

int hyper_set::print()
{
	//printf("PRINT HYPER_SET\n");
	/*
	for(MVII::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		int c = it->second;
		printf("hyper-edge (nodes), counts = %d, list = ( ", c); 
		printv(v);
		printf(")\n");
	}
	*/

	for(int i = 0; i < edges.size(); i++)
	{
		printf("hyper-edge (edges) %d: ( ", i);
		printv(edges[i]);
		printf(")\n");
	}
	return 0;
}
