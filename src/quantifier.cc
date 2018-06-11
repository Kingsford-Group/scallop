#include "quantifier.h"
#include <algorithm>

quantifier::quantifier(const vector<hit> &h, vector<path> &p)
	: hits(h), paths(p)
{}

int quantifier::quantify()
{
	if(paths.size() == 0) return 0;

	build_path_graph();
	build_super_hits();
	return 0;
}

int quantifier::build_path_graph()
{
	pgr.clear();
	if(paths.size() == 0) return 0;

	path &p = paths[0];
	assert(p.v.size() >= 2);
	assert(p.v[0] == 0);
	int n = p.v[p.v.size() - 1];

	for(int i = 0; i <= n; i++) pgr.add_vertex();

	for(int k = 0; k < paths.size(); k++)
	{
		p = paths[k];
		assert(p.v.size() >= 2);
		for(int i = 0; i < p.v.size() - 1; i++)
		{
			int s = p.v[i];
			int t = p.v[i + 1];
			assert(s >= 0 && s <= n);
			assert(t >= 0 && t <= n);
			pgr.add_edge(s, t);
		}
	}
	return 0;
}

int quantifier::build_super_hits()
{
	for(int k = 0; k < hits.size(); k++)
	{
		const set<int> &s = hits[k].phasing;
		
		if(check_valid_phasing(s) == false) continue;

		if(super_hits.find(s) == super_hits.end())
		{
			super_hit sh;
			sh.add_hit(k);
			super_hits.insert(PSSH(s, sh));
		}
		else
		{
			super_hits[s].add_hit(k);
		}
	}
	return 0;
}

bool quantifier::check_valid_phasing(const set<int> &s)
{
	if(s.size() <= 1) return true;

	vector<int> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	for(int k = 0; k < v.size() - 1; k++)
	{
		int s = v[k];
		int t = v[k + 1];
		PEB p = pgr.edge(s, t);
		if(p.second == false) return false;
	}
	return true;
}

int quantifier::print()
{
	int n = 0;
	if(paths.size() == 0) return 0;
	for(MSSH::iterator it = super_hits.begin(); it != super_hits.end(); it++)
	{
		n += it->second.hit_list.size();
	}
	printf("quantifier: total %lu paths, total %lu hits, total %lu super-hits, total %d valid phasing\n", paths.size(), hits.size(), super_hits.size(), n); 
	return 0;
}
