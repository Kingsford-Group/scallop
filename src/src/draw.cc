#include "draw.h"

int draw_header(ofstream & fout)
{
	fout<<"\\documentclass{llncs}\n";
	fout<<"\\usepackage{tikz}	\n";
	fout<<"\\usetikzlibrary{calc}\n";
	fout<<"\\usetikzlibrary{shapes.geometric}\n";
	fout<<"\\usetikzlibrary{fit}\n";
	fout<<"\\begin{document}\n";
	fout<<"{\\begin{tikzpicture}[mycircle/.style={draw, circle, minimum size=1.1em, inner sep = 0mm}, >=stealth]\n";
	fout<<"\\def\\cola{red}\n";
	fout<<"\\def\\colb{green}\n";
	fout<<"\\def\\colc{blue}\n";
	fout<<"\\def\\cold{gray}\n";
	fout<<"\\def\\colx{black}\n";
	return 0;
}

int draw_footer(ofstream & fout)
{
	fout<<"\\end{tikzpicture}}\n";
	fout<<"\\end{document}\n";
	return 0;
}
