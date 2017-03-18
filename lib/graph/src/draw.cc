/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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
	fout<<"\\def\\cola{blue}\n";
	fout<<"\\def\\colb{red}\n";
	fout<<"\\def\\colc{green}\n";
	fout<<"\\def\\cold{purple}\n";
	fout<<"\\def\\colx{black}\n";
	return 0;
}

int draw_footer(ofstream & fout)
{
	fout<<"\\end{tikzpicture}}\n";
	fout<<"\\end{document}\n";
	return 0;
}
