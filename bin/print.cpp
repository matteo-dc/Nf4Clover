#include "aliases.hpp"
#include <fstream>
#include <iostream>

// print of file
void print_internal(double t, ofstream &outfile) { outfile << t << endl; }
void print_internal(VectorXd &V, ofstream &outfile)
{
  for (int i = 0; i < V.size(); i++)
    print_internal(V(i), outfile);
}

void print_internal_bin(VectorXd &V, ofstream &outfile)
{
  for (int i = 0; i < V.size(); i++)
    print_internal(V(i), outfile);
}
