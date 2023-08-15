#ifndef BCENUMS_H
#define BCENUMS_H

#include <map>
#include "boost/assign.hpp"


enum BCcode {periodic, neumann, dirichlet, muiInterface, mpi}; 
enum BCdir {north=0, south, east, west, top, bottom};

extern std::map<std::string, BCcode> BCcodemap;
extern std::map<std::string, BCdir>  BCdirmap;

#endif
