#include "BC/BCenums.h"

std::map<std::string, BCcode> BCcodemap = boost::assign::map_list_of("periodic",periodic)("neumann",neumann)("dirichlet",dirichlet)("muiInterface",muiInterface)("mpi",mpi);

std::map<std::string, BCdir> BCdirmap  = boost::assign::map_list_of("north",north)("south",south)("east",east)("west",west)("top",top)("bottom",bottom);
