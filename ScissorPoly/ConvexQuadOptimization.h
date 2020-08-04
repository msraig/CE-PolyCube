#ifndef CONVEX_QUAD_OPRIMIZATION_H
#define CONVEX_QUAD_OPRIMIZATION_H

#include "mosek.h"

#include <vector>

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
	printf("%s",str);
}; /* printstr */

bool solveConvexQuadPorgramming_mosek(std::vector<MSKboundkeye>& bkc, std::vector<double>& blc, std::vector<double>& buc,
	std::vector<MSKboundkeye>& bkx, std::vector<double>& blx, std::vector<double>& bux,
	std::vector<MSKlidxt>& aptrb, std::vector<MSKidxt>& asub, std::vector<double>& aval,
	std::vector<MSKidxt>& qsubi, std::vector<MSKidxt>& qsubj, std::vector<double>& qval, std::vector<double>& c,
	std::vector<double>& XX);


//initial guess might need to be added
bool solveConvexQuadPorgramming_mosek_integer(std::vector<MSKboundkeye>& bkc, std::vector<double>& blc, std::vector<double>& buc,
	std::vector<MSKboundkeye>& bkx, std::vector<double>& blx, std::vector<double>& bux,
	std::vector<MSKlidxt>& aptrb, std::vector<MSKidxt>& asub, std::vector<double>& aval,
	std::vector<MSKidxt>& qsubi, std::vector<MSKidxt>& qsubj, std::vector<double>& qval, std::vector<double>& c,
	std::vector<double>& XX);

bool mip_test();



#endif