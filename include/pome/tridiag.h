// Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
// LICENCE available at https://www.gnu.org/software/gsl/doc/html/gpl.html
#ifndef TRIDIAG_H_
#define TRIDIAG_H_

#include <gsl/gsl_errno.h>

#include <cstdlib>
#include <iostream>
#include <vector>

namespace GSL {

int solve_tridiag_nonsym(const std::vector<double>& diag, const std::vector<double>& abovediag,
                         const std::vector<double>& belowdiag, const std::vector<double>& rhs, std::vector<double>& x,
                         size_t N);

int gsl_linalg_solve_tridiag(const std::vector<double>& diag, const std::vector<double>& abovediag,
                             const std::vector<double>& belowdiag, const std::vector<double>& rhs,
                             std::vector<double>& solution);

}  // namespace GSL

#endif  // TRIDIAG_H_