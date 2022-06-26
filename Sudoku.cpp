#include "Highs.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <string>
//#include <vector>

using std::cin;
using std::cout;
using std::endl;

const int square_dim = 3;
const int grid_dim = square_dim * square_dim;
const int num_cell = grid_dim * grid_dim;
const int num_var = num_cell * grid_dim;

const double inf = kHighsInf;

int varIndex(const int i, const int j, const int k);
HighsStatus defineLp(Highs& highs);

int main(int argc, char *argv[]) {
  // Check that there is a file specified and open it as std::ifstream
  assert(argc == 2);
  std::ifstream f(argv[1]);

  // Instantiate HiGHS
  Highs highs;
  // Define the LP variables and constraints
  HighsStatus return_status = defineLp(highs);
  assert(return_status == HighsStatus::kOk);

  // Set up const reference to internal info structure and LP
  const HighsInfo& info = highs.getInfo();
  const HighsLp& lp = highs.getLp();
  const int num_row = lp.num_row_;

  // Set up vectors of costs, bounds and RHS to be re-used, with upper being constant
  std::vector<double> cost;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> rhs;

  int problem_num = 0;
  // Read a sequence of Sudoku grids
  std::vector<int> grid;
  grid.assign(num_cell, 0);
  for (;;) {
    const bool og_read = false;
    std::string sudoku;
    if (og_read) {
      std::stringstream ss;
      ss << f.rdbuf();
      sudoku = std::move(ss.str());
    } else {
      std::getline(f, sudoku);
      printf("Read Sudoku line %s\n", sudoku.c_str());
    }

    int size = sudoku.length();
    if (size < num_cell) {
      if (size) printf("STRANGE: Line has size %d\n", size);
      break;
    } 
    
    char c;
    for (int k = 0; k < num_cell; k++) {
      c = sudoku[k];
      if (c != '.') {
	assert(isdigit(c));
	grid[k] = c - '0';
      }
      //    printf("%2d: %1d from %1c\n", k, grid[k], sudoku[k]);
    }
    // Set up costs
    cost.assign(num_var, 0);
    return_status = highs.changeColsCost(0, num_var-1, &cost[0]);
    assert(return_status == HighsStatus::kOk);
    // Set up bounds
    lower.assign(num_var, 0);
    upper.assign(num_var, 1);
    for (int i=0; i<grid_dim; i++)
      for (int j=0; j<grid_dim; j++) {
	int k = grid[j+grid_dim*i];
	if (k) {
	  k--;
	  int iVar = varIndex(i, j, k);
	  assert(iVar<num_var);
	  // Nonzero cells have bounds fixed at 1
	  lower[iVar] = 1;
	}
      }
    return_status = highs.changeColsBounds(0, num_var-1, &lower[0], &upper[0]);
    assert(return_status == HighsStatus::kOk);
    // Set up RHS
    rhs.assign(num_row, 1);
    return_status = highs.changeRowsBounds(0, num_row-1, &rhs[0], &rhs[0]);
    assert(return_status == HighsStatus::kOk);

    // Ensure lowest level of development logging
    return_status = highs.setOptionValue("log_dev_level", 1);
    assert(return_status == HighsStatus::kOk);

    // Solve the LP
    problem_num++;
    printf("Solving Sudoku %d\n", problem_num);
    // Ensure that any basis is cleared
    return_status = highs.setBasis();
    assert(return_status == HighsStatus::kOk);
    
    return_status = highs.run();
    assert(return_status == HighsStatus::kOk);
    HighsSolution solution = highs.getSolution();
    
    // Determine whether the solution is fractional, and set up costs
    // and bounds to investigate whether there exists an nontrivial
    // feasible direction
    //
    // Objective is set up so that its value is the negation of the
    // 1-norm of the direction
    int num_fractional = 0;
    const double tl_fractional = 1e-4;
    for (int iVar=0; iVar<num_var;iVar++) {
      if (solution.col_value[iVar] > tl_fractional &&
	  solution.col_value[iVar] < 1-tl_fractional) num_fractional++;
      if (solution.col_value[iVar] < 0.5) {
	// Zero point in feasible region, so component of feasible
	// direction must be non-negative, and cost is negative
	cost[iVar] = -1;
	lower[iVar] = 0;
	upper[iVar] = inf;
      } else {
	// Unit point in feasible region, so component of feasible
	// direction must be non-positive, and cost is positive
	cost[iVar] = 1;
	lower[iVar] = -inf;
	upper[iVar] = 0;
      }
    }
    if (num_fractional) {
      printf("Solution has %d fractional components\n", num_fractional);
      continue;
    }
    // Integer solution, so see whether there exists an nontrivial
    // feasible direction
    return_status = highs.changeColsCost(0, num_var-1, &cost[0]);
    assert(return_status == HighsStatus::kOk);
    return_status = highs.changeColsBounds(0, num_var-1, &lower[0], &upper[0]);
    assert(return_status == HighsStatus::kOk);
    // Set up RHS
    rhs.assign(num_row, 0);
    return_status = highs.changeRowsBounds(0, num_row-1, &rhs[0], &rhs[0]);
    assert(return_status == HighsStatus::kOk);
    // Ensure that any basis is cleared
    return_status = highs.setBasis();
    assert(return_status == HighsStatus::kOk);
    // Solve for a feasible direction
    return_status = highs.run();
    assert(return_status == HighsStatus::kOk);
    HighsSolution direction = highs.getSolution();
    printf("Solution of feasible direction problem has |d| = %g\n", -info.objective_function_value);
    //      if (info.objective_function_value < 0) {
    double alpha = inf;
    for (int iVar=0; iVar < num_var; iVar++) {
      if (direction.col_value[iVar] > 0) {
	alpha = std::min(1/direction.col_value[iVar], alpha);
      } else if (direction.col_value[iVar] < 0) {
	alpha = std::min(-1/direction.col_value[iVar], alpha);
      }
      //	  printf("%3d: x = %11.4g; d = %11.4g; alpha = %11.4g\n",
      //		 iVar, solution.col_value[iVar], direction.col_value[iVar], alpha);
    }
    assert(info.objective_function_value < 0 || alpha == inf);
    if (alpha >= inf) continue;
    // A finite step is possible
    std::vector<double> vertex = solution.col_value;
    num_fractional = 0;
    double max_bound_residual = 0;
    for (int iVar=0; iVar < num_var; iVar++) {
      vertex[iVar] += alpha * direction.col_value[iVar];
      if (vertex[iVar] > tl_fractional &&
	  vertex[iVar] < 1-tl_fractional) num_fractional++;
      // Get the max bound residual
      double bound_residual = std::min(vertex[iVar], 1-vertex[iVar]);
      max_bound_residual = std::max(-bound_residual, max_bound_residual);
      // Get the max row residual
      double max_row_residual = 0;
      std::vector<double> row_residual;
      row_residual.assign(num_row, 0);
      lp.a_matrix_.productQuad(row_residual, vertex);
      for (int iRow=0; iRow < num_row; iRow++)
	max_row_residual = std::max(fabs(row_residual[iRow]-1), max_row_residual);
      
      printf("Feasible direction yields step %g to vertex with %d fractional components and max bound / row residual = %g / %g\n",
	     alpha, num_fractional, max_bound_residual, max_row_residual);      
      
    }
  }
  return 0;
}

int varIndex(const int i, const int j, const int k) { return k + grid_dim*(j+grid_dim*i); }

HighsStatus defineLp(Highs& highs) {
  HighsStatus return_status = HighsStatus::kOk;
  // Add variables to HiGHS model
  for (int k=0; k<num_var; k++) {
    return_status = highs.addVar(0, 1);
    if (return_status != HighsStatus::kOk) return return_status;
  }

  // Add constraints to HiGHS model
  std::vector<int> index(grid_dim);
  std::vector<double> value(grid_dim);
  value.assign(grid_dim, 1);
  for (int i=0; i<grid_dim; i++)
    for (int j=0; j<grid_dim; j++) {
      for (int k=0; k<grid_dim; k++)
	index[k]=varIndex(i, j, k);
      return_status = highs.addRow(1, 1, grid_dim, &index[0], &value[0]);
    if (return_status != HighsStatus::kOk) return return_status;
    }
  for (int i=0; i<grid_dim; i++)
    for (int k=0; k<grid_dim; k++) {
      for (int j=0; j<grid_dim; j++)
	index[j]=varIndex(i, j, k);
      return_status = highs.addRow(1, 1, grid_dim, &index[0], &value[0]);
    if (return_status != HighsStatus::kOk) return return_status;
    }
  for (int k=0; k<grid_dim; k++)
    for (int j=0; j<grid_dim; j++) {
      for (int i=0; i<grid_dim; i++)
	index[i]=varIndex(i, j, k);
      return_status = highs.addRow(1, 1, grid_dim, &index[0], &value[0]);
    if (return_status != HighsStatus::kOk) return return_status;
    }
  for (int k=0; k<grid_dim; k++) {
    for (int istart = 0; istart<3; istart++) {
      for (int jstart = 0; jstart<3; jstart++) {
	int count = 0;
	for (int iloop = 0; iloop<3; iloop++) {
	  for (int jloop = 0; jloop<3; jloop++) {
	    int i = iloop + istart * square_dim;
	    int j = jloop + jstart * square_dim;
	    index[count++] = varIndex(i, j, k);
	  }
	}
	return_status = highs.addRow(1, 1, grid_dim, &index[0], &value[0]);
    if (return_status != HighsStatus::kOk) return return_status;
      }
    }
  }
  // Possibly solve as MIP
  const bool solve_as_mip = false;
  if (solve_as_mip) {
    // Set up integrality
    std::vector<HighsVarType> integrality;
    integrality.assign(num_var, HighsVarType::kInteger);
    return_status = highs.changeColsIntegrality(0, num_var-1, &integrality[0]);
    if (return_status != HighsStatus::kOk) return return_status;
  }
  return return_status;
}
