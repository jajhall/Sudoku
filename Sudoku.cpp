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

  // Set up const reference to internal solution structure
  const HighsSolution& solution = highs.getSolution();

  // Set up vectors of bounds to be re-used, with upper being constant
  std::vector<double> lower(num_var);
  std::vector<double> upper(num_var);
  upper.assign(num_var, 1);

  // Read a sequence of Sudoku grids
  std::vector<int> grid;
  grid.assign(num_cell, 0);
  for (;;) {
    std::stringstream ss;
    ss << f.rdbuf();
    
    std::string sudoku = std::move(ss.str());
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
    // Fix lower bounds of nonzero cells
    lower.assign(num_var, 0);
    for (int i=0; i<grid_dim; i++)
      for (int j=0; j<grid_dim; j++) {
	int k = grid[j+grid_dim*i];
	if (k) {
	  k--;
	  int iVar = varIndex(i, j, k);
	  assert(iVar<num_var);
	  lower[iVar] = 1;
	}
      }
    // Set the lower and upper bounds for this LP
    highs.changeColsBounds(0, num_var-1, &lower[0], &upper[0]);

    /*
    HighsLp lp = highs.getLp();
    HighsSparseMatrix& matrix = lp.a_matrix_;
    //    matrix.ensureRowwise();
    printf("LP has %d rows, %d cols and %d nonzeros\n",  lp.num_row_, lp.num_col_, matrix.numNz());
    */

    // Ensure lowest level of development logging
    highs.setOptionValue("log_dev_level", 1);

    // Solve the LP
    highs.run();
    
    // Determine whether the solution is fractional
    int num_fractional = 0;
    const double tl_fractional = 1e-4;
    for (int iVar=0; iVar<num_var;iVar++)
      if (solution.col_value[iVar] > tl_fractional &&
	  solution.col_value[iVar] < 1-tl_fractional) num_fractional++;
    if (num_fractional)
      printf("Solution has %d fractional components\n", num_fractional);
    //  highs.writeModel("Sudoku.mps");
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
