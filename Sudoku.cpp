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

// Set up variables to refer to the Sudoku grid and decision variables
// in the LP. These are defined "const" to stress that they don't
// change.
const int square_dim = 3;
const int grid_dim = square_dim * square_dim;
const int num_cell = grid_dim * grid_dim;
const int num_var = num_cell * grid_dim;

// kHighsInf is defined to be a logical infinity
//
// const double kHighsInf = std::numeric_limits<double>::infinity();
//
// rather than a very large number - or even the maximum floating
// point number. 
const double inf = kHighsInf;

int varIndex(const int i, const int j, const int k);
// Defined later in this file, varIndex translates the (i, j, k)
// triple into the corresponding variable number

HighsStatus defineLp(Highs& highs);
// Defined later in this file, defineLp sets up the variables and
// constraints of the LP

int main(int argc, char *argv[]) {
  // Check that there is a file specified and open it as std::ifstream
  assert(argc == 2);
  std::ifstream f(argv[1]);

  // Instantiate HiGHS
  Highs highs;
  // Define the LP variables and constraints
  HighsStatus return_status = defineLp(highs);
  assert(return_status == HighsStatus::kOk);

  // A pointer is a variable that holds the memory address of another
  // variable. It can change if, for example, the size of a vector is
  // increased and this is only possibly by using a new section of
  // memory. So, using copies of pointers can be dangerous
  //
  // A reference is an alias for an already existing variable. It
  // doesn't just point to the variable by storing its address, it's
  // equivalent to that variable so, for example, if that variable
  // moves in memory, the reference follows it.
  //
  // Set up const reference to internal info structure and LP
  //
  // There are variables in the HiGHS class [strictly they are called
  // class data members] that store data for the class. Users can get
  // access to the stored values using calls like the following
  // three. However, we can't allow users to modify the values, so
  // "const" references are returned. If calling a HiGHS method causes
  // the internal values to change, then the values available from
  // these "const" references also change - because they are
  // references to the original variables, not copies of them.
  const HighsInfo& info = highs.getInfo();
  const HighsLp& lp = highs.getLp();
  const int num_row = lp.num_row_;
  // Note that the following line
  //
  // HighsInfo local_info = highs.getInfo();
  //
  // takes a copy of the "HighsInfo" values in the HiGHS class. If
  // calling a HiGHS method causes the internal values to change, then
  // these copies are unaffected.

  // Set up vectors of bounds to be re-used, with upper being constant
  std::vector<double> lower;
  std::vector<double> upper;
  upper.assign(num_var, 1);

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
      } else {
	grid[k] = 0;
      }
      //    printf("%2d: %1d from %1c\n", k, grid[k], sudoku[k]);
    }
    // Set up bounds
    lower.assign(num_var, 0);
    for (int i=0; i<grid_dim; i++)
      for (int j=0; j<grid_dim; j++) {
	int k = grid[j+grid_dim*i];
	if (k) { // if (k) is true if k is nonzero
	  k--;
	  int iVar = varIndex(i, j, k);
	  assert(iVar<num_var);
	  // Nonzero cells have bounds fixed at 1
	  lower[iVar] = 1;
	}
      }
    // Change the bounds of the incumbent LP in the Highs class
    // instance. Highs::changeColsBounds takes pointers to the (first)
    // value in the lower and upper bound vectors so that the C
    // interface to HiGHS doesn't have to make std::vector<double>
    // copies of the bounds when called from C (which only has arrays
    // of memory).
    return_status = highs.changeColsBounds(0, num_var-1, &lower[0], &upper[0]);
    assert(return_status == HighsStatus::kOk);

    // Ensure lowest level of development logging
    return_status = highs.setOptionValue("log_dev_level", 1);
    assert(return_status == HighsStatus::kOk);

    // Solve the LP
    problem_num++;
    printf("Solving Sudoku %d\n", problem_num);
    // Ensure that any basis is cleared - otherwise HiGHS won't use
    // presolve when given the second and subsequent LPs
    return_status = highs.setBasis();
    assert(return_status == HighsStatus::kOk);
    
    // Highs::run() actually does the optimization!
    return_status = highs.run();
    assert(return_status == HighsStatus::kOk);
    HighsSolution solution = highs.getSolution();
    
    // Determine whether the solution is fractional
    int num_fractional = 0; // This has to change so isn't "const"
    const double tl_fractional = 1e-4; 
    for (int iVar=0; iVar<num_var;iVar++) {
      if (solution.col_value[iVar] > tl_fractional &&
	  solution.col_value[iVar] < 1-tl_fractional) num_fractional++;
    }
    if (num_fractional) {
      printf("Solution has %d fractional components\n", num_fractional);
      continue;
    }
  }
  return 0;
}

int varIndex(const int i, const int j, const int k) { return k + grid_dim*(j+grid_dim*i); }

HighsStatus defineLp(Highs& highs) {
  HighsStatus return_status = HighsStatus::kOk;
  // Add variables to HiGHS model
  std::vector<double> lower;
  std::vector<double> upper;
  // Set up bounds
  lower.assign(num_var, 0);
  upper.assign(num_var, 1);
  return_status = highs.addVars(num_var, &lower[0], &upper[0]);
  if (return_status != HighsStatus::kOk) return return_status;

  // Add constraints to HiGHS model. We know how many nonzeros are in
  // each row of the constraint matrix, so we can allocate the size of
  // index and value now
  std::vector<int> index(grid_dim);
  std::vector<double> value(grid_dim);
  // Indeed, we know that all the values are 1
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
