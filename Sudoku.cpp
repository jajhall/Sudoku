#include "Highs.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

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
const bool solve_as_mip = true;
const bool strong_formulation = false;

// kHighsInf is defined to be a logical infinity
//
// const double kHighsInf = std::numeric_limits<double>::infinity();
//
// rather than a very large number - or even the maximum floating
// point number. 
const double inf = kHighsInf;

// "Headers" of methods used in "main", but defined later
int varIndex(const int i, const int j, const int k);
// Defined later in this file, varIndex translates the (i, j, k)
// triple into the corresponding variable number

HighsStatus defineLp(Highs& highs);
// Defined later in this file, defineLp sets up the variables and
// constraints of the LP

void reportPresolveLog(Highs& highs);

void writeSudoku(const std::vector<double>& solution);

// The "main" method - what is called by the executable.
int main(int argc, char *argv[]) {
  // Check that there is a file specified and open it as std::ifstream
  assert(argc == 2);
  const std::string sudoku = argv[1];
  const std::string sudoku_dat = sudoku + ".dat";
  const std::string sudoku_mps = sudoku + ".mps";
  
  std::ifstream f(sudoku_dat);

  // Instantiate HiGHS
  Highs highs;
  // Define a variable to receive the return status from HiGHS
  // method. "HighsStatus" is a special type of idetifier that can
  // only take the following values: HighsStatus::kOk,
  // HighsStatus::kWarning, HighsStatus::kError
  HighsStatus return_status;

  // Define the level of output from HiGHS
  const bool allow_highs_dev_log = true;
    if (allow_highs_dev_log) {
    // Ensure lowest level of development logging
    return_status = highs.setOptionValue("log_dev_level", 1);
    // assert(bool) kills the executable when compiled with debugging
    // (-g flag in CompLinkFILE) if bool is false. In other words, the
    // programmer is "asserting" that something must be true.
    assert(return_status == HighsStatus::kOk);
  } else {
    // Force HiGHS to run quietly
    return_status = highs.setOptionValue("output_flag", false);
    assert(return_status == HighsStatus::kOk);
    highs.setOptionValue("presolve_log_report", true);
    highs.setOptionValue("doubleton_equation_sudoku_report", true);
  }

  // Define the LP variables and constraints without the starting grid
  return_status = defineLp(highs);
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
  // Set up const reference to internal LP, presolved LP, info structure, and solution
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
  const HighsLp& presolved_lp = highs.getPresolvedLp();
  const HighsSolution& solution = highs.getSolution();
  // Convenient short-hand for the number of rows
  const int num_row = lp.num_row_;
  // Note that the following line
  //
  // HighsInfo local_info = highs.getInfo();
  //
  // takes a copy of the "HighsInfo" values in the HiGHS class. If
  // calling a HiGHS method causes the internal values to change, then
  // these copies are unaffected.

  // The presolve_rule_off option will alow certain presolve rules to
  // be switched off
  int presolve_rule_off = 0;
  //  presolve_rule_off += 512;  // Doubleton equation rule off
  //  presolve_rule_off += 1024; // Dependent equations rule off
  //  presolve_rule_off += 2048; // Dependent free columns rule off (not implemented)
  //  presolve_rule_off += 4096; // Aggregator rule off
  //  presolve_rule_off += 8192; // Parallel rows and columns rule off
  if (presolve_rule_off) {
    return_status = highs.setOptionValue("presolve_rule_off", presolve_rule_off);
    assert(return_status == HighsStatus::kOk);
  }

  // Set up vectors of bounds to be re-used, with upper being constant
  std::vector<double> lower;
  std::vector<double> upper;
  upper.assign(num_var, 1);

  int problem_num = 0;
  // Read a sequence of Sudoku grids
  std::vector<int> grid;
  grid.assign(num_cell, 0);
  for (;;) {
    // This is an "infinite" loop, that terminates when a "break" is
    // encountered
    std::string sudoku;
    std::getline(f, sudoku);
    //      printf("Read Sudoku line %s\n", sudoku.c_str());

    int size = sudoku.length();
    if (size < num_cell) {
      if (size) printf("STRANGE: Line has size %d\n", size);
      break; // Leave the "for (;;)" loop
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
    if (strong_formulation) {
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
      return_status = highs.changeColsBounds(0, num_var-1, lower.data(), upper.data());
      assert(return_status == HighsStatus::kOk);
    } else {
      // Set up an objective that sums negation of all variables, and
      // single constraint that the sum of all variables corresponding
      // to initial values is the number of initial values
      std::vector<double> cost(num_var, -1);
      return_status = highs.changeColsCost(0, num_var-1, cost.data());
      assert(return_status == HighsStatus::kOk);
      std::vector<int> index;
      std::vector<double> value;
      for (int i=0; i<grid_dim; i++) {
	for (int j=0; j<grid_dim; j++) {
	  int k = grid[j+grid_dim*i];
	  if (k) { // if (k) is true if k is nonzero
	    k--;
	    int iVar = varIndex(i, j, k);
	    assert(iVar<num_var);
	    index.push_back(iVar);
	    value.push_back(1);
	  }
	}
      }
      const int num_start_value = index.size();
      const double start_value_rhs = num_start_value;
      return_status = highs.addRow(start_value_rhs, start_value_rhs, num_start_value, index.data(), value.data());
      assert(return_status == HighsStatus::kOk);
    }

    // Increment the problem counter
    problem_num++;
    printf("\nSolving Sudoku %d\n", problem_num);
    const bool just_presolve = false;
    if (just_presolve) {
      // Just run presolve 
      return_status = highs.presolve();
      assert(return_status == HighsStatus::kOk);
      // Flag up case where presolve does not reduce to empty
      if (presolved_lp.num_row_+presolved_lp.num_col_)
	printf("Presolved LP has %d rows and %d columns\n",
	       (int)presolved_lp.num_row_, (int)presolved_lp.num_col_);
    } else {
      // Solve the LP
      // Ensure that any basis is cleared - otherwise HiGHS won't use
      // presolve when given the second and subsequent LPs
      return_status = highs.setBasis();
      assert(return_status == HighsStatus::kOk);
    
      return_status = highs.writeModel(sudoku_mps);
      // Highs::run() actually does the optimization!
      return_status = highs.run();
      assert(return_status == HighsStatus::kOk);
      HighsModelStatus model_status = highs.getModelStatus();
      if (model_status != HighsModelStatus::kOptimal)
	printf("Model status is %s\n", highs.modelStatusToString(model_status).c_str());
      // Report whether simplex iterations were required
      if (info.simplex_iteration_count)
	printf("Solving problem required %d simplex iterations\n", info.simplex_iteration_count);
   
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
      } else {
	writeSudoku(solution.col_value);
      }
    }
    if (!solve_as_mip) {
      // Report the presolve log: only valid for LPs
      reportPresolveLog(highs);
    }
  }
  
  return 0;
}

int varIndex(const int i, const int j, const int k) { return k + grid_dim*(j+grid_dim*i); }

HighsStatus defineLp(Highs& highs) {
  HighsStatus return_status = HighsStatus::kOk;
  const double constraint_lower = strong_formulation ? 1 : -kHighsInf;
  
  // Add variables to HiGHS model
  std::vector<double> lower;
  std::vector<double> upper;
  // Set up bounds
  lower.assign(num_var, 0);
  upper.assign(num_var, 1);
  return_status = highs.addVars(num_var, lower.data(), upper.data());
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
      return_status = highs.addRow(constraint_lower, 1, grid_dim, index.data(), value.data());
      if (return_status != HighsStatus::kOk) return return_status;
    }
  for (int i=0; i<grid_dim; i++)
    for (int k=0; k<grid_dim; k++) {
      for (int j=0; j<grid_dim; j++)
	index[j]=varIndex(i, j, k);
      return_status = highs.addRow(constraint_lower, 1, grid_dim, index.data(), value.data());
      if (return_status != HighsStatus::kOk) return return_status;
    }
  for (int k=0; k<grid_dim; k++)
    for (int j=0; j<grid_dim; j++) {
      for (int i=0; i<grid_dim; i++)
	index[i]=varIndex(i, j, k);
      return_status = highs.addRow(constraint_lower, 1, grid_dim, index.data(), value.data());
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
	return_status = highs.addRow(constraint_lower, 1, grid_dim, index.data(), value.data());
	if (return_status != HighsStatus::kOk) return return_status;
      }
    }
  }
  // Possibly solve as MIP
  if (solve_as_mip) {
    // Set up integrality
    std::vector<HighsVarType> integrality;
    integrality.assign(num_var, HighsVarType::kInteger);
    return_status = highs.changeColsIntegrality(0, num_var-1, integrality.data());
    if (return_status != HighsStatus::kOk) return return_status;
  }
  return return_status;
}

void reportPresolveLog(Highs& highs) {
  const HighsPresolveLog& presolve_log = highs.getPresolveLog();
  // The presolve_rule_off option will alow certain presolve rules to
  // be switched off
  int presolve_rule_off = highs.getOptions().presolve_rule_off;
  int bit = 1;
  printf("\nRule  Bit| Call Row Col| Name\n");
  for (int rule_ix = 0; rule_ix < kPresolveRuleCount; rule_ix++) {
    bool allow = !(presolve_rule_off & bit);
    const HighsPresolveRuleLog& log = presolve_log.rule[rule_ix];
    if (log.call) printf("  %2d %4d|  %3d %3d %3d| %s\n", rule_ix, bit,
			 log.call, 
			 log.row_removed,
			 log.col_removed,
			 highs.presolveRuleTypeToString(rule_ix).c_str());
    bit *= 2;
  }
}
void writeSudoku(const std::vector<double>& solution) {
  std::vector<int>values(grid_dim);
  for (int i=0; i<grid_dim; i++) {
    for (int j=0; j<grid_dim; j++) {
      values[j] = 0;
      for (int k=0; k<grid_dim; k++) {
	if (solution[varIndex(i, j, k)] > 1e-3) {
	  values[j] = k+1;
	  break;
	}
      }
    }
    for (int j=0; j<grid_dim; j++)
      printf("%1d", values[j]);
    printf("\n");
  }
}
 
	

