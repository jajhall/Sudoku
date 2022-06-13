#include "Highs.h"
#include <cassert>
using std::cin;
using std::cout;
using std::endl;

const int grid_dim = 9;
const int num_var = grid_dim * grid_dim * grid_dim;

int varIndex(const int i, const int j, const int k);

int main() {
  std::string sudoku;
  std::vector<int> grid;
  std::vector<int> i_of_var(num_var);
  std::vector<int> j_of_var(num_var);
  std::vector<int> k_of_var(num_var);
  std::vector<int> gridStart;
  gridStart.push_back(0);
  gridStart.push_back(3);
  gridStart.push_back(6);

  for (int i=0; i<grid_dim; i++) {
    for (int j=0; j<grid_dim; j++) {
      cin >> sudoku;
      if (sudoku == ".") {
	grid.push_back(0);
      } else {
	const char* c = sudoku.c_str();
	grid.push_back(atoi(c));
      }
    }
  }
	
  /*
  cin >> sudoku;
  cout << sudoku << endl;

  for (string::iterator it=sudoku.begin(); it!=sudoku.end(); ++it) {

    if (*it == '.'){
      grid.push_back(0);
    } else {
      const char * c = (const char *)*it;
      grid.push_back(atoi(c));
    }
    int k = grid.size()-1;
    printf("%2d: %1d from %1c\n", k, grid[k], sudoku[k]);
  }
  */
  for (int i=0; i<grid_dim; i++)
    for (int j=0; j<grid_dim; j++) 
      for (int k=0; k<grid_dim; k++) {
	int iVar = varIndex(i, j, k);
	i_of_var[iVar] = i+1;
	j_of_var[iVar] = j+1;
	k_of_var[iVar] = k+1;
      }
  
	
  Highs highs;
  std::vector<HighsVarType> integrality;
  for (int k=0; k<num_var; k++) {
    highs.addVar(0, 1);
    integrality.push_back(HighsVarType::kInteger);
  }
  std::vector<int> index(grid_dim);
  std::vector<double> value(grid_dim);
  value.assign(grid_dim, 1);
  for (int i=0; i<grid_dim; i++)
    for (int j=0; j<grid_dim; j++) {
      for (int k=0; k<grid_dim; k++)
	index[k]=varIndex(i, j, k);
      highs.addRow(1, 1, grid_dim, &index[0], &value[0]);
    }
  for (int i=0; i<grid_dim; i++)
    for (int k=0; k<grid_dim; k++) {
      for (int j=0; j<grid_dim; j++)
	index[j]=varIndex(i, j, k);
      highs.addRow(1, 1, grid_dim, &index[0], &value[0]);
    }
  for (int k=0; k<grid_dim; k++)
    for (int j=0; j<grid_dim; j++) {
      for (int i=0; i<grid_dim; i++)
	index[i]=varIndex(i, j, k);
      highs.addRow(1, 1, grid_dim, &index[0], &value[0]);
    }
  for (int k=0; k<grid_dim; k++) {
    for (int istart = 0; istart<3; istart++) {
      for (int jstart = 0; jstart<3; jstart++) {
	int count = 0;
	for (int iloop = 0; iloop<3; iloop++) {
	  for (int jloop = 0; jloop<3; jloop++) {
	    int i = iloop + gridStart[istart];
	    int j = jloop + gridStart[jstart];
	    index[count++] = varIndex(i, j, k);
	  }
	}
	highs.addRow(1, 1, grid_dim, &index[0], &value[0]);
      }
    }
  }
  for (int i=0; i<grid_dim; i++)
    for (int j=0; j<grid_dim; j++) {
      int k = grid[j+grid_dim*i];
      if (k) {
	k--;
	int iVar = varIndex(i, j, k);
	assert(iVar<num_var);
	highs.changeColBounds(iVar, 1, 1);
      }
    }
  HighsLp lp = highs.getLp();
  HighsSparseMatrix& matrix = lp.a_matrix_;
  matrix.ensureRowwise();
  printf("LP has %d rows, %d cols and %d nonzeros\n",  lp.num_row_, lp.num_col_, matrix.numNz());

  //  highs.changeColsIntegrality(0, num_var-1, &integrality[0]);
  highs.setOptionValue("log_dev_level", 1);
  highs.run();
  const HighsSolution& solution = highs.getSolution();
  int num_fractional = 0;
  const double tl_fractional = 1e-4;
  for (int iVar=0; iVar<num_var;iVar++)
    if (solution.col_value[iVar] > tl_fractional && solution.col_value[iVar] < 1-tl_fractional) num_fractional++;
  printf("Solution has %d fractional components\n", num_fractional);
  highs.writeModel("Sudoku.mps");
  return 0;
}

int varIndex(const int i, const int j, const int k) { return k + grid_dim*(j+grid_dim*i); }
