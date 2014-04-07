// Written by Peter Kutz.
// See Bridson's Fluids Notes, section 4.3.1 Putting Them In Matrix-Vector Form (on page 31, or PDF page 43)
// for an explanation of this symmetric sparse matrix data structure.

#ifndef GRID_DATA_VEC
#define GRID_DATA_VEC

#include "grid_data.h"

class GridDataVec {
private:
protected:
public:
	GridDataVec() {
		U.initialize();
		V.initialize();
		W.initialize();
	}
	GridData U;
	GridData V;
	GridData W;
};

#endif // GRID_VEC