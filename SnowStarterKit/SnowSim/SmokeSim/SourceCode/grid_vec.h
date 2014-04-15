// Modified by Tyler Barkin

#ifndef GRID_DATA_VEC
#define GRID_DATA_VEC

#include "grid_data.h"

class GridDataVec {
private:
protected:
public:
	GridDataVec() {
		mU.initialize();
		mV.initialize();
		mW.initialize();
	}
	GridData mU;
	GridData mV;
	GridData mW;
};

#endif // GRID_VEC