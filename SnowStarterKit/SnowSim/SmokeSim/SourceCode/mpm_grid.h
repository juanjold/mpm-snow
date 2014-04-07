//By Tyler Barkin

#ifndef MPMGrid_H_
#define MPMGrid_H_

#pragma warning(disable: 4244 4267 4996) //Not sure what these are

#include "open_gl_headers.h"
#include "vec.h"
#include "grid_data.h"
#include "grid_vec.h"

class Camera;

class MPMGrid
{

public:
	MPMGrid();
	~MPMGrid();
	MPMGrid(const MPMGrid& orig);
	MPMGrid& operator=(const MPMGrid& orig);

	void reset();

	void draw(const Camera& c);
	void updateSources();
	//void advectVelocity(double dt);
	//void addExternalForces(double dt);
	//void project(double dt);
	//void advectTemperature(double dt);
	//void advectDensity(double dt);

protected:

	// Setup:
	void initialize();

	// Simulation:
	//void computeBouyancy(double dt);
	//void computeVorticityConfinement(double dt);

	// Rendering:
	//struct Cube { vec3 pos; vec4 color; double dist; };
	//void drawWireGrid();
	//void drawSmokeCubes(const Camera& c);
	//void MACGrid::drawObjects(const Camera& c);
	//void drawSmoke(const Camera& c);
	//void drawCube(const MACGrid::Cube& c);
	//void drawFace(const MACGrid::Cube& c);
	//void drawVelocities();
	//vec4 getRenderColor(int i, int j, int k);
	//vec4 getRenderColor(const vec3& pt);
	//void drawZSheets(bool backToFront);
	//void drawXSheets(bool backToFront);

	// GridData accessors:
	enum Direction { X, Y, Z };
	//vec3 getVelocity(const vec3& pt);
	//double getVelocityX(const vec3& pt);
	//double getVelocityY(const vec3& pt);
	//double getVelocityZ(const vec3& pt);
	double getTemperature(const vec3& pt);
	double getDensity(const vec3& pt);
	vec3 getCenter(int i, int j, int k);

	// Sets up the A matrix:
	//void setUpAMatrix();

	// Conjugate gradient stuff:
	bool conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	double dotProduct(const GridData & vector1, const GridData & vector2);
	void add(const GridData & vector1, const GridData & vector2, GridData & result);
	void subtract(const GridData & vector1, const GridData & vector2, GridData & result);
	void multiply(const double scalar, const GridData & vector, GridData & result);
	double maxMagnitude(const GridData & vector);
	void apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result);
	bool isValidCell(int i, int j, int k);

	// Fluid grid cell properties:
	GridDataVec  V; //X,Y,Z Components of Velocity stored at grid nodes, size is dimX*dimY*dimZ
	GridData mP;  // Pressure, stored at grid nodes, size is dimX*dimY*dimZ
	GridData mD;  // Density, stored at grid nodes, size is dimX*dimY*dimZ
	GridData mT;  // Temperature, stored at grid nodes, size is dimX*dimY*dimZ

	

public:

	enum RenderMode { CUBES, SHEETS };
	static RenderMode theRenderMode;
	static bool theDisplayVel;
	
	// Saves smoke in CIS 460 volumetric format:
	void saveSmoke(const char* fileName);

};

#endif