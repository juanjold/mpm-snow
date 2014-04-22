#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "array3.h"
#include "vec.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "Eigen/Dense"
#include <vector>

class FluidSim {

public:
   void initialize(float width, int ni_, int nj_, int nk_);
   void set_boundary(float (*phi)(const Vec3f&));
   void set_liquid(float (*phi)(const Vec3f&));
   void add_particle(const Vec3f& pos);

   void advance(float dt);

   //Grid dimensions
   int ni,nj,nk;
   float dx;
   
   //Fluid velocity
   Array3f u, v, w;
   Array3f temp_u, temp_v, temp_w;
   Array3f uS, vS, wS;
   
   //Static geometry representation
   Array3f nodal_solid_phi;
   Array3f u_weights, v_weights, w_weights;
   Array3c u_valid, v_valid, w_valid;

   Array3f mass, density;

   std::vector<Vec3f> particles;
   std::vector<Vec3f> particlesVel;
   std::vector<Vec3f> particlesCell;
   std::vector<float> particlesMass;
   std::vector<float> particlesDensity;
   std::vector<float> particlesVolume;
   std::vector<Eigen::Matrix3f > particlesDeformation;
   std::vector<Eigen::Matrix3f > plasticDeformationGradient;
   std::vector<Eigen::Matrix3f > elasticDeformationGradient;
   
   //GridForces constants
    float thetaCompress;      //Crit Compression
    float thetaStretch;      //Crit Stretch
    float hardening;           //Hardening Coff
    float rho0;         //Inital Density
    float youngs0;       //Youngs Modulus 
    float poisson;          //poission Ratio
    float mu0;
    float lambda0;

   /*
   std::vector<float> JE;
   std::vector<float> JP;
   std::vector<Eigen::Matrix3f > RE;
   std::vector<Eigen::Matrix3f > SE;
   */

   float particle_radius;

   Array3f liquid_phi;

   //Data arrays for extrapolation
   Array3c valid, old_valid;

   //Solver data
   PCGSolver<double> solver;
   SparseMatrixd matrix;
   std::vector<double> rhs;
   std::vector<double> pressure;
   
   Vec3f get_velocity(const Vec3f& position);

private:

   Vec3f trace_rk2(const Vec3f& position, float dt);

   float cfl();

   void advect_particles(float dt);
   void advect(float dt);
   void add_force(float dt);
   void project(float dt);
   void constrain_velocity();

   //helpers for pressure projection
   void compute_weights();
   void solve_pressure(float dt);
   void compute_phi();

   // Need a better way to handle intepolation weights, might want to store they are used many many times
   void printParticles(); //Done lol
   void ParticleToGrid(); // (3/4) Done *Include Velocity to grid update (uncomment) 
   void ParticleVolumeDensity(); // Done 
   void ComputeGridForces(float dt); //Not Started *Equation (6) from paper
   void UpdateVelocityOnGrid(float dt); //Not Started
   void GridBasedBodyCollisions(); //Not Started
   void SolveLinearSystem(); //Not Started
   void UpdateDeformation(); //Not Started
   void UpdateParticleVel(); //Not Started
   void ParticleBasedBodyCollisions(); //Not Started
   void UpdateParticlePos(); //Not Started


};


#endif