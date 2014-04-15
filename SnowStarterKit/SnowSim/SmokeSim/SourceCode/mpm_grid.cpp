// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>


// Globals:
MACGrid target;
double density = 1;
int numberOfObjects=0;
vec3 object1S[6] = {vec3(4,3,1),vec3(2,3,1),vec3(8,0,0),vec3(11,0,3),vec3(0,9,0),vec3(4,9,0)};
vec3 object1E[6] = {vec3(5,11,7),vec3(4,5,7),vec3(8,4,5),vec3(11,4,8),vec3(2,9,8),vec3(12,9,8)};

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 






MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);

   setUpAMatrix();
}
void MACGrid::drawObjects(const Camera& c){
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
      MACGrid::Cube cube;
	  for(int l = 0; l<numberOfObjects; l++){
		 for(int i = object1S[l][0]; i <= object1E[l][0]; i++){
			for(int j = object1S[l][1]; j <= object1E[l][1]; j++){
			  for(int k = object1S[l][2]; k <= object1E[l][2]; k++){
	  cube.color = vec4(0.9,1.0,1.0,1.0);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
			  }
		  }
	  }
	}

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
	
    // TODO: Set initial values for density, temperature, and velocity.
	int SX = 6;
	int SY = 1;
	int SZ = 4;
	
	int SX2 = 6;
	int SY2 = 0;
	int SZ2 = 5;
	
	mT(SX-1,SY-1,SZ-1) = 30; 
	mD(SX-1,SY-1,SZ-1) = 1;
	mP(SX-1,SY-1,SZ-1) = 0;
	mU(SX,SY,SZ) = 0;
	mV(SX,SY,SZ) = 0.001;
	mW(SX,SY,SZ) = 0;

	//mT(SX2,SY2,SZ2) = 50; 
	//mD(SX2,SY2,SZ2) = 10;
	//mP(SX2,SY2,SZ2) = 100;
	mU(SX2,SY2,SZ2) = 0;
	mV(SX2,SY2,SZ2) = 0.001;
	mW(SX2,SY2,SZ2) = 0;

	/*for(int i=0;i<theDim[MACGrid::X];i++){
		for(int k=0;k<theDim[MACGrid::Z];k++){
			if(k==0 || k==1 || k == 2)
				mU(i,8,k)=-.25;
		}
	}*/

	FOR_EACH_CELL{
		if(j>5 && j<8){
			//mT(i,j,k) = mT(i,j,k)*(.999);
			//mD(i,j,k) = mD(i,j,k)*(.99);
		}
		if(j>=8 && j<=10){
			//mT(i,j,k) = mT(i,j,k)*(.995);
			//mD(i,j,k) = mD(i,j,k)*(1.001);
		}
		if(j>10){ //upper box cooling
			mT(i,j,k) = mT(i,j,k)-.05;
		}

	}

}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target.
	vec3 Xg,Yg,Zg,Xp,Yp,Zp;
	vec3 qx,qy,qz;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	

	FOR_EACH_FACE{ //Triple Loop with i j and k
		Xg = getCenter(i,j,k)-(vec3(1 , 0 , 0))*(theCellSize); //X face Current
		qx = getVelocity(Xg); //Current Velocity at Face
		Xp = Xg - dt*qx;  //Reverse projected point
		if(j<theDim[MACGrid::Y]+1 && k<theDim[MACGrid::Z]+1)
			target.mU(i,j,k)=getVelocityX(Xp)/2; //Current grid face set to Xp velocity U
		Xg = getCenter(i,j,k)-(vec3(0 , 0 , 0))*(theCellSize); //X face Current
		qx = getVelocity(Xg); //Current Velocity at Face
		Xp = Xg - dt*qx;  //Reverse projected point
		if(j<theDim[MACGrid::Y]+1 && k<theDim[MACGrid::Z]+1)
			target.mU(i,j,k)+=getVelocityX(Xp)/2; //Current grid face set to Xp velocity U

		Yg = getCenter(i,j,k)-(vec3(0 , 1 , 0))*(theCellSize); //Y face Current
		qy = getVelocity(Yg); //Current Velocity at Face
		Yp = Yg - dt*qy; //Reverse projected point
		if(i<theDim[MACGrid::X]+1 && k<theDim[MACGrid::Z]+1)
			target.mV(i,j,k)=getVelocityY(Yp)/2;; //Current grid face set to Yp velocity V
		Yg = getCenter(i,j,k)-(vec3(0 , 0 , 0))*(theCellSize); //Y face Current
		qy = getVelocity(Yg); //Current Velocity at Face
		Yp = Yg - dt*qy; //Reverse projected point
		if(i<theDim[MACGrid::X]+1 && k<theDim[MACGrid::Z]+1)
			target.mV(i,j,k)+=getVelocityY(Yp)/2;; //Current grid face set to Yp velocity V

		Zg = getCenter(i,j,k)-(vec3(0 , 0 , 1))*(theCellSize); //Z face Current
		qz = getVelocity(Zg); //Current Velocity at Face
		Zp = Zg - dt*qz; //Reverse projected point
		if(i<theDim[MACGrid::X]+1 && j<theDim[MACGrid::Y]+1)
			target.mW(i,j,k)=getVelocityZ(Zp)/2;; //Current grid face set to Zp velocity W
		Zg = getCenter(i,j,k)-(vec3(0 , 0 , 0))*(theCellSize); //Z face Current
		qz = getVelocity(Zg); //Current Velocity at Face
		Zp = Zg - dt*qz; //Reverse projected point
		if(i<theDim[MACGrid::X]+1 && j<theDim[MACGrid::Y]+1)
			target.mW(i,j,k)+=getVelocityZ(Zp)/2;; //Current grid face set to Zp velocity W

		/*if(i==3 && j==2 && k==2){
			vec3 zerozerozero=getCenter(1,1,1);
			printf("%d,%d,%d\n",i,j,k);
			printf("%f,%f,%f\n",zerozerozero[0],zerozerozero[1],zerozerozero[2]);
			//printf("%f,%f,%f\n",Xg,Yg,Zg);
			printf("Qx : %f,%f,%f\n",qx[0],qx[1],qx[2]);
			//Can be replaced with RK2 or RK4 from the first project.... Maybe Later
			//printf("%f,%f,%f\n\n",Xp,Yp,Zp);
		}*/
	}

    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
	
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.
	vec3 Xg,Xp,q1,q2;
	target.mT = mT;
	FOR_EACH_CELL{
		Xg = getCenter(i,j,k); //Current center of grid
		q1 = getVelocity(Xg); //Current Velocity
		Xp = Xg - dt*q1; //Reverse projected point, Can be replaced with RK2 or RK4 from the first project.... Maybe Later
		target.mT(i,j,k) = getTemperature(Xg)*0.2 + getTemperature(Xp)*0.8;
	}
    // Then save the result to our object.
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target.

	vec3 Xg,Xp,q1,q2;
	target.mD = mD;
	FOR_EACH_CELL{
		Xg = getCenter(i,j,k); //Current center of grid
		q1 = getVelocity(Xg); //Current Velocity
		Xp = Xg - dt*q1; //Reverse projected point, Can be replaced with RK2 or RK4 from the first project.... Maybe Later
		target.mD(i,j,k) = getDensity(Xg)*0.2 + getDensity(Xp)*0.8;
	}
    // Then save the result to our object.
    mD = target.mD;
	
}

void MACGrid::computeBouyancy(double dt)
{	target.mV = mV;
	double alpha = 7;
	double beta = 7;
	double Tamb = 0;
	GridData f_buoy;
	f_buoy.initialize();
	FOR_EACH_CELL{
		f_buoy(i,j,k) = -alpha*mD(i,j,k) + beta*(mT(i,j,k)-Tamb);

		if(j<=theDim[MACGrid::Y] && j>0)
			target.mV(i,j,k)=dt*(f_buoy(i,j,k)+f_buoy(i,j-1,k))/2;
		else 
			target.mV(i,j,k)=mV(i,j,k);
	}
	
   // Then save the result to our object.
   mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	// Apply the forces to the current velocity and store the result in target.
	double epsilon=1;
	
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	GridData omega[3];
	omega[0].initialize();
	omega[1].initialize();
	omega[2].initialize();
	GridData Gradomega[3];
	Gradomega[0].initialize();
	Gradomega[1].initialize();
	Gradomega[2].initialize();
	GridData NormGradomega[3];
	NormGradomega[0].initialize();
	NormGradomega[1].initialize();
	NormGradomega[2].initialize();
	GridData Vforce[3];
	Vforce[0].initialize();
	Vforce[1].initialize();
	Vforce[2].initialize();
	FOR_EACH_CELL{
		if(i<theDim[MACGrid::X]-1 && i>1 && j<theDim[MACGrid::Y]-1 && j>1 && k<theDim[MACGrid::Z]-1 && k>1){
			omega[0](i,j,k) = (getVelocityZ(getCenter(i,j+1,k))-getVelocityZ(getCenter(i,j-1,k)))/(2*theCellSize) - (getVelocityY(getCenter(i,j,k+1)) - getVelocityY(getCenter(i,j,k-1)))/(2*theCellSize);
			omega[1](i,j,k) = (getVelocityX(getCenter(i,j,k+1))-getVelocityX(getCenter(i,j,k-1)))/(2*theCellSize) - (getVelocityZ(getCenter(i+1,j,k)) - getVelocityZ(getCenter(i-1,j,k)))/(2*theCellSize);
			omega[2](i,j,k) = (getVelocityY(getCenter(i+1,j,k))-getVelocityY(getCenter(i-1,j,k)))/(2*theCellSize) - (getVelocityX(getCenter(i,j+1,k)) - getVelocityX(getCenter(i,j-1,k)))/(2*theCellSize);
		//printf("Omega1 : %f,%f,%f\n",omega[0](i,j,k),omega[1](i,j,k),omega[2](i,j,k));
		}else{
			omega[0](i,j,k) = 0;
			omega[1](i,j,k) = 0;
			omega[2](i,j,k) = 0;
		}
	}
	FOR_EACH_CELL{
		if(i<theDim[MACGrid::X]-1 && i>1 && j<theDim[MACGrid::Y]-1 && j>1 && k<theDim[MACGrid::Z]-1 && k>1){
			double test = std::pow(2.,2.);
			Gradomega[0](i,j,k) = (sqrt(pow(omega[0](i+1,j,k),2)+pow(omega[1](i+1,j,k),2)+pow(omega[2](i+1,j,k),2)) - sqrt(pow(omega[0](i-1,j,k),2)+pow(omega[1](i-1,j,k),2)+pow(omega[2](i-1,j,k),2)))/(2*theCellSize);
			Gradomega[1](i,j,k) = (sqrt(pow(omega[0](i,j+1,k),2)+pow(omega[1](i,j+1,k),2)+pow(omega[2](i,j+1,k),2)) - sqrt(pow(omega[0](i,j-1,k),2)+pow(omega[1](i,j-1,k),2)+pow(omega[2](i,j-1,k),2)))/(2*theCellSize);
			Gradomega[2](i,j,k) = (sqrt(pow(omega[0](i,j,k+1),2)+pow(omega[1](i,j,k+1),2)+pow(omega[2](i,j,k+1),2)) - sqrt(pow(omega[0](i,j,k-1),2)+pow(omega[1](i,j,k-1),2)+pow(omega[2](i,j,k-1),2)))/(2*theCellSize);
		//printf("N1 : %f,%f,%f\n",NormGradomega[0](i,j,k),NormGradomega[1](i,j,k),NormGradomega[2](i,j,k));
		}else{
			Gradomega[0](i,j,k)=0;
			Gradomega[1](i,j,k)=0;
			Gradomega[2](i,j,k)=0;
		}
		NormGradomega[0](i,j,k) = Gradomega[0](i,j,k)/sqrt(pow(Gradomega[0](i,j,k),2)+pow(Gradomega[1](i,j,k),2)+pow(Gradomega[2](i,j,k),2));
		NormGradomega[1](i,j,k) = Gradomega[1](i,j,k)/sqrt(pow(Gradomega[0](i,j,k),2)+pow(Gradomega[1](i,j,k),2)+pow(Gradomega[2](i,j,k),2));
		NormGradomega[2](i,j,k) = Gradomega[2](i,j,k)/sqrt(pow(Gradomega[0](i,j,k),2)+pow(Gradomega[1](i,j,k),2)+pow(Gradomega[2](i,j,k),2));

		Vforce[0](i,j,k)=epsilon*theCellSize*(NormGradomega[1](i,j,k)*omega[2](i,j,k)-NormGradomega[2](i,j,k)*omega[1](i,j,k));
		Vforce[1](i,j,k)=epsilon*theCellSize*(NormGradomega[2](i,j,k)*omega[0](i,j,k)-NormGradomega[0](i,j,k)*omega[2](i,j,k));
		Vforce[2](i,j,k)=epsilon*theCellSize*(NormGradomega[0](i,j,k)*omega[1](i,j,k)-NormGradomega[1](i,j,k)*omega[0](i,j,k));
		
		if(!(Vforce[0](i,j,k) != Vforce[0](i,j,k))){
			//printf("F : %f,%f,%f\n",Vforce[0](i,j,k),Vforce[1](i,j,k),Vforce[2](i,j,k));
			target.mU(i,j,k)=(mU(i,j,k) + dt*Vforce[0](i,j,k))*.95;
			target.mV(i,j,k)=(mV(i,j,k) + dt*Vforce[1](i,j,k))*.95;
			target.mW(i,j,k)=(mW(i,j,k) + dt*Vforce[2](i,j,k))*.95;
		}
	}



	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
	
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
   target.mV = mV; 
	FOR_EACH_FACE{ //gravity
		target.mV(i,j,k) =mV(i,j,k)-mD(i,j,k)*pow(theCellSize,3)*98.1*dt;
	}
	mV=target.mV;
   //int x=0,y=0,state=0;
   //glutMouseFunc(*editVels(x,y,mU,mV,mW),GLUT_LEFT_BUTTON,state,x,y)

}
//void editVels(int dx,int dy, GridData mU,GridData mV,GridData mW){
//}



void MACGrid::project(double dt)
{
	// TODO: Solve Ap = d for pressure.
	// 1. Construct d

	double xp,xm,ym,yp,zm,zp;
	double u_solid=0;
	GridData d;
	d.initialize();

	FOR_EACH_CELL{

		if (i-1 >= 0)
			xm = mU(i-1,j,k);
		else 
			xm = u_solid;
		if (i+1 < theDim[MACGrid::X])
			xp = mU(i,j,k);
		else 
			xp = u_solid;
		if (j-1 >= 0)
			ym = mV(i,j-1,k);
		else 
			ym = u_solid;
		if (j+1 < theDim[MACGrid::Y])
			yp = mV(i,j,k);
		else 
			yp = u_solid;
		if (k-1 >= 0)
			zm = mW(i,j,k-1);
		else 
			zm = u_solid;
		if (k+1 < theDim[MACGrid::Z])
			zp = mW(i,j,k);
		else 
			zp = u_solid;
		for(int l = 0; l<numberOfObjects; l++){
		if(i+1>=object1S[l][0] && i+1<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
			xp =u_solid;	
		}
		if(i-1>=object1S[l][0] && i-1<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
			xm = u_solid;	
		}
		if(i>=object1S[l][0] && i<=object1E[l][0] && j+1>=object1S[l][1] && j+1<=object1E[l][1]  && k>=object1S[l][2] && k<=object1E[l][2] ){
			yp = u_solid;		
		}
		if(i>=object1S[l][0] && i<=object1E[l][0] && j-1>=object1S[l][1] && j-1<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
			ym =u_solid;		
		}
		if(i>=object1S[l][0] && i<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k+1>=object1S[l][2] && k+1<=object1E[l][2] ){
			zp = u_solid;	
		}
		if(i>=object1S[l][0] && i<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k-1>=object1S[l][2] && k-1<=object1E[l][2] ){
			zm = u_solid;	
		}
		if(i>=object1S[l][0] && i<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
			xm=u_solid;
			xp=u_solid;
			ym=u_solid;
			yp=u_solid;
			zm=u_solid;
			zp=u_solid;
		}
		}
			d(i,j,k)=-(((xp-xm)+(yp-ym)+(zp-zm)))*(theCellSize)/dt;

		}

		
		
	// 2. Construct A (This is where objects should be edited in
		
	// 3. Solve for p


	conjugateGradient(AMatrix,mP,d,500,0.001);


	// Subtract pressure from our velocity and save in target.
	target.mP = mP;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	vec3 Xg,Yg,Zg,Xp,Yp,Zp,qx,qy,qz;

	FOR_EACH_CELL{ //Triple Loop with i j and k
		
		if(i<=theDim[MACGrid::X] && i>0)
			target.mU(i,j,k)=mU(i,j,k) - (dt*(mP(i,j,k)-mP(i-1,j,k)))*theCellSize; //Current grid face set to Xp velocity U
		else 
			target.mU(i,j,k)=mU(i,j,k) - (dt*(mP(i,j,k)-mP(i,j,k)))*theCellSize;
		if(j<=theDim[MACGrid::Y] && j>0)
			target.mV(i,j,k)=mV(i,j,k) - (dt*(mP(i,j,k)-mP(i,j-1,k)))*theCellSize; //Current grid face set to Yp velocity V
		else 
			target.mV(i,j,k)=mV(i,j,k) - (dt*(mP(i,j,k)-mP(i,j,k)))*theCellSize;
		if(k<=theDim[MACGrid::Z] && k>0)
			target.mW(i,j,k)=mW(i,j,k) - (dt*(mP(i,j,k)-mP(i,j,k-1)))*theCellSize; //Current grid face set to Zp velocity W
		else
			target.mW(i,j,k)=mW(i,j,k) - (dt*(mP(i,j,k)-mP(i,j,k)))*theCellSize;
		for(int l = 0; l<numberOfObjects; l++){
		if(i>=object1S[l][0] && i<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
			target.mU(i,j,k)=0;//mU(i,j,k) - (dt*(mP(i,j,k)-mP(i,j,k)))*theCellSize;
			target.mV(i,j,k)=0;//mV(i,j,k) - (dt*(mP(i,j,k)-mP(i,j,k)))*theCellSize;
			target.mW(i,j,k)=0;//mW(i,j,k) - (dt*(mP(i,j,k)-mP(i,j,k)))*theCellSize;
		}
		}

	}


	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
	//IMPLEMENT THIS AS A SANITY CHECK: assert(checkDivergence());
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::setUpAMatrix() {
	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}

			
		// Set the diagonal:
			for(int l = 0; l<numberOfObjects; l++){
			if(i+1>=object1S[l][0] && i+1<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
				numFluidNeighbors--;
				AMatrix.plusI(i,j,k) = 0;
			}
			if(i-1>=object1S[l][0] && i-1<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
				numFluidNeighbors--;
				AMatrix.plusI(i-1,j,k) = 0;
			}
		
			if(i>=object1S[l][0] && i<=object1E[l][0] && j+1>=object1S[l][1] && j+1<=object1E[l][1]  && k>=object1S[l][2] && k<=object1E[l][2] ){
				numFluidNeighbors--;
				AMatrix.plusJ(i,j,k) = 0;
			}
			if(i>=object1S[l][0] && i<=object1E[l][0] && j-1>=object1S[l][1] && j-1<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
				numFluidNeighbors--;
				AMatrix.plusJ(i,j-1,k) = 0;
			}
			if(i>=object1S[l][0] && i<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k+1>=object1S[l][2] && k+1<=object1E[l][2] ){
				numFluidNeighbors--;
				AMatrix.plusK(i,j,k) = 0;
			}
			if(i>=object1S[l][0] && i<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k-1>=object1S[l][2] && k-1<=object1E[l][2] ){
				numFluidNeighbors--;
				AMatrix.plusK(i,j,k-1) = 0;
			}
			if(i>=object1S[l][0] && i<=object1E[l][0] && j>=object1S[l][1] && j<=object1E[l][1] && k>=object1S[l][2] && k<=object1E[l][2] ){
				numFluidNeighbors=0;
				AMatrix.plusI(i,j,k) = 0;
				AMatrix.plusJ(i,j,k) = 0;
				AMatrix.plusK(i,j,k) = 0;
				AMatrix.plusI(i-1,j,k) = 0;
				AMatrix.plusJ(i,j-1,k) = 0;
				AMatrix.plusK(i,j,k-1) = 0;
			}
			}
			AMatrix.diag(i,j,k) = numFluidNeighbors;//;
	}
}







/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	// For now, just bypass the preconditioner:
	z = r;

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// TODO: Apply a preconditioner here.
		// For now, just bypass the preconditioner:
		z = r;

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	//PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}




/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}





/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c)
{   
   drawObjects(c);
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
   
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
	double value2 = getVelocityY(getCenter(i,j,k));
	double value3 = getVelocityX(getCenter(i,j,k));
	double hot = mT(i,j,k);
    //return vec4(0.0, 1.4*value2, 1.1*hot, value);//Blue
	return vec4(1.5*hot, 0.5*value2, 0.3, value);//Red
	//return vec4(0.1, abs(0.3*value3*value2*value)+.1, 0.1, value);//Green

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
	double value2 = getVelocityY(pt);
	double value3 = getVelocityX(pt);
	double hot = getTemperature(pt);
    //return vec4(0.0, 1.4*value2, 1.1*hot, value);//Blue
	return vec4(1.5*hot, 0.5*value2, 0.3, value);//Red
	//return vec4(0.1,abs(0.3*value3*value2*value)+.1,0.1, value);//Green

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom
	
   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
