
#include <sstream>
#include <istream>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include "LBMConstants.h"
#include "Node.h"
#include "InputVariables.h"


using namespace std;

////Simulation Parameters
//double ObstR;//20.0;
//double ObstX;//10*ObstR;//40.0;
//double ObstY;//200*0.5;
//int xDim;//200;
//int yDim;//500;//1021; //(8.2*ObstR)+2; //43;
//double CharLength;//ObstR*2.0;
//bool D2Q9i;//true;
//bool InitCond;//1;
//int tMax;//10;
//double uMax;//0.08;
//double Re;//100;
//double Tolerance;//1.0e-6;


////Compute Omega
//double nu;// = uMax*(CharLength) / Re;
//double omega1;// = 1.0/(3.0*nu+0.5);
//std::vector<double> omega;


//Global vars
int ErrorFlag;


////LBM CONSTANTS

//D3Q19 const
double t [19] = {1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,
				1/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,
				1/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
int v [19][3] = {0,0,0,1,0,0,0,1,0,-1,0,0,0,-1,0,1,1,0,-1,1,0,-1,-1,0,1,-1,0,
				0,0,1, 1,0,1, 0,1,1, -1,0,1, 0,-1,1, 
				0,0,-1, 1,0,-1, 0,1,-1, -1,0,-1, 0,-1,-1,}; 
int opposite [19] = {0,3,4,1,2,7,8,5,6,
					14, 17, 18, 15, 16,
					9, 12, 13, 10, 11};
int ysym [19] = {0,1,4,3,2,8,7,6,5,
				9, 10, 13, 12, 11,
				14, 15, 18, 17, 16};
int xsym [19] = {0,3,2,1,4,6,5,8,7,
				9, 12, 11, 10, 13,
				14, 17, 16, 15, 18};
int zsym [19] = {0,1,2,3,4,5,6,7,8,
				14, 15, 16, 17, 18,
				9, 10, 11, 12, 13};

int halfstream [9] = {1,2,5,6,
						9,10,11,12,13};

double M [19][19] = {   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
					  -30,-11,-11,-11,-11,  8,  8,  8,  8,-11,  8,  8,  8,  8,-11,  8,  8,  8,  8,
					   12, -4, -4, -4, -4,  1,  1,  1,  1, -4,  1,  1,  1,  1, -4,  1,  1,  1,  1,
					    0,  1,  0, -1,  0,  1, -1, -1,  1,  0,  1,  0, -1,  0,  0,  1,  0, -1,  0,	//u0
					    0, -4,  0,  4,  0,  1, -1, -1,  1,  0,  1,  0, -1,  0,  0,  1,  0, -1,  0,
						0,  0,  1,  0, -1,  1,  1, -1, -1,  0,  0,  1,  0, -1,  0,  0,  1,  0, -1,	//u1
						0,  0, -4,  0,  4,  1,  1, -1, -1,  0,  0,  1,  0, -1,  0,  0,  1,  0, -1,
						0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1,	//u2
						0,  0,  0,  0,  0,  0,  0,  0,  0, -4,  1,  1,  1,  1,  4, -1, -1, -1, -1,
						0,  2, -1,  2, -1,  1,  1,  1,  1, -1,  1, -2,  1, -2, -1,  1, -2,  1, -2,
						0, -4,  2, -4,  2,  1,  1,  1,  1,  2,  1, -2,  1, -2,  2,  1, -2,  1, -2,
						0,  0,  1,  0,  1,  1,  1,  1,  1, -1, -1,  0, -1,  0, -1, -1,  0, -1,  0,
						0,  0, -2,  0, -2,  1,  1,  1,  1,  2, -1,  0, -1,  0,  2, -1,  0, -1,  0,
						0,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
						0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0, -1,  0,  1,
						0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0, -1,  0,  1,  0,  
						0,  0,  0,  0,  0,  1, -1, -1,  1,  0, -1,  0,  1,  0,  0, -1,  0,  1,  0,  
						0,  0,  0,  0,  0, -1, -1,  1,  1,  0,  0,  1,  0, -1,  0,  0,  1,  0, -1,  
						0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  0, -1,  1, -1,  1};  
double M_inv [19][19];
double psi_mag2 [19]; //square of the magnitude of vectors psi, constituting M

string casename = "test2_med";

double Spacing = 0.0;
double ObstR = 2.5;//5;//5.0;//10.0;//.0;
double ObstX = 50.5;//25;//70;//13.0*2.0*ObstR;
double ObstX2 = ObstX+Spacing*2.0*ObstR;//58;//70;//13.0*2.0*ObstR;
double ObstY = 2;//2;//50.5;//51.0*2.0*ObstR/2;//1022/2;
double ObstZ = 20.5;//100;//35.5;//70.5;//50.5;//51.0*2.0*ObstR/2;//1022/2;
int xDim = ObstX+15*2*ObstR+1+Spacing*2.0*ObstR;
int yDim = 1;//1;//4;//61;//11;//102;//281; //(8.2*ObstR)+2;
int zDim = 41;//71;//141;//102;//281; //(8.2*ObstR)+2;

int maxlevel = 2;

double Re = 100;//21400;
double uMax = 0.08;//1;
double rhoIn = 1.0;
double CharLength = ObstR*2.0;
double dPdx = 0;//5.56e-7;
double Tolerance = 1.0e-6;

bool D2Q9i = true;
bool MRT = 1;//1;//false;//true;
bool LES = 0;//true;//true;
bool AMR = 0;//true;//true;
bool CB = 1;//true;//true; //Curved boundary treatment
int  ForceB1 = CB1;//boundary for force calculations
int  ForceB2 = CB2;//boundary for force calculations

bool SaveState = 1;
bool Periodics = 0;
bool xPeriodics = 0;
bool yPeriodics = 1;
bool zPeriodics = 0;

bool InitPerturb = 0;
bool InitCond = 0;//false;//true;
//string initfile ("CC6_init1.state");
string initfile (casename+".state");

int tMax = 1000;// 20000;//10000;//10000;

int StartAvg =1000;
int StartRec =200000;//30;// 30000;//start recording velocity history
int StartF = 0;//5000;//start recording force history
int StartLR = 1;//start LR

int nCPU = 2;

double s1  = 1.0;//1.19;
double s2  = 1.0;//1.4;
double s4  = 1.0;//1.2;
double s10 = 1.0;//1.4;
double s16 = 1.0;//1.98;
bool scaleS = 0;//1; //scale s's for LR



int mp1x = xDim/4;
int mp1y = yDim/4;
int mp1z = zDim/4;
int mp2x = xDim-1-xDim/4;
int mp2y = yDim/4;
int mp2z = zDim/4;
int mp3x = xDim/4;
int mp3y = yDim-1-yDim/4;
int mp3z = zDim-1-zDim/4;
int mp4x = xDim-1-xDim/4;
int mp4y = yDim-1-yDim/4;
int mp4z = zDim-1-zDim/4;

vector<double> S1(0);
vector<double> S2(0);
vector<double> S4(0);
vector<double> S10(0);
vector<double> S16(0);


//Compute Omega
double nu;// = uMax*(CharLength) / Re;
double omega1;// = 1.0/(3.0*nu+0.5);
std::vector<double> omega;



vector<double> FXbuff(0);//pow(2,maxlevel));
vector<double> FYbuff(0);//pow(2,maxlevel));
vector<double> FZbuff(0);//pow(2,maxlevel));
vector<double> FX2buff(0);//pow(2,maxlevel));
vector<double> FY2buff(0);//pow(2,maxlevel));
vector<double> FZ2buff(0);//pow(2,maxlevel));



double S[19] = { 0, //s0 density
				 1.19,//.19,//s1 omega1,//1, //energy e
				 1.4,//.4,//s2 omega1,//1, //energy square, epsilon
				 0, //x momentum
				 1.2,//.2,//s4 omega1,//1, //x energy flux
				 0, //y momentum
				 1.2,//.2,//s4 omega1,//1, //y energy flux
				 0, //z momentum
				 1.2,//.2,//s4 omega1,//1, //z energy flux
				 omega1, //s9 diagonal stress (3pxx)
				 1.4,//.4, //s10 (3pixx)
				 omega1, //s9 (pww)
				 1.4,//.4, //s10 (piww)
				 omega1, //s13 (pxy)
				 omega1, //s13 (pyz)
				 omega1, //s13 (pxz)
				 1.98,//.98, //s16 (mx)
				 1.98,//.98, //s16 (my)
				 1.98};//.98}; //s16 (mz)




// global functions for image and basic interpolations

double PoisProf (double x){
	double radius = (yDim-1-1)/2.0;
	double a = 2.0*radius;
	double result = -1.0*(((1.0-(x-0.5)/radius))*((1.0-(x-0.5)/radius))-1.0);
	return (result);
}

bool CylinderImage(double x,double y,double obX,double obY,double obR){
	return ((x-obX)*(x-obX)+(y-obY)*(y-obY) <= obR*obR);
}
bool CylinderImage2(int x,int y){
	bool output;
	return ((x-ObstX)*(x-ObstX)+(y-ObstY)*(y-ObstY) <= (ObstR+5)*(ObstR+5)
			&& (x-ObstX)*(x-ObstX)+(y-ObstY)*(y-ObstY) >= (ObstR-2)*(ObstR-2));
}


int Geometry(double x,double y,double z)
{
	int image = fluid;
	if (CylinderImage(x,z,ObstX,ObstZ,ObstR) == true) image = CB1;
//	if (CylinderImage(x,z,ObstX2,ObstZ,ObstR) == true) image = CB2;
//	if(x>=ObstX-ObstR && x<=ObstX+ObstR &&
//		z>=ObstZ-ObstR && z<=ObstZ+ObstR){
//		image = BB1;
//	}
	if (z < 0.5) image = zsymmetry;
	if (z > zDim-1.5) image = zsymmetry;
	if (x < 0.1) image = FullyDevelopedIn;
	if (x > xDim-1.1) image = FullyDevelopedOut;
	return (image);
}

double HermitSpline(double f1,double f2,double f3,double f4)
{
	double fDir1 = (f3-f1)/2.0;
	double fDir2 = (f4-f2)/2.0;
	double s = 0.5;
	double result = (2.0*s*s*s-3.0*s*s+1.0)*f2+(s*s*s-2.0*s*s+s)*fDir1+
					(-2.0*s*s*s+3.0*s*s)*f3+(s*s*s-s*s)*fDir2;
	return (result);
}

void CountElements (Node * Point,int TargetLevel,int& nodes, int& elements){
	Node * Kid;
	for(int k = 0; k<8; k++){
		Kid = Point->child[k];
		if(Kid->level == TargetLevel){
			nodes++;
			if(Kid->nb[1] != NULL &&
				Kid->nb[2] != NULL &&
				Kid->nb[5] != NULL &&
				Kid->nb[9] != NULL &&
				Kid->nb[10] != NULL &&
				Kid->nb[11] != NULL &&
				(xPeriodics == 0 || (Kid->xcoord<(xDim-1)+pow(2,-Kid->level)))&&
				(yPeriodics == 0 || (Kid->ycoord<(yDim-1)+pow(2,-Kid->level)))
//				(Kid->image >= 0 || (Kid->xcoord<(xDim-1)))&&
//				(Kid->image >= 0 || (Kid->ycoord<(yDim-1)))
				//(Kid->image >= 0 || Kid->xcoord<1)
				){
				//if(Kid->nb[5]->nb[9] != NULL){//26
				//if(Kid->nb[2]->nb[10] != NULL){//24
				//if(Kid->nb[1]->nb[11] != NULL){//26
				if(Kid->nb[9]->nb[5] != NULL){//26
					elements++;
				}
			}
		}
		else if (Kid->child[0] != NULL){
			CountElements(Kid,TargetLevel,nodes,elements);
		}
		else{
			cout<<"something wrong\n";
			return;
		}
	}
}
void WriteElements (ofstream& output,Node * Point,int TargetLevel){
	Node * Kid;
	for(int k = 0; k<8; k++){
		Kid = Point->child[k];
		if(Kid->level == TargetLevel){
			if(Kid->nb[1] != NULL &&
				Kid->nb[2] != NULL &&
				Kid->nb[5] != NULL &&
				Kid->nb[9] != NULL &&
				Kid->nb[10] != NULL &&
				Kid->nb[11] != NULL &&
				(xPeriodics == 0 || (Kid->xcoord<(xDim-1)+pow(2,-Kid->level)))&&
				(yPeriodics == 0 || (Kid->ycoord<(yDim-1)+pow(2,-Kid->level)))
				//(Kid->image >= 0 || Kid->ycoord<1) &&
				//(Kid->image >= 0 || Kid->xcoord<1)
				){
				if(Kid->nb[5]->nb[9] != NULL){
				output<<Kid->outnumb<<" ";
				output<<Kid->nb[1]->outnumb<<" ";
				output<<Kid->nb[5]->outnumb<<" ";
				output<<Kid->nb[2]->outnumb<<" ";
				output<<Kid->nb[9]->outnumb<<" ";
				output<<Kid->nb[10]->outnumb<<" ";
				output<<Kid->nb[5]->nb[9]->outnumb<<" ";
				output<<Kid->nb[11]->outnumb<<" ";
				output<<"\n";
				}
			}
//			if(Kid->nb[1] != NULL &&
//				Kid->nb[2] != NULL &&
//				Kid->nb[5] != NULL){
//				output<<Kid->outnumb<<" ";
//				output<<Kid->nb[1]->outnumb<<" ";
//				output<<Kid->nb[5]->outnumb<<" ";
//				output<<Kid->nb[2]->outnumb<<" ";
//				output<<"\n";
//			}
		}
		else if (Kid->child[0] != NULL){
			WriteElements(output,Kid,TargetLevel);
		}
		else{
			cout<<"something wrong\n";
			return;
		}
	}
}
void WriteData (ofstream& output,Node * Point,int TargetLevel,int& numb){
	Node * Kid;
	for(int k = 0; k<8; k++){
		Kid = Point->child[k];
		if(Kid->level == TargetLevel){
			output<<Kid->xcoord<<", ";
			output<<Kid->ycoord<<", ";
			output<<Kid->zcoord<<", ";
			output<<Kid->u0<<", ";
			output<<Kid->u1<<", ";
			output<<Kid->u2<<", ";
			output<<Kid->rho<<","<<Kid->u0avg<<","<<Kid->u1avg<<","<<Kid->u2avg;
			output<<","<<0<<","<<0<<","<<0<<","<<Kid->Smag;
			output<<",\n";
			Kid->outnumb = numb;
			numb++;
		}
		else if (Kid->child[0] != NULL){
			WriteData(output,Kid,TargetLevel,numb);
		}
		else{
			cout<<"something wrong\n";
			return;
		}
	}
}
void uXHist (ofstream& u_output, ofstream& v_output, Node * Point,int TargetLevel){
	Node * Kid;
	for(int k = 7; k>5; k--){//only want to output child [3] and [2]
		Kid = Point->child[k];
		if(Kid->level == TargetLevel){
			u_output<<Kid->u0<<", ";
			v_output<<Kid->u2<<", ";
		}
		else if (Kid->child[0] != NULL){
			uXHist(u_output,v_output,Kid,TargetLevel);
		}
		else{
			cout<<"something wrong\n";
			return;
		}
	}
}
void uYHist (ofstream& u_output, ofstream& v_output, Node * Point,int TargetLevel){
	Node * Kid;
	for(int k = 5; k<7; k++){//only want to output child [3] and [2]
		Kid = Point->child[k];
		if(Kid->level == TargetLevel){
			u_output<<Kid->u0<<", ";
			v_output<<Kid->u1<<", ";
		}
		else if (Kid->child[0] != NULL){
			uYHist(u_output,v_output,Kid,TargetLevel);
		}
		else{
			cout<<"something wrong\n";
			return;
		}
	}
}
void uZHist (ofstream& u_output, ofstream& w_output,
				Node * Point,int TargetLevel){
	Node * Kid;
	for(int k = 2; k<7; k+=4){//only want to output child [3] and [2]
		Kid = Point->child[k];
		if(Kid->level == TargetLevel){
			u_output<<Kid->u0<<", ";
			w_output<<Kid->u2<<", ";
		}
		else if (Kid->child[0] != NULL){
			uZHist(u_output,w_output,Kid,TargetLevel);
		}
		else{
			cout<<"something wrong\n";
			return;
		}
	}
}



void CoarsenLoop (Node * Point, int TargetLevel){
/*
Input "Point" should be the coarsest level grid (Grid[i][j]), and target level is the resulting mesh level. Function checks the level of the "current" grid, and if it is less than the target level, it updates the pointer to Point->child[k] to look at the next refinement level. This is continued until the grid level of the pointer is equal to the target level. At this point, the function will delete all children nodes held by the current node, and ensures the children are NULL.
*/
	Node * Kid;
	Kid = Point;
		if(Kid->level < TargetLevel && Kid->child[0] != NULL){
		for(int k = 0; k<8; k++){
			Kid = Point->child[k];
			CoarsenLoop(Kid,TargetLevel);
			}
		}
		else if(Kid->level == TargetLevel){
			Kid->DeleteNode();
			for(int k = 0; k<8; k++)
				Kid->child[k] = NULL;
		}
}
//void NeighborLoop (Node * Point){//, int TargetLevel){
///*
//If Point has children that are not NULL, ChildNeighbor is used on Point, to define the neighbor relations of its children. Then, each child in Point are checked, until the NULL child is found.
//*/
//	Node * Kid;
//	Kid = Point->child[0];
//	if(Kid != NULL){
//		Point->ChildNeighbor();
//		for(int k = 0; k<4; k++){
//			Kid = Point->child[k];
//			NeighborLoop(Kid);
//		}
//	}
//	else if(Kid == NULL){
//		return;
//	}
//	else
//		cout<<"something wrong in NeighborLoop\n";
//}
//
void OmegaValues(){
	if(omega1 < 0.5 || omega1 > 2.0)
		cout<<"omega is out of acceptable range\n";
	while(omega1 > 0.50 && omega1 < 2.00) {
		omega.push_back(omega1);
		S1 .push_back(s1 );
		S2 .push_back(s2 );
		S4 .push_back(s4 );
		S10.push_back(s10);
		S16.push_back(s16);
		omega1 = 2.0/(1.0+2.0*(2.0/omega1-1.0));
		if(scaleS == true){
			s1  = 2.0/(1.0+2.0*(2.0/s1 -1.0));
			s2  = 2.0/(1.0+2.0*(2.0/s2 -1.0));
			s4  = 2.0/(1.0+2.0*(2.0/s4 -1.0));
			s10 = 2.0/(1.0+2.0*(2.0/s10-1.0));
			s16 = 2.0/(1.0+2.0*(2.0/s16-1.0));
			if(s1  < 1.0) s1  = 1.0;
			if(s2  < 1.0) s2  = 1.0;
			if(s4  < 1.0) s4  = 1.0;
			if(s10 < 1.0) s10 = 1.0;
			if(s16 < 1.0) s16 = 1.0;
		}
		if(omega.size() > 20){
			cout<<"omega error\n";
			break;
		}
	}
}

//void MarchLoop1(Node * Point,int GridLevel,int Target)//,void Node::function())
void MarchLoop1(Node * Point,int Target)//,void Node::function())
{
	Node * Kid;
	if(Point->level == Target) // && Point->edgeP == 0)
		Point->Stream1();
	else if (Point->child[0] != NULL)//(GridLevel > Point->level) //
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			MarchLoop1(Kid,Target);//,function);
		}
	}
} 

void MarchLoop2(Node * Point,int Target)
{
	Node * Kid;
	if(Point->level == Target){
		Point->Stream2();
		Point->BoundaryMacro();
//		Point->InletOutlet();
		Point->RegularizedBoundary();
		Point->Symmetry();
		Point->ComputeMacros();
		Point->ComputeFeq();
		Point->Collide();
		Point->BounceBack();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			MarchLoop2(Kid,Target);
		}
	}
} 
void MarchLoop3(Node * Point,int Target)
{
	Node * Kid;
	if(Point->level == Target){
		if(CB == true) Point->CurvedB();
//		Point->Symmetry();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			MarchLoop3(Kid,Target);
		}
	}
} 
void RefineLoop(Node * Point,int GridLevel,int Target)
{
/*
Target is the level of the parent. If grid is to be refined to make a level 2 mesh, the target would be 1.
*/
	Node * Kid;
	if(Point->level == Target-1 && GridLevel >= Target){
		if(Point->child[0] == NULL){
			Point->Refine();
		}
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			RefineLoop(Kid,GridLevel,Target);
		}
	}
} 

void SpatialInterpLoop1(Node * Point,int GridLevel,int Target)
{
/*
Interpolate f values to obtain child f values. 
*/
	Node * Kid;
	if(Point->level == Target-1 && GridLevel >= Target){
//		for(int k = 0;k<4;k++)
//			Point->child[k]->edge = 1;
		Point->SpatialInterp1();
//		for(int i = 0;i<4;i++)
//			Point->child[i]->ComputeMacros();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			SpatialInterpLoop1(Kid,GridLevel,Target);
		}
	}
} 
//void SpatialInterpLoop3(Node * Point,int GridLevel,int Target)
//{
///*
//Interpolate f values to obtain child f values. 
//*/
//	Node * Kid;
//	if(Point->level == Target-1 && GridLevel >= Target){
////		for(int k = 0;k<4;k++)
////			Point->child[k]->edge = 1;
//		Point->SpatialInterp1();
////		for(int i = 0;i<4;i++)
////			Point->child[i]->ComputeMacros();
//	}
//	else if (Point->child[0] != NULL)
//	{
//		for(int k = 0;k<4;k++){
//			Kid = Point->child[k];
//			SpatialInterpLoop3(Kid,GridLevel,Target);
//		}
//	}
//} 
void SpatialInterpLoop2(Node * Point,int GridLevel,int Target)
{
/*
Interpolate f values to obtain child f values at 2 timesteps ahead. Store in ftemp
*/
	Node * Kid;
	if(Point->level == Target-1 && GridLevel >= Target){
		Point->SpatialInterp2();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			SpatialInterpLoop2(Kid,GridLevel,Target);
		}
	}
} 
void SpatialAverageLoop(Node * Point,int GridLevel,int Target)
{
/*
Extract f values from children
*/
	Node * Kid;
	if(Point->level == Target-1 && GridLevel >= Target){
		Point->SpatialAverage();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			SpatialAverageLoop(Kid,GridLevel,Target);
		}
	}
} 
void TemporalInterpLoop(Node * Point,int GridLevel,int Target)
{
/*
Update ftemp so its value is the interpolated value between two timesteps. Parent level
*/
	Node * Kid;
	if(Point->level == Target-1 && GridLevel >= Target){
		Point->TemporalInterp();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			TemporalInterpLoop(Kid,GridLevel,Target);
		}
	}
} 
void TemporalInterpLoop2(Node * Point,int GridLevel,int Target)
{
/*
Copy ftemp to f, for edge nodes
*/
	Node * Kid;
	if(Point->level == Target-1 && GridLevel >= Target){
		Point->TemporalInterp2();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			TemporalInterpLoop2(Kid,GridLevel,Target);
		}
	}
} 

void ChildNeighborLoop(Node * Point,int GridLevel,int Target)
{
/*
*/
	Node * Kid;
	if(Point->level == Target-1 && GridLevel >= Target){
		Point->ChildNeighbor();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			ChildNeighborLoop(Kid,GridLevel,Target);
		}
	}
} 
void MacroLoop(Node * Point,int Target)
{
	Node * Kid;
	if(Point->level == Target){
		Point->ComputeMacros();
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			MacroLoop(Kid,Target);
		}
	}
} 
void DeltaLoop(Node * Point,int Target)
{
	Node * Kid;
	if(Point->level == Target){
		if((Point->image == CB1 || Point->image == CB2) && Point->edge == 0){
		Point->ComputeDelta(ObstX,ObstZ,ObstX2,ObstZ,ObstR);
		}
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			DeltaLoop(Kid,Target);
		}
	}
} 
void BBForceLoop(Node * Point,int Target,double& ForceX, double& ForceY, double& ForceZ, int BBimage)
{
	Node * Kid;
	if(Point->level == Target){
		Point->BBForce(ForceX,ForceY,ForceZ,BBimage);
	}
	else if (Point->child[0] != NULL)
	{
		for(int k = 0;k<8;k++){
			Kid = Point->child[k];
			BBForceLoop(Kid,Target,ForceX,ForceY,ForceZ,BBimage);
		}
	}
} 







void LRLBM(Node **** Grid,int *** GridLevel,int Target,
						int minx,int miny,int minz, int maxx,int maxy,int maxz,
						int t,vector<ofstream*> u_output1, vector<ofstream*> v_output1,
						vector<ofstream*> u_output2, vector<ofstream*> v_output2)
{
	bool check = 0;
	int i,j,k,MinX,MinY,MinZ,MaxX,MaxY,MaxZ;
	MinY = yDim;
	MinX = xDim;
	MinZ = zDim;
	MaxY = 0;
	MaxX = 0;
	MaxZ = 0;
	#pragma omp parallel 
	{
	//interpolate target mesh to get finer mesh values at edges
	#pragma omp for collapse(3) 
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] > Target){
			SpatialInterpLoop1(Grid[i][j][k],GridLevel[i][j][k],Target+1);
			check = 1;
			#pragma omp critical
			{
			if(i<MinX) MinX = i;
			if(j<MinY) MinY = j;
			if(k<MinZ) MinZ = k;
			if(i+1>MaxX) MaxX = i+1;
			if(j+1>MaxY) MaxY = j+1;
			if(k+1>MaxZ) MaxZ = k+1;
			}
		}
	}}}
	//March 1 for target mesh
	#pragma omp for collapse(3) 
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] >= Target)
			MarchLoop1(Grid[i][j][k],Target);
	}}}
	//March 2 for target mesh
	#pragma omp for collapse(3) 
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] >= Target)
			MarchLoop2(Grid[i][j][k],Target);
	}}}
	//March 3 for target mesh
	#pragma omp for collapse(3) 
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] >= Target)
			MarchLoop3(Grid[i][j][k],Target);
	}}}
	//obtain f from target mesh and store in ftemp of finer mesh
	#pragma omp for collapse(3) 
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] > Target)
			SpatialInterpLoop2(Grid[i][j][k],GridLevel[i][j][k],Target+1);
	}}}
	//interpolate current f with ftemp, and store in ftemp
	#pragma omp for collapse(3) 
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] > Target) //changed to >= from >
			TemporalInterpLoop(Grid[i][j][k],GridLevel[i][j][k],Target+1); 
	}}}
	}//end of omp parallel
	
	//force calculations, store to buffer
	if(Target == maxlevel && t >= StartF){
		double Fx = 0;
		double Fy = 0;
		double Fz = 0;
		double Fx2 = 0;
		double Fy2 = 0;
		double Fz2 = 0;
		for(i = minx;i<maxx;i++){
		for(j = miny;j<maxy;j++){
		for(k = minz;k<maxz;k++){
			if(GridLevel[i][j][k] >= Target){ //changed to >= from >
				BBForceLoop(Grid[i][j][k],Target,Fx,Fy,Fz,ForceB1);
				BBForceLoop(Grid[i][j][k],Target,Fx2,Fy2,Fz2,ForceB2);
			}
		}}}
		FXbuff.push_back(Fx);
		FYbuff.push_back(Fy);
		FZbuff.push_back(Fz);
		FX2buff.push_back(Fx2);
		FY2buff.push_back(Fy2);
		FZ2buff.push_back(Fz2);
	}

	//write u v to file (different file for each refinement level)
	if(t >= StartRec){
	*u_output1[Target]<<t<<",";
	*v_output1[Target]<<t<<",";
	*u_output2[Target]<<t<<",";
	*v_output2[Target]<<t<<",";
	for(i = minx;i<maxx;i++){
		if(GridLevel[i][int(ObstY)][int(ObstZ)] >= Target)
			uXHist(*u_output1[Target],*v_output1[Target],Grid[i][int(ObstY)][int(ObstZ)],Target); 
	}
	for(i = miny;i<maxy;i++){
		if(GridLevel[int(ObstX)][i][int(ObstZ)] >= Target)
			uYHist(*u_output2[Target],*v_output2[Target],Grid[int(ObstX)][i][int(ObstZ)],Target); 
	}
	*u_output1[Target]<<endl;
	*v_output1[Target]<<endl;
	*u_output2[Target]<<endl;
	*v_output2[Target]<<endl;
	}



	if(check == 1){
		LRLBM(Grid,GridLevel,Target+1,MinX,MinY,MinZ,MaxX,MaxY,MaxZ,t,u_output1,v_output1,u_output2,v_output2);
	}
	#pragma omp parallel 
	{
	//apply ftemp to f on edges of current mesh
	#pragma omp for collapse(3)
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] >= Target)
			TemporalInterpLoop2(Grid[i][j][k],GridLevel[i][j][k],Target);
	}}}
	#pragma omp for collapse(3)
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] > Target){
			SpatialInterpLoop1(Grid[i][j][k],GridLevel[i][j][k],Target+1);
			check = 1;
		}
	}}}
	#pragma omp for collapse(3)
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] >= Target)
			MarchLoop1(Grid[i][j][k],Target);
	}}}
	#pragma omp for collapse(3)
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] >= Target)
			MarchLoop2(Grid[i][j][k],Target);
	}}}
	#pragma omp for collapse(3)
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] >= Target)
			MarchLoop3(Grid[i][j][k],Target);
	}}}

	//obtain f from target mesh and store in ftemp of finer mesh
	#pragma omp for collapse(3)
	for(i = minx; i<maxx; i++){
	for(j = miny; j<maxy; j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] > Target)
			SpatialInterpLoop2(Grid[i][j][k],GridLevel[i][j][k],Target+1);//1 to 2
	}}}
	//interpolate current f with ftemp, and store in ftemp
	#pragma omp for collapse(3)
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] > Target) //changed to >= from >
			TemporalInterpLoop(Grid[i][j][k],GridLevel[i][j][k],Target+1); 
	}}}

	}//end of omp parallel


	//force calculations, store to buffer
	if(Target == maxlevel && t >= StartF){
		double Fx = 0;
		double Fy = 0;
		double Fz = 0;
		double Fx2 = 0;
		double Fy2 = 0;
		double Fz2 = 0;
		for(i = minx;i<maxx;i++){
		for(j = miny;j<maxy;j++){
		for(k = minz;k<maxz;k++){
			if(GridLevel[i][j][k] >= Target){ //changed to >= from >
				BBForceLoop(Grid[i][j][k],Target,Fx,Fy,Fz,ForceB1);
				BBForceLoop(Grid[i][j][k],Target,Fx2,Fy2,Fz2,ForceB2);
			}
		}}}
		FXbuff.push_back(Fx);
		FYbuff.push_back(Fy);
		FZbuff.push_back(Fz);
		FX2buff.push_back(Fx2);
		FY2buff.push_back(Fy2);
		FZ2buff.push_back(Fz2);
	}

//	//write u v to file
//	if(t >= StartRec){
//	*u_output[Target]<<t+(1.0/pow(2.0,Target))<<",";
//	*v_output[Target]<<t+(1.0/pow(2.0,Target))<<",";
//	for(i = minx;i<maxx;i++){
//		if(GridLevel[i][int(ObstY)] >= Target)
//			uvHist(*u_output[Target],*v_output[Target],Grid[i][int(ObstY)][int(ObstZ)],Target); 
//	}
//	*u_output[Target]<<endl;
//	*v_output[Target]<<endl;
//	}
	//write u v to file (different file for each refinement level)
	if(t >= StartRec){
	*u_output1[Target]<<t+(1.0/pow(2.0,Target))<<",";
	*v_output1[Target]<<t+(1.0/pow(2.0,Target))<<",";
	*u_output2[Target]<<t+(1.0/pow(2.0,Target))<<",";
	*v_output2[Target]<<t+(1.0/pow(2.0,Target))<<",";
	for(i = minx;i<maxx;i++){
		if(GridLevel[i][int(ObstY)][int(ObstZ)] >= Target)
			uXHist(*u_output1[Target],*v_output1[Target],Grid[i][int(ObstY)][int(ObstZ)],Target); 
	}
//	for(i = miny;i<maxy;i++){
//		if(GridLevel[xDim/2][int(ObstY)][i] >= Target)
//			uZHist(*u_output2[Target],*v_output2[Target],Grid[xDim/2][int(ObstY)][i],Target); 
//	}
	for(i = miny;i<maxy;i++){
		if(GridLevel[int(ObstX)][i][int(ObstZ)] >= Target)
			uYHist(*u_output2[Target],*v_output2[Target],Grid[int(ObstX)][i][int(ObstZ)],Target); 
	}
	*u_output1[Target]<<endl;
	*v_output1[Target]<<endl;
	*u_output2[Target]<<endl;
	*v_output2[Target]<<endl;
	}
	if(check == 1){
		LRLBM(Grid,GridLevel,Target+1,MinX,MinY,MinZ,MaxX,MaxY,MaxZ,t+(1.0/pow(2.0,Target)),u_output1,v_output1,u_output2,v_output2);
	}
	//extract finer mesh data and apply to current mesh	
	#pragma omp parallel for collapse(3)
	for(i = minx;i<maxx;i++){
	for(j = miny;j<maxy;j++){
	for(k = minz;k<maxz;k++){
		if(GridLevel[i][j][k] > Target){
			SpatialInterpLoop1(Grid[i][j][k],GridLevel[i][j][k],Target+1);
			SpatialAverageLoop(Grid[i][j][k],GridLevel[i][j][k],Target+1);
			MacroLoop(Grid[i][j][k],Target+1);
		}
	}}}
}

void Remap(int *** GridLevel,Node **** Grid,
					double Crit[],int Max)
{
	int i,j,k,n,m,a,b,c,nbBoundary,l;
	for(i = 0; i<xDim; i++){
	for(j = 0; j<yDim; j++){
	for(k = 0; k<zDim; k++){
		GridLevel[i][j][k] = 0;
	}}}
	for(n = 0; n<Max; n++){
//		#pragma omp parallel for collapse(2)
		for(i = 1; i<xDim-1; i++){
		for(j = 1; j<yDim-1; j++){
		for(k = 1; k<zDim-1; k++){
			nbBoundary = 0;
			for(l = 1; l<19; l++) 
				nbBoundary += Grid[i][j][k]->nb[l]->image;//refine nodes next to boundaries
			if(((Grid[i][j][k]->Smag+Grid[i][j][k]->image+nbBoundary) > Crit[Max-n] && GridLevel[i][j][k] < Max-n &&
				i >= Max-n && i <= xDim-1-(Max-n) &&
				j >= Max-n && j <= yDim-1-(Max-n) &&
				k >= Max-n && k <= zDim-1-(Max-n))) { //|| 
				GridLevel[i][j][k] = Max-n;
				for(m = (Max-n-1);m>0;m--){
				for(a = i-m;a<=i+m;a++){
				for(b = j-m;b<=j+m;b++){
				for(c = k-m;c<=k+m;c++){
					if(GridLevel[a][b][c] < Max-n-m)
					GridLevel[a][b][c] = Max-n-m;
				}}}}
			}
		}}}
	}	
}
void vort_iCalc(Node **** Grid,double *** vort)
{
	int i,j,k;
//	#pragma omp parallel for collapse(2)
	for(i = 1; i<xDim-1; i++){
	for(j = 1; j<yDim-1; j++){
	for(k = 1; k<zDim-1; k++){
		if (Grid[i][j][k]->image == fluid){
			vort[i][j][k] = (Grid[i][j+1][k]->u2-Grid[i][j-1][k]->u2+
						Grid[i][j][k+1]->u1-Grid[i][j][k-1]->u1)/2.0;
		}
		else{
			vort[i][j][k] = 0;
		}
	}}}
}
void vort_jCalc(Node **** Grid,double *** vort)
{
	int i,j,k;
//	#pragma omp parallel for collapse(2)
	for(i = 1; i<xDim-1; i++){
	for(j = 1; j<yDim-1; j++){
	for(k = 1; k<zDim-1; k++){
		if (Grid[i][j][k]->image == fluid){
			vort[i][j][k] = (Grid[i+1][j][k]->u2-Grid[i-1][j][k]->u2+
						Grid[i][j][k+1]->u0-Grid[i][j][k-1]->u0)/2.0;
		}
		else{
			vort[i][j][k] = 0;
		}
	}}}
}
void vort_kCalc(Node **** Grid,double *** vort)
{
	int i,j,k;
//	#pragma omp parallel for collapse(2)
	for(i = 1; i<xDim-1; i++){
	for(j = 1; j<yDim-1; j++){
	for(k = 1; k<zDim-1; k++){
		if (Grid[i][j][k]->image == fluid){
			int check = 0;
//			if (CylinderImage2(i,j) == true)
//				check = 1; 
			if(check == 0){
				vort[i][j][k] = (Grid[i+1][j][k]->u1-Grid[i-1][j][k]->u1+
							Grid[i][j+1][k]->u0-Grid[i][j-1][k]->u0)/2.0;
			}
			else
				vort[i][j][k] = 0.05;
		}
		else{
			int check = 0;
//			if (CylinderImage2(i,j) == true)
//				vort[i][j][k] = 0.05;
//			else
				vort[i][j][k] = 0;
		}
	}}}
}
void vort_iCalc_Avg(Node **** Grid,double *** vort)
{
	int i,j,k;
//	#pragma omp parallel for collapse(2)
	for(i = 1; i<xDim-1; i++){
	for(j = 1; j<yDim-1; j++){
	for(k = 1; k<zDim-1; k++){
		if (Grid[i][j][k]->image == fluid){
			vort[i][j][k] = (Grid[i][j+1][k]->u2avg-Grid[i][j-1][k]->u2avg+
						Grid[i][j][k+1]->u1avg-Grid[i][j][k-1]->u1avg)/2.0;
		}
		else{
			vort[i][j][k] = 0;
		}
	}}}
}
void vort_jCalc_Avg(Node **** Grid,double *** vort)
{
	int i,j,k;
//	#pragma omp parallel for collapse(2)
	for(i = 1; i<xDim-1; i++){
	for(j = 1; j<yDim-1; j++){
	for(k = 1; k<zDim-1; k++){
		if (Grid[i][j][k]->image == fluid){
			vort[i][j][k] = (Grid[i+1][j][k]->u2avg-Grid[i-1][j][k]->u2avg+
						Grid[i][j][k+1]->u0avg-Grid[i][j][k-1]->u0avg)/2.0;
		}
		else{
			vort[i][j][k] = 0;
		}
	}}}
}
void vort_kCalc_Avg(Node **** Grid,double *** vort)
{
	int i,j,k;
//	#pragma omp parallel for collapse(2)
	for(i = 1; i<xDim-1; i++){
	for(j = 1; j<yDim-1; j++){
	for(k = 1; k<zDim-1; k++){
		if (Grid[i][j][k]->image == fluid){
			vort[i][j][k] = (Grid[i+1][j][k]->u1avg-Grid[i-1][j][k]->u1avg+
						Grid[i][j+1][k]->u0avg-Grid[i][j-1][k]->u0avg)/2.0;
		}
		else{
			vort[i][j][k] = 0;
		}
	}}}
}

void SmagCalc(Node **** Grid,double *** smag)
{
	int i,j,k;
	double S11, S22, S33, S12, S13, S23;
//	#pragma omp parallel for collapse(2)
	for(i = 1; i<xDim-1; i++){
	for(j = 1; j<yDim-1; j++){
	for(k = 1; k<zDim-1; k++){
		S11 = (Grid[i+1][j][k]->u0-Grid[i-1][j][k]->u0)*0.5;
		S22 = (Grid[i][j+1][k]->u1-Grid[i][j-1][k]->u1)*0.5;
		S33 = (Grid[i][j][k+1]->u2-Grid[i][j][k-1]->u2)*0.5;
		S12 = (Grid[i+1][j][k]->u1-Grid[i-1][j][k]->u1+
				Grid[i][j+1][k]->u0-Grid[i][j-1][k]->u0)*0.25;
		S23 = (Grid[i][j+1][k]->u2-Grid[i][j-1][k]->u2+
				Grid[i][j][k+1]->u1-Grid[i][j][k-1]->u1)*0.25;
		S13 = (Grid[i][j][k+1]->u0-Grid[i][j][k-1]->u0+
				Grid[i+1][j][k]->u2-Grid[i-1][j][k]->u2)*0.25;
		smag[i][j][k] = sqrt(S11*S11+S22*S22+S33*S33+S12*S12+S13*S13+S23*S23);
	}}}
}int main ()
{
	ErrorFlag = 0;
	int ConvFlag = 0;
	// For computing force in X and Y
	double FX = 0;
	double FY = 0;
	double FZ = 0;
	double FX2 = 0;
	double FY2 = 0;
	double FZ2 = 0;


	//Compute Omega
	nu = uMax*(CharLength)/Re;
	omega1 = 1.0/(3.0*nu+0.5);


	int i,j,k,l,n,Target;
	Node **** Grid;
	Grid = new Node *** [xDim];
	for(i = 0; i<xDim; i++){
		Grid[i] = new Node ** [yDim];
		for(j = 0; j<yDim; j++){
			Grid[i][j] = new Node * [zDim];
			for(k = 0; k<zDim; k++){
				Grid[i][j][k] = new Node();
			}
		}
	}
	int *** GridLevel;
	GridLevel = new int ** [xDim];
	for(i = 0; i<xDim; i++){
		GridLevel[i] = new int * [yDim];
		for(j = 0; j<yDim; j++){
			GridLevel[i][j] = new int [zDim];
			for(k = 0; k<zDim; k++){
				GridLevel[i][j][k] = 0;
			}
		}
	}
	double *** vort_i;
	double *** vort_j;
	double *** vort_k;
	vort_i = new double ** [xDim];
	vort_j = new double ** [xDim];
	vort_k = new double ** [xDim];
	for(i = 0; i<xDim; i++){
		vort_i[i] = new double * [yDim];
		vort_j[i] = new double * [yDim];
		vort_k[i] = new double * [yDim];
		for(j = 0; j<yDim; j++){
			vort_i[i][j] = new double [zDim];
			vort_j[i][j] = new double [zDim];
			vort_k[i][j] = new double [zDim];
			for(k = 0; k<zDim; k++){
				vort_i[i][j][k] = 0;
				vort_j[i][j][k] = 0;
				vort_k[i][j][k] = 0;
			}
		}
	}

	double *** smag;
	smag = new double ** [xDim];
	for(i = 0; i<xDim; i++){
		smag[i] = new double * [yDim];
		for(j = 0; j<yDim; j++){
			smag[i][j] = new double [zDim];
			for(k = 0; k<zDim; k++){
				smag[i][j][k] = 0;
			}
		}
	}


//	ofstream Minvout;
//	Minvout.open("Minvout.dat");
	//Compute M inverse by transposing M and dividing by psi_mag2
	//Then multiply by S to get M_inv*S, and store as M_inv
	for(i = 0; i<19; i++){
		psi_mag2[i] = 0;
		for(j = 0; j<19; j++){
			M_inv[j][i] = M[i][j];
			psi_mag2[i] += M[i][j]*M[i][j];
		}
	}
	cout<<endl;
	for(i = 0; i<19; i++){
		for(j = 0; j<19; j++){
			M_inv[i][j] = M_inv[i][j]/psi_mag2[j];
//			Minvout<<", "<<setw(10)<<setprecision(8)<<M_inv[i][j];
		}
//		Minvout<<endl;
	}
//	Minvout.close();



//	double value = 0;
//	for(i = 0; i<19; i++){
//		for(j = 0; j<19; j++){
//			value = 0;
//			for(k = 0; k<19; k++){
//				value += M[i][k]*M_inv[k][j];
//			}
//			cout<<value<<", ";
//		}
//		cout<<endl;
//	}

	OmegaValues();
//	int maxlevel = omega.size()-1;
//	maxlevel = 3;
	for(i = 0; i<omega.size();i++)
		cout<<"omega"<<i<<" = "<<omega[i]<<"\n";
	for(i = 0; i<S1.size();i++){
		cout<<"S1 "<<i<<" = "<<S1[i]<<"\n";
		cout<<"S2 "<<i<<" = "<<S2 [i]<<"\n";
		cout<<"S4 "<<i<<" = "<<S4 [i]<<"\n";
		cout<<"S10"<<i<<" = "<<S10[i]<<"\n";
		cout<<"S16"<<i<<" = "<<S16[i]<<"\n";
	}
//	double Crit[4+1] = {0,0.0025,0.005,0.0075,0.01};
	double Crit[3] = {0,0.005,0.01};

	cout<<"\nCASE NAME: "<<casename<<"\n\n";
	cout<<"\nSIMULATION PARAMETERS:\n\n";
	cout<<"xDim:\t"<<xDim<<"\n";
	cout<<"yDim:\t"<<yDim<<"\n";
	cout<<"zDim:\t"<<zDim<<"\n";
	cout<<"ObstR:\t"<<ObstR<<"\n";
	cout<<"ObstX:\t"<<ObstX<<"\n";
	cout<<"ObstY:\t"<<ObstY<<"\n";
	cout<<"ObstZ:\t"<<ObstZ<<"\n";
	cout<<"maxlevel:\t"<<maxlevel<<"\n\n";

	cout<<"Re:\t"<<Re<<"\n";
	cout<<"uMax:\t"<<uMax<<"\n";
	cout<<"CharLength:\t"<<CharLength<<"\n";
	cout<<"dp/dx:\t"<<dPdx<<"\n";
	cout<<"Tolerance:\t"<<Tolerance<<"\n\n";

	cout<<"D2Q9i:\t"<<D2Q9i<<"\n";
	cout<<"MRT:\t"<<MRT<<"\n";
	cout<<"LES:\t"<<LES<<"\n";
	cout<<"SaveState:\t"<<SaveState<<"\n";
	cout<<"Periodics:\t"<<Periodics<<"\n";
	cout<<"InitPerturb:\t"<<InitPerturb<<"\n";
	cout<<"InitCond:\t"<<InitCond<<"\n";
	if(InitCond == true){
		cout<<"InitFile:\t"<<initfile<<"\n\n";
	}
	else
		cout<<endl;

	cout<<"tMax:\t"<<tMax<<"\n";
	cout<<"StartAvg:\t"<<StartAvg<<"\n";
	cout<<"StartRec:\t"<<StartRec<<"\n";
	cout<<"StartF:\t"<<StartF<<"\n\n";

	cout<<"omega:\t"<<omega[0]<<"\n\n";



	//read header lines for input file
	fstream infile;
	
	if(InitCond == true){
	infile.open(initfile.c_str());
//	char dummy[256];
//	infile.getline(dummy,256);
//	infile.getline(dummy,256);
	cout<<"Loading initial conditions file\n";
	}
	string line;

    for(i = 0; i<xDim; i++){
	for(j = 0; j<yDim; j++){
	for(k = 0; k<zDim; k++){
		Grid[i][j][k]->SetImage(Geometry(i,j,k));
		if(InitCond == true){
			//input file
	    	getline(infile, line); 
	        istringstream tokenizer(line);
	        string token;
	
	        getline(tokenizer, token, ',');
	        istringstream int_iss(token);
	        int xin;
	        int_iss >> xin;
	
	        getline(tokenizer, token, ',');
	        istringstream int_iss2(token);
	        int yin;
	        int_iss2 >> yin;
	
	        getline(tokenizer, token, ',');
	        istringstream int_iss3(token);
	        int zin;
	        int_iss3 >> zin;

			Grid[xin][yin][zin]->Initialize(0,0,0,1.0,xin,yin,zin);
			for(l = 0; l<19; l++){
		        getline(tokenizer, token, ',');
		        istringstream double_iss1(token);
		        double f_in;
		        double_iss1 >> f_in;
				Grid[xin][yin][zin]->f[l] = f_in;
			}
		}
		else
		   	Grid[i][j][k]->Initialize(uMax,0,0,1.0,i,j,k);
		GridLevel[i][j][k] = 0;
	}}}
	if(InitCond == true)
		infile.close();


	if(InitPerturb == true){
		Grid[int(xDim/8)][int(yDim/2)][int(zDim/2)]->f[11] += 0.05;
	}

//	for(i = ObstX-0.5*ObstR; i<ObstX+0.5*ObstR; i++){
//	for(j = ObstY-0.5*ObstR; j<ObstY+0.5*ObstR; j++){
//		GridLevel[i][j] = 2;
//	}}
	


//	//Define neighbor nodes for coarse grid
	int nextX,nextY,nextZ;
//	#pragma omp parallel for collapse(3)
	for(i = 0; i<xDim; i++){
	for(j = 0; j<yDim; j++){
	for(k = 0; k<zDim; k++){
//		Grid[i][j]->nb[0] = NULL;
		for (l = 1; l<19; l++){
			nextX = i+v[l][0];	
			nextY = j+v[l][1];
			nextZ = k+v[l][2];
			if ((0<=nextX) && (nextX<xDim) && (0<=nextY) && (nextY<yDim)
				&& (0<=nextZ) && (nextZ<zDim)){
				Grid[i][j][k]->nb[l] = (Grid[nextX][nextY][nextZ]);	
			}
			else
				Grid[i][j][k]->nb[l] = NULL;
		}
	}}}


	if(yPeriodics == true){
	for(i = 0; i<xDim; i++){
	for(k = 0; k<zDim; k++){
//		if(Grid[i][yDim-1][k]->image == yperiodic){
		for (l = 1; l<19; l++){
			if(Grid[i][yDim-1][k]->nb[l] == NULL){
				nextX = i+v[l][0];	
				nextZ = k+v[l][2];
//				if(nextX > xDim-1) nextX -= xDim;
				if(nextZ > zDim-1) nextZ -= zDim;
//				if(nextX < 0) nextX += xDim;
				if(nextZ < 0) nextZ += zDim;
				if ((0<=nextX) && (nextX<xDim) && (0<=nextZ) && (nextZ<zDim)){
				Grid[i][yDim-1][k]->nb[l] = Grid[nextX][0][nextZ];
				}
			}
		}
//		}
//		if(Grid[i][0][k]->image == yperiodic){
		for (l = 1; l<19; l++){
			if(Grid[i][0][k]->nb[l] == NULL){
				nextX = i+v[l][0];	
				nextZ = k+v[l][2];
//				if(nextX > xDim-1) nextX -= xDim;
				if(nextZ > zDim-1) nextZ -= zDim;
//				if(nextX < 0) nextX += xDim;
				if(nextZ < 0) nextZ += zDim;
				if ((0<=nextX) && (nextX<xDim) && (0<=nextZ) && (nextZ<zDim)){
				Grid[i][0][k]->nb[l] = Grid[nextX][yDim-1][nextZ];
				}
			}
		}
//		}
	}}
	
	}
	if(xPeriodics == true){
	for(j = 0; j<yDim; j++){
	for(k = 0; k<zDim; k++){
//		if(Grid[xDim-1][j][k]->image <0){//== xperiodic){
		for (l = 1; l<19; l++){
			if(Grid[xDim-1][j][k]->nb[l] == NULL){
				nextY = j+v[l][1];	
				nextZ = k+v[l][2];
//				if(nextY > yDim-1) nextY -= yDim;
//				if(nextZ > zDim-1) nextZ -= zDim;
//				if(nextY < 0) nextY += yDim;
//				if(nextZ < 0) nextZ += zDim;
				if ((0<=nextY) && (nextY<yDim) && (0<=nextZ) && (nextZ<zDim)){
				Grid[xDim-1][j][k]->nb[l] = Grid[0][nextY][nextZ];
				}
			}
		}
//		}
//		if(Grid[0][j][k]->image <0){//== xperiodic){
		for (l = 1; l<19; l++){
			if(Grid[0][j][k]->nb[l] == NULL){
				nextY = j+v[l][1];	
				nextZ = k+v[l][2];
//				if(nextY > yDim-1) nextY -= yDim;
//				if(nextZ > zDim-1) nextZ -= zDim;
//				if(nextY < 0) nextY += yDim;
//				if(nextZ < 0) nextZ += zDim;
				if ((0<=nextY) && (nextY<yDim) && (0<=nextZ) && (nextZ<zDim)){
				Grid[0][j][k]->nb[l] = Grid[xDim-1][nextY][nextZ];
				}
			}
		}
//		}

	}}
	}
	if(zPeriodics == true){
	for(i = 0; i<xDim; i++){
	for(j = 0; j<yDim; j++){
//		if(Grid[xDim-1][j][k]->image <0){//== xperiodic){
		for (l = 1; l<19; l++){
			if(Grid[i][j][zDim-1]->nb[l] == NULL){
				nextX = i+v[l][0];	
				nextY = j+v[l][1];	
				if(nextY > yDim-1) nextY -= yDim;
//				if(nextZ > zDim-1) nextZ -= zDim;
				if(nextY < 0) nextY += yDim;
//				if(nextZ < 0) nextZ += zDim;
				if ((0<=nextY) && (nextY<yDim) && (0<=nextX) && (nextX<xDim)){
				Grid[i][j][zDim-1]->nb[l] = Grid[nextX][nextY][0];
				}
			}
		}
//		}
//		if(Grid[0][j][k]->image <0){//== xperiodic){
		for (l = 1; l<19; l++){
			if(Grid[i][j][0]->nb[l] == NULL){
				nextX = i+v[l][0];	
				nextY = j+v[l][1];	
				if(nextY > yDim-1) nextY -= yDim;
//				if(nextZ > zDim-1) nextZ -= zDim;
				if(nextY < 0) nextY += yDim;
//				if(nextZ < 0) nextZ += zDim;
				if ((0<=nextY) && (nextY<yDim) && (0<=nextX) && (nextX<xDim)){
				Grid[i][j][0]->nb[l] = Grid[nextX][nextY][xDim-1];
				}
			}
		}
//		}

	}}
	}


	for(Target = 0; Target<maxlevel; Target++){
		#pragma omp parallel
		{
//		#pragma omp parallel for collapse(2)
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > Target)
				RefineLoop(Grid[i][j][k],GridLevel[i][j][k],1+Target);
		}}}
//		#pragma omp parallel for collapse(2)
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > Target)
				SpatialInterpLoop1(Grid[i][j][k],GridLevel[i][j][k],1+Target);
		}}}
//		#pragma omp parallel for collapse(2)
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > Target)
				ChildNeighborLoop(Grid[i][j][k],GridLevel[i][j][k],1+Target);
		}}}
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > Target){
				MacroLoop(Grid[i][j][k],1+Target);
			}
		}}}
		}//end of parallel
	}
	#pragma omp parallel for collapse(3)
	for(i = 0; i<xDim; i++){
	for(j = 0; j<yDim; j++){
	for(k = 0; k<zDim; k++){
		Grid[i][j][k]->ComputeDelta(ObstX,ObstZ,ObstX2,ObstZ,ObstR);
	}}}
//	
//    double Residual = 1000.0;
	ofstream outputF;
	string ForceOut = casename;
	ForceOut+=".force";
	outputF.open (ForceOut.c_str());
//	
////	ofstream output_Vel1;
////	string VelOut = casename;
////	VelOut+="_uVelocityHist.dat";
////	output_Vel1.open (VelOut.c_str());
////	ofstream output_Vel2;
////	string VelOut2 = casename;
////	VelOut2+="_vVelocityHist.dat";
////	output_Vel2.open (VelOut2.c_str());
//
	
	
	ofstream outputVelhist1;
//	outputVelhist1.open((casename+".velhist1").c_str());

	vector<ofstream*> u_outputs1(1+maxlevel);
	vector<ofstream*> v_outputs1(1+maxlevel);
	vector<ofstream*> w_outputs1(1+maxlevel);
	vector<ofstream*> u_outputs2(1+maxlevel);
	vector<ofstream*> v_outputs2(1+maxlevel);
	vector<ofstream*> w_outputs2(1+maxlevel);
	vector<ofstream*> u_outputs3(1+maxlevel);
	vector<ofstream*> v_outputs3(1+maxlevel);
	vector<ofstream*> w_outputs3(1+maxlevel);
	string outputname;
	for(i = 0; i<maxlevel+1; i++){
//		outputname = casename;
//		outputname += "_lv";
//		stringstream lv;
//		lv<<i;	
//		outputname += lv.str();
		u_outputs1[i] = new ofstream;
		v_outputs1[i] = new ofstream;
		w_outputs1[i] = new ofstream;
//		u_outputs1[i]->open((outputname+".uVelocityHist1").c_str());
//		v_outputs1[i]->open((outputname+".vVelocityHist1").c_str());
//		w_outputs1[i]->open((outputname+".wVelocityHist1").c_str());
		u_outputs2[i] = new ofstream;
		v_outputs2[i] = new ofstream;
		w_outputs2[i] = new ofstream;
//		u_outputs2[i]->open((outputname+".uVelocityHist2").c_str());
//		v_outputs2[i]->open((outputname+".vVelocityHist2").c_str());
//		w_outputs2[i]->open((outputname+".wVelocityHist2").c_str());
		u_outputs3[i] = new ofstream;
		v_outputs3[i] = new ofstream;
		w_outputs3[i] = new ofstream;
//		u_outputs3[i]->open((outputname+".uVelocityHist3").c_str());
//		v_outputs3[i]->open((outputname+".vVelocityHist3").c_str());
//		w_outputs3[i]->open((outputname+".wVelocityHist3").c_str());
	}




	double sumError;
	double sumVel;
	int MinY, MinX, MinZ, MaxY, MaxX, MaxZ;

//////////////Start of time iteration///////////
	time_t start,end,current,previous,previous2,start2,current2;
	cout<<"Starting time iteration\n";
	time (&start);
	time (&previous);
	time (&previous2);

    for (int tStep = 0; tStep<tMax+1; tStep++){
//		//Remesh
		if(tStep == StartLR && maxlevel > 0){
////		if(tStep == 20000){// || tStep == 85 || tStep == 100){
////			time (&start2);
//////			VortCalc(Grid,vort);
			if(AMR == 1)
				Remap(GridLevel,Grid,Crit,maxlevel);

////			#pragma omp parallel for collapse(2)
////			for(i = 0; i<xDim; i++){
////			for(j = 0; j<yDim; j++){
////				CoarsenLoop(Grid[i][j],GridLevel[i][j]);
////			}}
////			

			if(maxlevel == 1){
			for(i = ObstX-15; i<1+ObstX+15+Spacing*2.0*ObstR; i++){
			for(j = 0; j<yDim; j++){
			for(k = ObstZ-10; k<1+ObstZ+10; k++){
				GridLevel[i][j][k] = 1;
			}}}
			}

			if(maxlevel == 2){
			for(i = ObstX-16; i<1+ObstX+16+Spacing*2.0*ObstR; i++){
			for(j = 0; j<yDim; j++){
			for(k = ObstZ-11; k<1+ObstZ+11; k++){
			//for(k = ObstZ-16; k<1+ObstZ+16; k++){
				GridLevel[i][j][k] = 1;
			}}}
			for(i = ObstX-15; i<1+ObstX+15+Spacing*2.0*ObstR; i++){
			for(j = 0; j<yDim; j++){
			for(k = ObstZ-10; k<1+ObstZ+10; k++){
			//for(k = ObstZ-15; k<1+ObstZ+15; k++){
				GridLevel[i][j][k] = 2;
			}}}
			}

			if(maxlevel == 3){
			for(i = ObstX-17; i<1+ObstX+17+Spacing*2.0*ObstR; i++){
			for(j = 0; j<yDim; j++){
			for(k = ObstZ-12; k<1+ObstZ+12; k++){
				GridLevel[i][j][k] = 1;
			}}}
			for(i = ObstX-16; i<1+ObstX+16+Spacing*2.0*ObstR; i++){
			for(j = 0; j<yDim; j++){
			for(k = ObstZ-11; k<1+ObstZ+11; k++){
				GridLevel[i][j][k] = 2;
			}}}
			for(i = ObstX-15; i<1+ObstX+15+Spacing*2.0*ObstR; i++){
			for(j = 0; j<yDim; j++){
			for(k = ObstZ-10; k<1+ObstZ+10; k++){
				GridLevel[i][j][k] = 3;
			}}}
			}

			for(Target = 0; Target<maxlevel; Target++){
				#pragma omp parallel
				{
				#pragma omp for collapse(3)
				for(i = 0; i<xDim; i++){
				for(j = 0; j<yDim; j++){
				for(k = 0; k<zDim; k++){
					if(GridLevel[i][j][k] > Target)
						RefineLoop(Grid[i][j][k],GridLevel[i][j][k],1+Target);
				}}}
				#pragma omp for collapse(3)
				for(i = 0; i<xDim; i++){
				for(j = 0; j<yDim; j++){
				for(k = 0; k<zDim; k++){
					if(GridLevel[i][j][k] > Target)
						SpatialInterpLoop1(Grid[i][j][k],GridLevel[i][j][k],1+Target);
				}}}
				#pragma omp for collapse(3)
				for(i = 0; i<xDim; i++){
				for(j = 0; j<yDim; j++){
				for(k = 0; k<zDim; k++){
					if(GridLevel[i][j][k] > Target)
						ChildNeighborLoop(Grid[i][j][k],GridLevel[i][j][k],1+Target);
				}}}
				#pragma omp for collapse(3)
				for(i = 0; i<xDim; i++){
				for(j = 0; j<yDim; j++){
				for(k = 0; k<zDim; k++){
					if(GridLevel[i][j][k] > Target){
						MacroLoop(Grid[i][j][k],1+Target);
					}
				}}}
				}
  			}
			for(Target = 0; Target<maxlevel; Target++){
				#pragma omp parallel for collapse(3)
				for(i = 0; i<xDim; i++){
				for(j = 0; j<yDim; j++){
				for(k = 0; k<zDim; k++){
					if(GridLevel[i][j][k] > Target){
						DeltaLoop(Grid[i][j][k],1+Target);
					}
				}}}
			}
		}
//
		sumError = 0.0;
	    sumVel = 0.0;
		FX = 0;//Force in -x
		FY = 0;//Force in -y
		FZ = 0;//Force in -z
		FX2 = 0;//Force in -x
		FY2 = 0;//Force in -y
		FZ2 = 0;//Force in -z
		MinY = yDim;
		MinX = xDim;
		MinZ = zDim;
		MaxY = 0;
		MaxX = 0;
		MaxZ = 0;


		if(InitPerturb == true){
		if(tStep == 1000 || tStep == 1200 || tStep == 1400 ||
			tStep == 2000 || tStep == 2200 || tStep == 2400){
			for(j = 0; j<yDim; j++){
				Grid[int(ObstX/2)][j][int(zDim/2)]->f[9] += 0.05;
			}
		}
		}

		//interpolate and apply boundary for fine mesh
		#pragma omp parallel //reduction(+:sumError,sumVel,FX,FY) 
		{
		if(maxlevel > 0){
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > 0){
				SpatialInterpLoop1(Grid[i][j][k],GridLevel[i][j][k],1);
				#pragma omp critical
				{
				if(i<MinX) MinX = i;
				if(j<MinY) MinY = j;
				if(k<MinZ) MinZ = k;
				if(i+1>MaxX) MaxX = i+1;
				if(j+1>MaxY) MaxY = j+1;
				if(k+1>MaxZ) MaxZ = k+1;
				}
			}
		}}}
		}
		
		int pos;

		//Part 1 for coarse mesh march
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
//			if(GridLevel[i][j][k] >= 0){
			//	MarchLoop1(Grid[i][j][k],0);
			Grid[i][j][k]->Stream1();
//			}
		}}}
		
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
//			if(GridLevel[i][j][k] >= 0){
			//MarchLoop2(Grid[i][j][k],0);
			Grid[i][j][k]->Stream2();
//			}
		}}}

//		#pragma omp for collapse(2)
//		for(j = 0; j<yDim; j++){
//		for(k = 0; k<zDim; k++){
//			Grid[0][j][k]->BoundaryMacro();
//			Grid[xDim-1][j][k]->BoundaryMacro();
//		}}


		//Fully Developed outlet
		#pragma omp for collapse(2)
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(Grid[xDim-1][j][k]->image == FullyDeveloped){
			for(i = 0; i<19; i++){
			Grid[xDim-1][j][k]->f[i] = Grid[xDim-2][j][k]->f[i];
			}
			}
			if(Grid[0][j][k]->image == FullyDeveloped){
			for(i = 0; i<19; i++){
			Grid[0][j][k]->f[i] = Grid[1][j][k]->f[i];
			}
			}
		}}	

		//Part 3 for coarse mesh march
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
//			if(GridLevel[i][j][k] >= 0){
//				MarchLoop3(Grid[i][j][k],0);
				Grid[i][j][k]->BoundaryMacro();
				Grid[i][j][k]->RegularizedBoundary();
				Grid[i][j][k]->Symmetry();
				Grid[i][j][k]->ComputeMacros();
				Grid[i][j][k]->ComputeFeq();
				Grid[i][j][k]->Collide();
				Grid[i][j][k]->BounceBack();
//			}
		}}}
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(CB == true) Grid[i][j][k]->CurvedB();
		}}}

		if(maxlevel > 0){
		//Interpolate and store boundary for t+dt in ftemp
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > 0)
				SpatialInterpLoop2(Grid[i][j][k],GridLevel[i][j][k],1);
		}}}
		//Temporal interpolation between f and ftemp, store in ftemp
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > 0)
				TemporalInterpLoop(Grid[i][j][k],GridLevel[i][j][k],1);	
		}}}
		}

		if(tStep >= StartAvg){
			//calculate average velocities on coarse mesh
			#pragma omp for collapse(3)
			for(i = 0; i<xDim; i++){
			for(j = 0; j<yDim; j++){
			for(k = 0; k<zDim; k++){
				Grid[i][j][k]->VelAverage(tStep);
			}}}
		}


		//Fully Developed outlet2
		#pragma omp for collapse(2)
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(Grid[xDim-1][j][k]->image == FullyDevelopedOut){
				Grid[xDim-2][j][k]->ComputeM();//compute post collision m
				for(i = 0; i<19; i++){
					if(i == 0) Grid[xDim-1][j][k]->m[i] = rhoIn;
					else{
					Grid[xDim-1][j][k]->m[i] = Grid[xDim-2][j][k]->m[i];
					}
				}
				Grid[xDim-1][j][k]->MToF();
			}
			if(Grid[0][j][k]->image == FullyDevelopedIn){
				Grid[1][j][k]->ComputeM();//compute post collision m
				for(i = 0; i<19; i++){
					if(i == 3) Grid[0][j][k]->m[i] = uMax;
					else if(i == 5) Grid[0][j][k]->m[i] = 0.0;
					else if(i == 7) Grid[0][j][k]->m[i] = 0.0;
					else{
					Grid[0][j][k]->m[i] = Grid[1][j][k]->m[i];
					}
				}
				Grid[0][j][k]->MToF();
			}
		}}	

		}//end omp parallel

		//write u v to file
		if(tStep >= StartRec){
			if(tStep == StartRec){
				outputVelhist1.open((casename+".velhist1").c_str());
				for(i = 0; i<maxlevel+1; i++){
					outputname = casename;
					outputname += "_lv";
					stringstream lv;
					lv<<i;	
					outputname += lv.str();
					u_outputs1[i]->open((outputname+".uVelocityHist1").c_str());
					v_outputs1[i]->open((outputname+".vVelocityHist1").c_str());
					w_outputs1[i]->open((outputname+".wVelocityHist1").c_str());
					u_outputs2[i]->open((outputname+".uVelocityHist2").c_str());
					v_outputs2[i]->open((outputname+".vVelocityHist2").c_str());
					w_outputs2[i]->open((outputname+".wVelocityHist2").c_str());
					u_outputs3[i]->open((outputname+".uVelocityHist3").c_str());
					v_outputs3[i]->open((outputname+".vVelocityHist3").c_str());
					w_outputs3[i]->open((outputname+".wVelocityHist3").c_str());
				}
			}
			*u_outputs1[0]<<tStep<<",";
			*v_outputs1[0]<<tStep<<",";
			*w_outputs1[0]<<tStep<<",";
			*u_outputs2[0]<<tStep<<",";
			*v_outputs2[0]<<tStep<<",";
			*w_outputs2[0]<<tStep<<",";
			*u_outputs3[0]<<tStep<<",";
			*v_outputs3[0]<<tStep<<",";
			*w_outputs3[0]<<tStep<<",";
			for(i = 0;i<xDim;i++){
				(*(u_outputs1[0]))<<Grid[i][int((yDim-1)/2)][int((zDim-1)/2)]->u0<<",";
				(*(v_outputs1[0]))<<Grid[i][int((yDim-1)/2)][int((zDim-1)/2)]->u1<<",";
				(*(w_outputs1[0]))<<Grid[i][int((yDim-1)/2)][int((zDim-1)/2)]->u2<<",";
			}
			for(j = 0;j<yDim;j++){
				(*(u_outputs2[0]))<<Grid[int((xDim-1)/2)][j][int((zDim-1)/2)]->u0<<",";
				(*(v_outputs2[0]))<<Grid[int((xDim-1)/2)][j][int((zDim-1)/2)]->u1<<",";
				(*(w_outputs2[0]))<<Grid[int((xDim-1)/2)][j][int((zDim-1)/2)]->u2<<",";
			}
			for(k = 0;k<zDim;k++){
				(*(u_outputs3[0]))<<Grid[int((xDim-1)/2)][int((yDim-1)/2)][k]->u0<<",";
				(*(v_outputs3[0]))<<Grid[int((xDim-1)/2)][int((yDim-1)/2)][k]->u1<<",";
				(*(w_outputs3[0]))<<Grid[int((xDim-1)/2)][int((yDim-1)/2)][k]->u2<<",";
			}
			*u_outputs1[0]<<endl;
			*v_outputs1[0]<<endl;
			*w_outputs1[0]<<endl;
			*u_outputs2[0]<<endl;
			*v_outputs2[0]<<endl;
			*w_outputs2[0]<<endl;
			*u_outputs3[0]<<endl;
			*v_outputs3[0]<<endl;
			*w_outputs3[0]<<endl;
		}
		
		if(MaxX > MinX && MaxY > MinY)
			LRLBM(Grid,GridLevel,1,MinX,MinY,MinZ,MaxX,MaxY,MaxZ,tStep,u_outputs1,v_outputs1,u_outputs2,v_outputs2);
		
		#pragma omp parallel for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > 0){
				SpatialAverageLoop(Grid[i][j][k],GridLevel[i][j][k],1);
			}
		}}}

		if(tStep >= StartF){// && tStep%1 == 0){
			#pragma omp master
			{
			if(tStep < StartLR || maxlevel == 0){
			for(i = 2; i<xDim-2; i++){
			for(j = 0; j<yDim; j++){
			for(k = 2; k<zDim-2; k++){
				BBForceLoop(Grid[i][j][k],0,FX,FY,FZ,ForceB1);
				BBForceLoop(Grid[i][j][k],0,FX2,FY2,FZ2,ForceB2);
			}}}
			}
			for(i = 0; i<FXbuff.size(); i++){
				outputF<<tStep<<", "<<FXbuff[i]/(0.5*uMax*uMax*ObstR*2.0*pow(2.0,2.0*maxlevel)*(yDim))
							  <<", "<<FYbuff[i]/(0.5*uMax*uMax*ObstR*2.0*pow(2.0,2.0*maxlevel)*(yDim))
							  <<", "<<FZbuff[i]/(0.5*uMax*uMax*ObstR*2.0*pow(2.0,2.0*maxlevel)*(yDim))
							  <<", "<<FX2buff[i]/(0.5*uMax*uMax*ObstR*2.0*pow(2.0,2.0*maxlevel)*(yDim))
							  <<", "<<FY2buff[i]/(0.5*uMax*uMax*ObstR*2.0*pow(2.0,2.0*maxlevel)*(yDim))
							  <<", "<<FZ2buff[i]/(0.5*uMax*uMax*ObstR*2.0*pow(2.0,2.0*maxlevel)*(yDim))<<"\n";
			}
			if(FXbuff.size() == 0){
				outputF<<tStep<<", "<<FX/(0.5*uMax*uMax*ObstR*2.0*(yDim))
							  <<", "<<FY/(0.5*uMax*uMax*ObstR*2.0*(yDim))
							  <<", "<<FZ/(0.5*uMax*uMax*ObstR*2.0*(yDim))
							  <<", "<<FX2/(0.5*uMax*uMax*ObstR*2.0*(yDim))
							  <<", "<<FY2/(0.5*uMax*uMax*ObstR*2.0*(yDim))
							  <<", "<<FZ2/(0.5*uMax*uMax*ObstR*2.0*(yDim))<<"\n";
			}
			FXbuff.resize(0,0);
			FYbuff.resize(0,0);
			FZbuff.resize(0,0);
			FX2buff.resize(0,0);
			FY2buff.resize(0,0);
			FZ2buff.resize(0,0);
			}
		}

		//added maco loop for data writing
		if(tStep == tMax){
		#pragma omp parallel //reduction(+:FX,FY)
		{
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > 0)
				SpatialInterpLoop1(Grid[i][j][k],GridLevel[i][j][k],1);
		}}}
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(GridLevel[i][j][k] > 0){
				SpatialAverageLoop(Grid[i][j][k],GridLevel[i][j][k],1);
				MacroLoop(Grid[i][j][k],1);
			}
		}}}
		#pragma omp for collapse(3)
		for(i = 0; i<xDim; i++){
		for(j = 0; j<yDim; j++){
		for(k = 0; k<zDim; k++){
			if(Grid[i][j][k]->image >= BB && Grid[i][j][k]->image <= CB2){
				Grid[i][j][k]->u0 = 0;
				Grid[i][j][k]->u1 = 0;
				Grid[i][j][k]->u2 = 0;
				Grid[i][j][k]->rho = 1.0;
			}
			else{
				Grid[i][j][k]->ComputeMacros();
			}
		}}}
		}//end of parallel
		}

//		}//end of parallel
//
		
		//velocity history monitor
		if(tStep >= StartAvg){
			outputVelhist1<<Grid[mp1x][mp1y][mp1z]->usqr<<", "
						  <<Grid[mp2x][mp2y][mp2z]->usqr<<", "
						  <<Grid[mp3x][mp3y][mp3z]->usqr<<", "
						  <<Grid[mp4x][mp4y][mp4z]->usqr<<endl;
		}	
		
		
		if((tStep != 0) && (tStep%250 == 0)){
			time (&current);
	        cout<<"Finished timestep: "<<tStep<<endl;//<<", Residual = "<<Residual<<"\n";
//			cout<<"Cd = "<<FX/(0.5*uMax*uMax*ObstR*2.0)<<", Fx:"<<FX<<"\n";
//			cout<<"Cl = "<<-FY/(0.5*uMax*uMax*ObstR*2.0)<<", Fy:"<<FY<<"\n";
			cout<<"Current wall clock time: "<<difftime(current,start)<<"\n";
			cout<<"Estimated time to completion: "
					<<(tMax-tStep)*difftime(current,previous2)/(100*3600)<<"hrs\n";
			time (&previous2);
		}
		if(tStep%100 == 0){
			ofstream outcheck;
			time (&current);
			string outcheckname = "outcheck_"+casename+".dat";
			outcheck.open(outcheckname.c_str());
	        	outcheck<<"Finished timestep: "<<tStep<<"\n";
//			outcheck<<"Cd = "<<FX/(0.5*uMax*uMax*ObstR*2.0)<<", Fx:"<<FX<<"\n";
//			outcheck<<"Cl = "<<-FY/(0.5*uMax*uMax*ObstR*2.0)<<", Fy:"<<FY<<"\n";
			outcheck<<"Current wall clock time: "<<difftime(current,start)<<"\n";
			outcheck<<"Estimated time to completion: "
					<<(tMax-tStep)*difftime(current,previous)/(100*3600)<<"hrs\n";
			time (&previous);
			outcheck.close();
		}	
////	    if(tStep > 100 && Residual<Tolerance){
////	        ConvFlag = 1;
////	        cout<<"Converged solution obtained at timestep: "<<tStep<<"\n";
////	        cout<<"Residual = "<<Residual<<"\n";
////	        tStep = tMax;//break;
////	    }
	    if(ErrorFlag == 1){
	        cout<<"Error has occured. Exiting at timestep: "<<tStep<<"\n";
	        tStep = tMax;//break;
	    }
//	    Residual = sumError/sumVel;
//
//
//		//Write output
		if(tStep == tMax || tStep == tMax/2 || tStep == tMax/4 || tStep == 1000){

			int num = tStep;
			string s;
			stringstream out;
			out<<num;
			s = out.str();
			string filename = casename+"_"+s+".dat";
//			VortCalc(Grid,vort);
			vort_iCalc(Grid,vort_i);
			vort_jCalc(Grid,vort_j);
			vort_kCalc(Grid,vort_k);
			SmagCalc(Grid,smag);
		    ofstream output;
		    output.open (filename.c_str() );
//		    output<<"VARIABLES = \"X\",\"Y\",\"Z\",\"u\",\"v\",\"w\",\"rho\",\"u_avg\",\"v_avg\",\"w_avg\",\"k_vort\"\n";
		    output<<"VARIABLES = \"X\",\"Y\",\"Z\",\"u\",\"v\",\"w\",\"rho\",\"u_avg\",\"v_avg\",\"w_avg\",\"i_vort\",\"j_vort\",\"k_vort\",\"Smag\"\n";
		    output<<"ZONE F=POINT, I = "<<yDim<<", J = "<<xDim<<", K = "<<zDim<<"\n";
			//Write coarse mesh data
		    for(k = 0; k<zDim; k++){
		    for(i = 0; i<xDim; i++){
		    for(j = 0; j<yDim; j++){	
		        output<<i<<", "<<j<<", "<<k<<", "<<
					Grid[i][j][k]->u0<<", "<<Grid[i][j][k]->u1<<", "<<Grid[i][j][k]->u2;
		        output<<", "<<Grid[i][j][k]->rho;
		        output<<", "<<Grid[i][j][k]->u0avg<<", "<<Grid[i][j][k]->u1avg
						<<", "<<Grid[i][j][k]->u2avg<<", "<<vort_i[i][j][k]<<", "
						<<vort_j[i][j][k]<<", "<<vort_k[i][j][k]<<", "
						<<Grid[i][j][k]->Smag;
				output<<"\n";
//		        output<<i<<", "<<j<<", "<<", "<<k<<", ";
//				for(int a = 8; a<19; a++)
//					output<<Grid[i][j][k]->delta[a]<<", ";
//				output<<"\n";
		    }}}

			if(SaveState == true && tStep == tMax){
		    	ofstream outputstate;
		    	outputstate.open((casename+".state").c_str());
			    for(i = 0; i<xDim; i++){
			    for(j = 0; j<yDim; j++){
			    for(k = 0; k<zDim; k++){
					outputstate<<i<<","<<j<<","<<k<<",";
			    	for(l = 0; l<19; l++){
						outputstate<<Grid[i][j][k]->f[l]<<",";
					}
					outputstate<<endl;
				}}}
				outputstate.close();
			}

			for(n = 1; n<maxlevel+2; n++){
			//Count nodes and elements in refinement region
			int nodes = 0;
			int elements = 0;
		    for(k = 0; k<zDim; k++){
			for(i = 0; i<xDim; i++){
			for(j = 0; j<yDim; j++){
				if(GridLevel[i][j][k] >= n)
					CountElements((Grid[i][j][k]),n,nodes,elements);
			}}}
			if(nodes == 0) break;
//			output<<"VARIABLES = \"X\",\"Y\",\"u\",\"v\",\"rho\"\n";
//		    output<<"VARIABLES = \"X\",\"Y\",\"Z\",\"u\",\"v\",\"w\",\"rho\",\"u_avg\",\"v_avg\",\"w_avg\",\"k_vort\"\n";
		    output<<"VARIABLES = \"X\",\"Y\",\"Z\",\"u\",\"v\",\"w\",\"rho\",\"u_avg\",\"v_avg\",\"w_avg\",\"i_vort\",\"j_vort\",\"k_vort\",\"Smag\"\n";
			output<<"ZONE N="<<nodes<<", E="<<elements<<
					", DATAPACKING=POINT, ZONETYPE=FEBRICK\n";
			int numb = 1;//initialize data line number count
			for(i = 0; i<xDim; i++){
			for(j = 0; j<yDim; j++){
		    for(k = 0; k<zDim; k++){
				if(GridLevel[i][j][k] >= n)
					WriteData(output,(Grid[i][j][k]),n,numb);
			}}}
			for(i = 0; i<xDim; i++){
			for(j = 0; j<yDim; j++){
		    for(k = 0; k<zDim; k++){
				if(GridLevel[i][j][k] >= n)
					WriteElements(output,(Grid[i][j][k]),n);
			}}}
			}
			output.close();
		}//end of output 
////		}//end of omp master
////		}//end of omp parallel
    } 
////////////////End of time iteration///////////

//	for(i = 1; i<xDim-1; i++){
//	for(k = 1; k<zDim-1; k++){
//	for(l = 1; l<19; l++){
//		cout<<Grid[i][0][k]->nb[l]->ycoord<<",";
////		cout<<Grid[0][0][2]->nb[l]->xcoord<<","<<Grid[0][0][2]->nb[l]->ycoord<<","<<Grid[0][0][2]->nb[l]->zcoord<<endl;
//	}
//	cout<<endl;
//	}}
//
///	time (&end);
//	cout<<"Finished run in: "<<difftime(end,start)<<"s\n";
	outputF.close();
////	for(i = 0; i<maxlevel+1; i++){
////		u_outputs1[i]->close();
////		v_outputs1[i]->close();
////	}
//	
	
//	cout<<Grid[5][5][1]->u1<<endl;
//
//	cout<<endl;
//
//	for(i = 0; i< 19; i++){
//		cout<<i<<", "<<Grid[5][5][1]->m[i]<<endl;
//	}

	outputVelhist1.close();

	for(i = 0; i<maxlevel+1; i++){
		u_outputs1[i]->close();
		v_outputs1[i]->close();
		w_outputs1[i]->close();
		u_outputs2[i]->close();
		v_outputs2[i]->close();
		w_outputs2[i]->close();
		u_outputs3[i]->close();
		v_outputs3[i]->close();
		w_outputs3[i]->close();
	}

	for(i = 0; i<xDim; i++){
	for(j = 0; j<yDim; j++){
	for(k = 0; k<zDim; k++){
		CoarsenLoop(Grid[i][j][k],0);
	}}}
	for(i = 0; i<xDim; i++){
	for(j = 0; j<yDim; j++){
	for(k = 0; k<zDim; k++){
		Grid[i][j][k]->DeleteNode();
	}
	delete [] Grid[i][j];
	delete [] GridLevel[i][j];
	delete [] vort_i[i][j];
	delete [] vort_j[i][j];
	delete [] vort_k[i][j];
	}
	delete [] Grid[i];
	delete [] GridLevel[i];
	delete [] vort_i[i];
	delete [] vort_j[i];
	delete [] vort_k[i];
	}
	delete [] Grid;
	delete [] GridLevel;
	delete [] vort_i;
	delete [] vort_j;
	delete [] vort_k;

	return 0;
}

