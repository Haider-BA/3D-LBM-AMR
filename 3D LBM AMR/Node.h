#ifndef NODE_H
#define NODE_H

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
#include <iomanip>

//Class Definitions
class Node
{
public: 
    double u0,u1,u2,usqr,rho,xcoord,ycoord,zcoord,nu0;
//    double u0_prev,u1_prev,u_error,u_sum;
//	double m1eq,m2eq,m4eq,m6eq,m7eq,m8eq;//moments (starting from 0, to 8)
//	double m1,m2,m4,m6,m7,m8;//moments (starting from 0, to 8)
//	double Q11,Q12,Q21,Q22,Q,Smag,OmegaStar,Cs;
	double Smag,OmegaStar,Cs;
//	double PI11, PI22, PI33, PI12, PI13, PI23;
//	double PI11,PI22,PI12;
	double u0avg, u1avg, u2avg;
    int image,level,outnumb,edge;
	bool feqflag;//1 if MeqToFeq was used
	Node * Parent;
    std::vector<double> f,feq,ftemp,delta,m,meq;
	std::vector<Node *> nb;
    Node * child[8];
	Node();
	Node(Node& parent);
    void Initialize(double u0_ini,double u1_ini,double u2_ini,double rho_ini,int x,int y,int z);
    void SetImage(int Image_val);
	void ComputeM();//compute M from f and meq
    void ComputeFeq();
    void Collide();
    void ComputeMacros();
	void VelAverage(int tStep);//compute time averaged velocity
    void Symmetry();
	void Periodic(Node *);
    void BounceBack();
    void CurvedB();
	void BBForce(double& ForceX,double& ForceY,double& ForceZ,int BBimage);
    void BoundaryMacro();
    void InletOutlet();
    void RegularizedBoundary();
    void DeleteNode();
    void Refine();
	void ChildNeighbor();
	void Stream1();
	void Stream2();
	void SpatialInterp1();
	void SpatialInterp2();
	void SpatialAverage();
	void TemporalInterp();
	void TemporalInterp2();
	void MeqToFeq();
	void MToF();
	void ComputeDelta(double obX,double obY,double obX2,double obY2,double R);
//	void Outlet2();

};


#endif