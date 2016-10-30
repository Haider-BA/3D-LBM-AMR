#ifndef INPUTVARIABLES_H
#define INPUTVARIABLES_H

#include <vector>


//Simulation Parameters
extern std::string casename;

extern double Spacing;//0.0;
extern double ObstR;//2.5;//5;//5.0;//10.0;//.0;
extern double ObstX;//50.5;//25;//70;//13.0*2.0*ObstR;
extern double ObstX2;//ObstX+Spacing*2.0*ObstR;//58;//70;//13.0*2.0*ObstR;
extern double ObstY;//2;//2;//50.5;//51.0*2.0*ObstR/2;//1022/2;
extern double ObstZ;//100.5;//100;//35.5;//70.5;//50.5;//51.0*2.0*ObstR/2;//1022/2;
extern int xDim;//ObstX+15*2*ObstR+1+Spacing*2.0*ObstR;
extern int yDim;//1;//1;//4;//61;//11;//102;//281; //(8.2*ObstR)+2;
extern int zDim;//201;//71;//141;//102;//281; //(8.2*ObstR)+2;

extern int maxlevel;//2;

extern double Re;//100;//21400;
extern double uMax;//0.08;//1;
extern double rhoIn;//1.0;
extern double CharLength;//ObstR*2.0;
extern double dPdx;//0;//5.56e-7;
extern double Tolerance;//1.0e-6;

extern bool D2Q9i;//true;
extern bool MRT;//1;//1;//false;//true;
extern bool LES;//0;//true;//true;
extern bool AMR;//0;//true;//true;
extern bool CB;//1;//true;//true; //Curved boundary treatment
extern int  ForceB1;//CB1;//boundary for force calculations
extern int  ForceB2;//CB2;//boundary for force calculations

extern bool SaveState;//1;
extern bool Periodics;//0;
extern bool xPeriodics;//0;
extern bool yPeriodics;//1;
extern bool zPeriodics;//0;

extern bool InitPerturb;//0;
extern bool InitCond;//1;//false;//true;
//string initfile ("CC6_init1.state");

extern int tMax;//20000;//10000;//10000;

extern int StartAvg;// =1000;
extern int StartRec;// = 200000;//30;// 30000;//start recording velocity history
extern int StartF;//0;//5000;//start recording force history
extern int StartLR;//1;//start LR

extern int nCPU;//2;

extern double s1 ;//1.0;//1.19;
extern double s2 ;//1.0;//1.4;
extern double s4 ;//1.0;//1.2;
extern double s10;//1.0;//1.4;
extern double s16;//1.0;//1.98;
extern bool scaleS;//0;//1; //scale s's for LR



extern int mp1x;//xDim/4;
extern int mp1y;//yDim/4;
extern int mp1z;//zDim/4;
extern int mp2x;//xDim-1-xDim/4;
extern int mp2y;//yDim/4;
extern int mp2z;//zDim/4;
extern int mp3x;//xDim/4;
extern int mp3y;//yDim-1-yDim/4;
extern int mp3z;//zDim-1-zDim/4;
extern int mp4x;//xDim-1-xDim/4;
extern int mp4y;//yDim-1-yDim/4;
extern int mp4z;//zDim-1-zDim/4;

extern std::vector<double> S1;// (0);
extern std::vector<double> S2;//(0);
extern std::vector<double> S4;//(0);
extern std::vector<double> S10;//(0);
extern std::vector<double> S16;//(0);





//Compute Omega
extern double nu;// = uMax*(CharLength) / Re;
extern double omega1;// = 1.0/(3.0*nu+0.5);
extern std::vector<double> omega;




extern std::vector<double> FXbuff;//pow(2,maxlevel));
extern std::vector<double> FYbuff;//pow(2,maxlevel));
extern std::vector<double> FZbuff;//pow(2,maxlevel));
extern std::vector<double> FX2buff;//pow(2,maxlevel));
extern std::vector<double> FY2buff;//pow(2,maxlevel));
extern std::vector<double> FZ2buff;//pow(2,maxlevel));



//relaxation rates for the 19 moments
extern double S[19];

//Global vars
extern int ErrorFlag;



////LBM CONSTANTS

//D2Q9 extern
extern double t [19];//{4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
extern int v [19][3];//{0,0,1,0,0,1,-1,0,0,-1,1,1,-1,1,-1,-1,1,-1}; 
extern int opposite [19];//{0,3,4,1,2,7,8,5,6};
extern int xsym [19];//{0,1,4,3,2,8,7,6,5};
extern int ysym [19];//{0,3,2,1,4,6,5,8,7};

extern int halfstream [9];//{1,2,5,6};

extern double M[19][19];
extern double M_inv [19][19];
extern double psi_mag2 [19]; //square of the magnitude of vectors psi, constituting M


// global functions for image and basic interpolations

double PoisProf (double x);

bool CylinderImage(double x, double y);

bool CylinderImage2(int x, int y);

int Geometry(double x, double y, double z);

double HermitSpline(double f1, double f2, double f3, double f4);


#endif
