#include "LBMConstants.h"
#include "InputVariables.h"
#include "Node.h"

using namespace std;


Node::Node (){
	level = 0;
	for(int i = 0; i<8; i++)
        child[i] = NULL;
}

Node::Node (Node& parent){
	Parent = &parent;
	level = Parent->level+1;
	nb.resize(19);
	nb[0] = NULL;
	f.resize(19,0);
    feq.resize(19,0);
    ftemp.resize(19,0);
	if(MRT == true){
	    m.resize(19,0);
	    meq.resize(19,0);
		feqflag = 1;//1 if MeqToFeq was used
	}
	image = fluid;
	u0avg = 0;
	u1avg = 0;
	u2avg = 0;
	edge = 1;//newly refined nodes need to be interpolated
	for(int i = 0; i<8; i++)
        child[i] = NULL;
}
void Node::Refine()
{
    double offset = 0.25*(pow(2.0,-double(level)));//0.25 of current grid size

	child[0] = new Node(*this);//(level+1);//(2);
	child[0]->xcoord = xcoord-offset;
	child[0]->ycoord = ycoord-offset;
	child[0]->zcoord = zcoord-offset;
	child[0]->image = Geometry(child[0]->xcoord,child[0]->ycoord,child[0]->zcoord);
	
	child[1] = new Node(*this);//(3);
    child[1]->xcoord = xcoord+offset;
	child[1]->ycoord = ycoord-offset;
	child[1]->zcoord = zcoord-offset;
	child[1]->image = Geometry(child[1]->xcoord,child[1]->ycoord,child[1]->zcoord);

	child[2] = new Node(*this);//(3);
    child[2]->xcoord = xcoord+offset;
	child[2]->ycoord = ycoord+offset;
	child[2]->zcoord = zcoord-offset;
	child[2]->image = Geometry(child[2]->xcoord,child[2]->ycoord,child[2]->zcoord);
	
	child[3] = new Node(*this);//(5);
	child[3]->xcoord = xcoord-offset;
	child[3]->ycoord = ycoord+offset;
	child[3]->zcoord = zcoord-offset;
	child[3]->image = Geometry(child[3]->xcoord,child[3]->ycoord,child[3]->zcoord);

	child[4] = new Node(*this);//(level+1);//(2);
	child[4]->xcoord = xcoord-offset;
	child[4]->ycoord = ycoord-offset;
	child[4]->zcoord = zcoord+offset;
	child[4]->image = Geometry(child[4]->xcoord,child[4]->ycoord,child[4]->zcoord);
	
	child[5] = new Node(*this);//(3);
    child[5]->xcoord = xcoord+offset;
	child[5]->ycoord = ycoord-offset;
	child[5]->zcoord = zcoord+offset;
	child[5]->image = Geometry(child[5]->xcoord,child[5]->ycoord,child[5]->zcoord);

	child[6] = new Node(*this);//(3);
    child[6]->xcoord = xcoord+offset;
	child[6]->ycoord = ycoord+offset;
	child[6]->zcoord = zcoord+offset;
	child[6]->image = Geometry(child[6]->xcoord,child[6]->ycoord,child[6]->zcoord);
	
	child[7] = new Node(*this);//(5);
	child[7]->xcoord = xcoord-offset;
	child[7]->ycoord = ycoord+offset;
	child[7]->zcoord = zcoord+offset;
	child[7]->image = Geometry(child[7]->xcoord,child[7]->ycoord,child[7]->zcoord);

	for(int i = 0; i<8; i++){
		if(child[i]->image == CB1 || child[i]->image == CB2) child[i]->delta.resize(19,0);
	}
}

void Node::ChildNeighbor()
{

	child[0]->nb[1] = child[1];
	child[0]->nb[2] = child[3];
	child[0]->nb[5] = child[2];
	child[0]->nb[9] = child[4];
	child[0]->nb[10] = child[5];
	child[0]->nb[11] = child[7];
	child[0]->edge = 0;
	if(nb[3]->child[0] != NULL){
		child[0]->nb[3] = nb[3]->child[1];
		child[0]->nb[6] = nb[3]->child[2];
		child[0]->nb[12] = nb[3]->child[5];}
	else{
		child[0]->nb[3] = NULL;
		child[0]->nb[6] = NULL;
		child[0]->nb[12] = NULL;
		child[0]->edge = 1;}

	if(nb[4]->child[0] != NULL){
		child[0]->nb[4] = nb[4]->child[3];
		child[0]->nb[8] = nb[4]->child[2];
		child[0]->nb[13] = nb[4]->child[7];}
	else{
		child[0]->nb[4] = NULL;
		child[0]->nb[8] = NULL;
		child[0]->nb[13] = NULL;
		child[0]->edge = 1;}

	if(nb[7]->child[0] != NULL)
		child[0]->nb[7] = nb[7]->child[2];
	else{
		child[0]->nb[7] = NULL;
		child[0]->edge = 1;}

	if(nb[14]->child[0] != NULL){
		child[0]->nb[14] = nb[14]->child[4];
		child[0]->nb[15] = nb[14]->child[5];
		child[0]->nb[16] = nb[14]->child[7];}
	else{
		child[0]->nb[14] = NULL;
		child[0]->nb[15] = NULL;
		child[0]->nb[16] = NULL;
		child[0]->edge = 1;}

	if(nb[17]->child[0] != NULL)
		child[0]->nb[17] = nb[17]->child[5];
	else{
		child[0]->nb[17] = NULL;
		child[0]->edge = 1;}
	if(nb[18]->child[0] != NULL)
		child[0]->nb[18] = nb[18]->child[7];
	else{
		child[0]->nb[18] = NULL;
		child[0]->edge = 1;}

	child[1]->nb[2] = child[2];
	child[1]->nb[3] = child[0];
	child[1]->nb[6] = child[3];
	child[1]->nb[9] = child[5];
	child[1]->nb[11] = child[6];
	child[1]->nb[12] = child[4];
	child[1]->edge = 0;
	if(nb[1]->child[0] != NULL){
		child[1]->nb[1] = nb[1]->child[0];
		child[1]->nb[5] = nb[1]->child[3];
		child[1]->nb[10] = nb[1]->child[4];}
	else{
		child[1]->nb[1] = NULL;
		child[1]->nb[5] = NULL;
		child[1]->nb[10] = NULL;
		child[1]->edge = 1;}

	if(nb[4]->child[0] != NULL){
		child[1]->nb[4] = nb[4]->child[2];
		child[1]->nb[7] = nb[4]->child[3];
		child[1]->nb[13] = nb[4]->child[6];}
	else{
		child[1]->nb[4] = NULL;
		child[1]->nb[7] = NULL;
		child[1]->nb[13] = NULL;
		child[1]->edge = 1;}

	if(nb[8]->child[0] != NULL)
		child[1]->nb[8] = nb[8]->child[3];
	else{
		child[1]->nb[8] = NULL;
		child[1]->edge = 1;}

	if(nb[14]->child[0] != NULL){
		child[1]->nb[14] = nb[14]->child[5];
		child[1]->nb[16] = nb[14]->child[6];
		child[1]->nb[17] = nb[14]->child[4];}
	else{
		child[1]->nb[14] = NULL;
		child[1]->nb[16] = NULL;
		child[1]->nb[17] = NULL;
		child[1]->edge = 1;}

	if(nb[15]->child[0] != NULL)
		child[1]->nb[15] = nb[15]->child[4];
	else{
		child[1]->nb[15] = NULL;
		child[1]->edge = 1;}
	if(nb[18]->child[0] != NULL)
		child[1]->nb[18] = nb[18]->child[6];
	else{
		child[1]->nb[18] = NULL;
		child[1]->edge = 1;}

	child[2]->nb[3] = child[3];
	child[2]->nb[4] = child[1];
	child[2]->nb[7] = child[0];
	child[2]->nb[9] = child[6];
	child[2]->nb[12] = child[7];
	child[2]->nb[13] = child[5];
	child[2]->edge = 0;
	if(nb[1]->child[0] != NULL){
		child[2]->nb[1] = nb[1]->child[3];
		child[2]->nb[8] = nb[1]->child[0];
		child[2]->nb[10] = nb[1]->child[7];}
	else{
		child[2]->nb[1] = NULL;
		child[2]->nb[8] = NULL;
		child[2]->nb[10] = NULL;
		child[2]->edge = 1;}
	if(nb[2]->child[0] != NULL){
		child[2]->nb[2] = nb[2]->child[1];
		child[2]->nb[6] = nb[2]->child[0];
		child[2]->nb[11] = nb[2]->child[5];}
	else{
		child[2]->nb[2] = NULL;
		child[2]->nb[6] = NULL;
		child[2]->nb[11] = NULL;
		child[2]->edge = 1;}
	if(nb[5]->child[0] != NULL)
		child[2]->nb[5] = nb[5]->child[0];
	else{
		child[2]->nb[5] = NULL;
		child[2]->edge = 1;}

	if(nb[15]->child[0] != NULL)
		child[2]->nb[15] = nb[15]->child[7];
	else{
		child[2]->nb[15] = NULL;
		child[2]->edge = 1;}
	if(nb[16]->child[0] != NULL)
		child[2]->nb[16] = nb[16]->child[5];
	else{
		child[2]->nb[16] = NULL;
		child[2]->edge = 1;}

	if(nb[14]->child[0] != NULL){
		child[2]->nb[14] = nb[14]->child[6];
		child[2]->nb[17] = nb[14]->child[7];
		child[2]->nb[18] = nb[14]->child[5];}
	else{
		child[2]->nb[14] = NULL;
		child[2]->nb[17] = NULL;
		child[2]->nb[18] = NULL;
		child[2]->edge = 1;}

	child[3]->nb[1] = child[2];
	child[3]->nb[4] = child[0];
	child[3]->nb[8] = child[1];
	child[3]->nb[9] = child[7];
	child[3]->nb[10] = child[6];
	child[3]->nb[13] = child[4];
	child[3]->edge = 0;
	if(nb[2]->child[0] != NULL){
		child[3]->nb[2] = nb[2]->child[0];
		child[3]->nb[5] = nb[2]->child[1];
		child[3]->nb[11] = nb[2]->child[4];}
	else{
		child[3]->nb[2] = NULL;
		child[3]->nb[5] = NULL;
		child[3]->nb[11] = NULL;
		child[3]->edge = 1;}
	if(nb[3]->child[0] != NULL){
		child[3]->nb[3] = nb[3]->child[2];
		child[3]->nb[7] = nb[3]->child[1];
		child[3]->nb[12] = nb[3]->child[6];}
	else{
		child[3]->nb[3] = NULL;
		child[3]->nb[7] = NULL;
		child[3]->nb[12] = NULL;
		child[3]->edge = 1;}
	if(nb[6]->child[0] != NULL)
		child[3]->nb[6] = nb[6]->child[1];
	else{
		child[3]->nb[6] = NULL;
		child[3]->edge = 1;}

	if(nb[16]->child[0] != NULL)
		child[3]->nb[16] = nb[16]->child[4];
	else{
		child[3]->nb[16] = NULL;
		child[3]->edge = 1;}
	if(nb[17]->child[0] != NULL)
		child[3]->nb[17] = nb[17]->child[6];
	else{
		child[3]->nb[17] = NULL;
		child[3]->edge = 1;}

	if(nb[14]->child[0] != NULL){
		child[3]->nb[14] = nb[14]->child[7];
		child[3]->nb[15] = nb[14]->child[6];
		child[3]->nb[18] = nb[14]->child[4];}
	else{
		child[3]->nb[14] = NULL;
		child[3]->nb[15] = NULL;
		child[3]->nb[18] = NULL;
		child[3]->edge = 1;}

	child[4]->nb[1] = child[5];
	child[4]->nb[2] = child[7];
	child[4]->nb[5] = child[6];
	child[4]->nb[14] = child[0];
	child[4]->nb[15] = child[1];
	child[4]->nb[16] = child[3];
	child[4]->edge = 0;
	if(nb[3]->child[0] != NULL){
		child[4]->nb[3] = nb[3]->child[5];
		child[4]->nb[6] = nb[3]->child[6];
		child[4]->nb[17] = nb[3]->child[1];}
	else{
		child[4]->nb[3] = NULL;
		child[4]->nb[6] = NULL;
		child[4]->nb[17] = NULL;
		child[4]->edge = 1;}

	if(nb[4]->child[0] != NULL){
		child[4]->nb[4] = nb[4]->child[7];
		child[4]->nb[8] = nb[4]->child[6];
		child[4]->nb[18] = nb[4]->child[3];}
	else{
		child[4]->nb[4] = NULL;
		child[4]->nb[8] = NULL;
		child[4]->nb[18] = NULL;
		child[4]->edge = 1;}

	if(nb[7]->child[0] != NULL)
		child[4]->nb[7] = nb[7]->child[6];
	else{
		child[4]->nb[7] = NULL;
		child[4]->edge = 1;}

	if(nb[9]->child[0] != NULL){
		child[4]->nb[9] = nb[9]->child[0];
		child[4]->nb[10] = nb[9]->child[1];
		child[4]->nb[11] = nb[9]->child[3];}
	else{
		child[4]->nb[9] = NULL;
		child[4]->nb[10] = NULL;
		child[4]->nb[11] = NULL;
		child[4]->edge = 1;}

	if(nb[12]->child[0] != NULL)
		child[4]->nb[12] = nb[12]->child[1];
	else{
		child[4]->nb[12] = NULL;
		child[4]->edge = 1;}
	if(nb[13]->child[0] != NULL)
		child[4]->nb[13] = nb[13]->child[3];
	else{
		child[4]->nb[13] = NULL;
		child[4]->edge = 1;}

	child[5]->nb[2] = child[6];
	child[5]->nb[3] = child[4];
	child[5]->nb[6] = child[7];
	child[5]->nb[14] = child[1];
	child[5]->nb[16] = child[2];
	child[5]->nb[17] = child[0];
	child[5]->edge = 0;
	if(nb[1]->child[0] != NULL){
		child[5]->nb[1] = nb[1]->child[4];
		child[5]->nb[5] = nb[1]->child[7];
		child[5]->nb[15] = nb[1]->child[0];}
	else{
		child[5]->nb[1] = NULL;
		child[5]->nb[5] = NULL;
		child[5]->nb[15] = NULL;
		child[5]->edge = 1;}

	if(nb[4]->child[0] != NULL){
		child[5]->nb[4] = nb[4]->child[6];
		child[5]->nb[7] = nb[4]->child[7];
		child[5]->nb[18] = nb[4]->child[2];}
	else{
		child[5]->nb[4] = NULL;
		child[5]->nb[7] = NULL;
		child[5]->nb[18] = NULL;
		child[5]->edge = 1;}

	if(nb[8]->child[0] != NULL)
		child[5]->nb[8] = nb[8]->child[7];
	else{
		child[5]->nb[8] = NULL;
		child[5]->edge = 1;}

	if(nb[9]->child[0] != NULL){
		child[5]->nb[9] = nb[9]->child[1];
		child[5]->nb[11] = nb[9]->child[2];
		child[5]->nb[12] = nb[9]->child[0];}
	else{
		child[5]->nb[9] = NULL;
		child[5]->nb[11] = NULL;
		child[5]->nb[12] = NULL;
		child[5]->edge = 1;}

	if(nb[10]->child[0] != NULL)
		child[5]->nb[10] = nb[10]->child[0];
	else{
		child[5]->nb[10] = NULL;
		child[5]->edge = 1;}
	if(nb[13]->child[0] != NULL)
		child[5]->nb[13] = nb[13]->child[2];
	else{
		child[5]->nb[13] = NULL;
		child[5]->edge = 1;}

	child[6]->nb[3] = child[7];
	child[6]->nb[4] = child[5];
	child[6]->nb[7] = child[4];
	child[6]->nb[14] = child[2];
	child[6]->nb[17] = child[3];
	child[6]->nb[18] = child[1];
	child[6]->edge = 0;
	if(nb[1]->child[0] != NULL){
		child[6]->nb[1] = nb[1]->child[7];
		child[6]->nb[8] = nb[1]->child[4];
		child[6]->nb[15] = nb[1]->child[3];}
	else{
		child[6]->nb[1] = NULL;
		child[6]->nb[8] = NULL;
		child[6]->nb[15] = NULL;
		child[6]->edge = 1;}
	if(nb[2]->child[0] != NULL){
		child[6]->nb[2] = nb[2]->child[5];
		child[6]->nb[6] = nb[2]->child[4];
		child[6]->nb[16] = nb[2]->child[1];}
	else{
		child[6]->nb[2] = NULL;
		child[6]->nb[6] = NULL;
		child[6]->nb[16] = NULL;
		child[6]->edge = 1;}
	if(nb[5]->child[0] != NULL)
		child[6]->nb[5] = nb[5]->child[4];
	else{
		child[6]->nb[5] = NULL;
		child[6]->edge = 1;}

	if(nb[11]->child[0] != NULL)
		child[6]->nb[11] = nb[11]->child[1];
	else{
		child[6]->nb[11] = NULL;
		child[6]->edge = 1;}
	if(nb[10]->child[0] != NULL)
		child[6]->nb[10] = nb[10]->child[3];
	else{
		child[6]->nb[12] = NULL;
		child[6]->edge = 1;}

	if(nb[9]->child[0] != NULL){
		child[6]->nb[9] = nb[9]->child[2];
		child[6]->nb[12] = nb[9]->child[3];
		child[6]->nb[13] = nb[9]->child[1];}
	else{
		child[6]->nb[9] = NULL;
		child[6]->nb[12] = NULL;
		child[6]->nb[13] = NULL;
		child[6]->edge = 1;}

	child[7]->nb[1] = child[6];
	child[7]->nb[4] = child[4];
	child[7]->nb[8] = child[5];
	child[7]->nb[14] = child[3];
	child[7]->nb[15] = child[2];
	child[7]->nb[18] = child[0];
	child[7]->edge = 0;
	if(nb[2]->child[0] != NULL){
		child[7]->nb[2] = nb[2]->child[4];
		child[7]->nb[5] = nb[2]->child[5];
		child[7]->nb[16] = nb[2]->child[0];}
	else{
		child[7]->nb[2] = NULL;
		child[7]->nb[5] = NULL;
		child[7]->nb[16] = NULL;
		child[7]->edge = 1;}
	if(nb[3]->child[0] != NULL){
		child[7]->nb[3] = nb[3]->child[6];
		child[7]->nb[7] = nb[3]->child[5];
		child[7]->nb[17] = nb[3]->child[2];}
	else{
		child[7]->nb[3] = NULL;
		child[7]->nb[7] = NULL;
		child[7]->nb[17] = NULL;
		child[7]->edge = 1;}
	if(nb[6]->child[0] != NULL)
		child[7]->nb[6] = nb[6]->child[5];
	else{
		child[7]->nb[6] = NULL;
		child[7]->edge = 1;}

	if(nb[11]->child[0] != NULL)
		child[7]->nb[11] = nb[11]->child[0];
	else{
		child[7]->nb[11] = NULL;
		child[7]->edge = 1;}
	if(nb[12]->child[0] != NULL)
		child[7]->nb[12] = nb[12]->child[2];
	else{
		child[7]->nb[12] = NULL;
		child[7]->edge = 1;}

	if(nb[9]->child[0] != NULL){
		child[7]->nb[9] = nb[9]->child[3];
		child[7]->nb[10] = nb[9]->child[2];
		child[7]->nb[13] = nb[9]->child[0];}
	else{
		child[7]->nb[9] = NULL;
		child[7]->nb[10] = NULL;
		child[7]->nb[13] = NULL;
		child[7]->edge = 1;}
}

//initialize u0,u1,rho,f=feq of node
void Node::Initialize(double u0_ini,double u1_ini,double u2_ini,double rho_ini,int x,int y,int z){
    for(int i = 0; i<8; i++)
        child[i] = NULL;
    level = 0;
    xcoord = x;
    ycoord = y;
    zcoord = z;
    u0 = u0_ini;
    u1 = u1_ini;
    u2 = u2_ini;
    usqr = u0*u0+u1*u1+u2*u2;
    rho = rho_ini;
    f.resize(19,0);
    feq.resize(19,0);
    if(image == CB1 || image == CB2) delta.resize(19,0);
	nb.resize(19);
	if(MRT == true){
	m.resize(19,0);
	meq.resize(19,0);
	}
	nb[0] = NULL;
	u0avg = 0;
	u1avg = 0;
	u2avg = 0;
	OmegaStar = omega[level];
    double uxy;

	if(MRT == true){
		meq[1] = -11.0*rho+19.0*usqr;
		meq[2] = -475.0/63.0*usqr; //epsilon_eq (uses rho, Yu)
		meq[4] = -2.0/3.0*u0;//qx_eq
		meq[6] = -2.0/3.0*u1;//qx_eq
		meq[8] = -2.0/3.0*u2;//qx_eq
		meq[9] = (2.0*u0*u0-(u1*u1+u2*u2));///3.0;//pxx_eq
		meq[10] = 0.0;//pixx
		meq[11] = u1*u1-u2*u2;//pww_eq
		meq[12] = 0.0;//piww
		meq[13] = u0*u1;//pxy_eq
		meq[14] = u1*u2;//pyz_eq
		meq[15] = u0*u2;//pxz_eq
		meq[16] = 0.0;//mx_eq
		meq[17] = 0.0;//my_eq
		meq[18] = 0.0;//mz_eq
		m[1] = meq[1];
		m[2] = meq[2];
		m[4] = meq[4];
		m[6] = meq[6];
		m[8] = meq[8];
		m[9] = meq[9];
		m[10] = meq[10];
		m[11] = meq[11];
		m[12] = meq[12];
		m[13] = meq[13];
		m[14] = meq[14];
		m[15] = meq[15];
		m[16] = meq[16];
		m[17] = meq[17];
		m[18] = meq[18];
		for(int i = 0; i<19; i++){
			f[i] += M_inv[i][0]*rho;
			f[i] += M_inv[i][1]*m[1];
			f[i] += M_inv[i][2]*m[2];
			f[i] += M_inv[i][3]*u0;
			f[i] += M_inv[i][4]*m[4];
			f[i] += M_inv[i][5]*u1;
			f[i] += M_inv[i][6]*m[6];
			f[i] += M_inv[i][7]*u2;
			f[i] += M_inv[i][8]*m[8];
			f[i] += M_inv[i][9]*m[9];
			f[i] += M_inv[i][10]*m[10];
			f[i] += M_inv[i][11]*m[11];
			f[i] += M_inv[i][12]*m[12];
			f[i] += M_inv[i][13]*m[13];
			f[i] += M_inv[i][14]*m[14];
			f[i] += M_inv[i][15]*m[15];
			f[i] += M_inv[i][16]*m[16];
			f[i] += M_inv[i][17]*m[17];
			f[i] += M_inv[i][18]*m[18];
		}
	}
	else{
		for(int i = 0; i<19; i++){ 
			uxy = u0*v[i][0]+u1*v[i][1]+u2*v[i][2];
			if(D2Q9i == true){
				feq[i] = t[i]*(rho+3.0*uxy+4.5*uxy*uxy-1.5*usqr);//D2Q9i
			}
			else
				feq[i] = t[i]*rho*(1.0+3.0*uxy+4.5*uxy*uxy-1.5*usqr);//D2Q9
			f[i] = feq[i];
		}
	}
}

void Node::SetImage(int Image_val)
{
    image = Image_val;
}

void Node::ComputeM()
{
	if(D2Q9i == true && MRT == true){//MRT -- compute moments
		m[0] = meq[0];
		m[3] = meq[3];
		m[5] = meq[5];
		m[7] = meq[7];
		m[1]  = -30*f[0]+-11*f[1]+-11*f[2]+-11*f[3]+-11*f[4]+  8*f[5]+  8*f[6]+  8*f[7]+  8*f[8]+-11*f[9]+  8*f[10]+  8*f[11]+  8*f[12]+  8*f[13]+-11*f[14]+  8*f[15]+  8*f[16]+  8*f[17]+  8*f[18];
		m[2]  =  12*f[0]+ -4*f[1]+ -4*f[2]+ -4*f[3]+ -4*f[4]+    f[5]+    f[6]+    f[7]+  1*f[8]+ -4*f[9]+    f[10]+  1*f[11]+    f[12]+    f[13]+ -4*f[14]+    f[15]+    f[16]+    f[17]+    f[18];
		m[4]  =           -4*f[1]         +  4*f[3]         +    f[5]+ -  f[6]+ -  f[7]+    f[8]         +    f[10]          + -  f[12]                    +    f[15]          + -  f[17]          ;
		m[6]  =                    -4*f[2]         +  4*f[4]+    f[5]+    f[6]+ -  f[7]+ -  f[8]                   +    f[11]          + -  f[13]                    +    f[16]          + -  f[18];
		m[8]  =                                                                                 + -4*f[9]+    f[10]+    f[11]+    f[12]+    f[13]+  4*f[14]+ -  f[15]+ -  f[16]+ -  f[17]+ -  f[18];
		m[9]  =            2*f[1]+ -  f[2]+  2*f[3]+ -  f[4]+    f[5]+    f[6]+    f[7]+    f[8]+ -  f[9]+    f[10]+ -2*f[11]+    f[12]+ -2*f[13]+ -  f[14]+    f[15]+ -2*f[16]+    f[17]+ -2*f[18];
		m[10] =           -4*f[1]+  2*f[2]+ -4*f[3]+  2*f[4]+    f[5]+    f[6]+    f[7]+    f[8]+  2*f[9]+    f[10]+ -2*f[11]+    f[12]+ -2*f[13]+  2*f[14]+    f[15]+ -2*f[16]+    f[17]+ -2*f[18];
		m[11] =                       f[2]         +    f[4]+    f[5]+    f[6]+    f[7]+    f[8]+ -  f[9]+ -  f[10]          + -  f[12]          + -  f[14]+ -  f[15]          + -  f[17]          ;
		m[12] =                    -2*f[2]           -2*f[4]+    f[5]+    f[6]+    f[7]+    f[8]+  2*f[9]+ -  f[10]          + -  f[12]          +  2*f[14]+ -  f[15]          + -  f[17]          ;
		m[13] =                                                  f[5]+ -  f[6]+    f[7]+ -  f[8]                                                                                                   ;
		m[14] =                                                                                                         f[11]          + -  f[13]                    + -  f[16]          +    f[18];
		m[15] =                                                                                               f[10]          + -  f[12]                    + -  f[15]          +    f[17]          ;  
		m[16] =                                                  f[5]+ -  f[6]+ -  f[7]+    f[8]           -  f[10]          +    f[12]                    + -  f[15]          +    f[17]          ;  
		m[17] =                                               -  f[5]+ -  f[6]+    f[7]+    f[8]                   +    f[11]          + -  f[13]                    +    f[16]          + -  f[18];  
		m[18] =                                                                                               f[10]+ -  f[11]+    f[12]+ -  f[13]          + -  f[15]+    f[16]+ -  f[17]+    f[18];
	}
}



void Node::MToF()
{
f[ 0] = 0.052631579*m[0]+  -0.012531328*m[1]+  0.047619048*m[2]                                                                        ;
f[ 1] = 0.052631579*m[0]+ -0.0045948204*m[1]+ -0.015873016*m[2]+  0.1*m[3]+   -0.1*m[4]                                                ;
f[ 2] = 0.052631579*m[0]+ -0.0045948204*m[1]+ -0.015873016*m[2]                        +  0.1*m[5]+   -0.1*m[6]                        ;
f[ 3] = 0.052631579*m[0]+ -0.0045948204*m[1]+ -0.015873016*m[2]+ -0.1*m[3]+    0.1*m[4]                                                ;
f[ 4] = 0.052631579*m[0]+ -0.0045948204*m[1]+ -0.015873016*m[2]                        + -0.1*m[5]+    0.1*m[6]                        ;
f[ 5] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+  0.1*m[3]+  0.025*m[4]+  0.1*m[5]+  0.025*m[6]                        ;
f[ 6] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+ -0.1*m[3]+ -0.025*m[4]+  0.1*m[5]+  0.025*m[6]                        ;
f[ 7] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+ -0.1*m[3]+ -0.025*m[4]+ -0.1*m[5]+ -0.025*m[6]                        ;
f[ 8] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+  0.1*m[3]+  0.025*m[4]+ -0.1*m[5]+ -0.025*m[6]                        ;
f[ 9] = 0.052631579*m[0]+ -0.0045948204*m[1]+ -0.015873016*m[2]                                                +  0.1*m[7]+   -0.1*m[8];
f[10] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+  0.1*m[3]+  0.025*m[4]                        +  0.1*m[7]+  0.025*m[8];
f[11] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]                        +  0.1*m[5]+  0.025*m[6]+  0.1*m[7]+  0.025*m[8];
f[12] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+ -0.1*m[3]+ -0.025*m[4]                        +  0.1*m[7]+  0.025*m[8];
f[13] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]                        + -0.1*m[5]+ -0.025*m[6]+  0.1*m[7]+  0.025*m[8];
f[14] = 0.052631579*m[0]+ -0.0045948204*m[1]+ -0.015873016*m[2]                                                + -0.1*m[7]+    0.1*m[8];
f[15] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+  0.1*m[3]+  0.025*m[4]                        + -0.1*m[7]+ -0.025*m[8];
f[16] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]                        +  0.1*m[5]+  0.025*m[6]+ -0.1*m[7]+ -0.025*m[8];
f[17] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+ -0.1*m[3]+ -0.025*m[4]                        + -0.1*m[7]+ -0.025*m[8];
f[18] = 0.052631579*m[0]+  0.0033416876*m[1]+  0.003968254*m[2]+    0*m[3]+      0*m[4]+ -0.1*m[5]+ -0.025*m[6]+ -0.1*m[7]+ -0.025*m[8];

//feq[ 0] +=            0*meq[9]+            0*meq[10]+            0*meq[11]+            0*meq[12]+     0*meq[13]+     0*meq[14]+     0*meq[15]+      0*meq[16]+      0*meq[17]+      0*meq[18];
f[ 1] +=  0.055555556*m[9]+ -0.055555556*m[10];//+            meq[11]+            meq[12]+     meq[13]+     meq[14]+     meq[15]+      meq[16]+      meq[17]+      meq[18];
f[ 2] += -0.027777778*m[9]+  0.027777778*m[10]+  0.083333333*m[11]+ -0.083333333*m[12];//+     meq[13]+     meq[14]+     meq[15]+      meq[16]+      meq[17]+      meq[18];
f[ 3] +=  0.055555556*m[9]+ -0.055555556*m[10];//+            meq[11]+            meq[12]+     meq[13]+     meq[14]+     meq[15]+      meq[16]+      meq[17]+      meq[18];
f[ 4] += -0.027777778*m[9]+  0.027777778*m[10]+  0.083333333*m[11]+ -0.083333333*m[12];//+     meq[13]+     meq[14]+     meq[15]+      meq[16]+      meq[17]+      meq[18];
f[ 5] +=  0.027777778*m[9]+  0.013888889*m[10]+  0.083333333*m[11]+  0.041666667*m[12]+  0.25*m[13]                          +  0.125*m[16]+ -0.125*m[17]              ;
f[ 6] +=  0.027777778*m[9]+  0.013888889*m[10]+  0.083333333*m[11]+  0.041666667*m[12]+ -0.25*m[13]                          + -0.125*m[16]+ -0.125*m[17]              ;
f[ 7] +=  0.027777778*m[9]+  0.013888889*m[10]+  0.083333333*m[11]+  0.041666667*m[12]+  0.25*m[13]                          + -0.125*m[16]+  0.125*m[17]              ;
f[ 8] +=  0.027777778*m[9]+  0.013888889*m[10]+  0.083333333*m[11]+  0.041666667*m[12]+ -0.25*m[13]                          +  0.125*m[16]+  0.125*m[17]              ;
f[ 9] += -0.027777778*m[9]+  0.027777778*m[10]+ -0.083333333*m[11]+  0.083333333*m[12]                                                                                 ;
f[10] +=  0.027777778*m[9]+  0.013888889*m[10]+ -0.083333333*m[11]+ -0.041666667*m[12]                          +  0.25*m[15]+ -0.125*m[16]              +  0.125*m[18];
f[11] += -0.055555556*m[9]+ -0.027777778*m[10]                                                     +  0.25*m[14]                           +  0.125*m[17]+ -0.125*m[18];
f[12] +=  0.027777778*m[9]+  0.013888889*m[10]+ -0.083333333*m[11]+ -0.041666667*m[12]                          + -0.25*m[15]+  0.125*m[16]              +  0.125*m[18];
f[13] += -0.055555556*m[9]+ -0.027777778*m[10]                                                       -0.25*m[14]                           + -0.125*m[17]+ -0.125*m[18];
f[14] += -0.027777778*m[9]+  0.027777778*m[10]+ -0.083333333*m[11]+  0.083333333*m[12]                                                                                 ;
f[15] +=  0.027777778*m[9]+  0.013888889*m[10]+ -0.083333333*m[11]+ -0.041666667*m[12]                          + -0.25*m[15]+ -0.125*m[16]              + -0.125*m[18];
f[16] += -0.055555556*m[9]+ -0.027777778*m[10]                                                     + -0.25*m[14]                           +  0.125*m[17]+  0.125*m[18];
f[17] +=  0.027777778*m[9]+  0.013888889*m[10]+ -0.083333333*m[11]+ -0.041666667*m[12]                          +  0.25*m[15]+  0.125*m[16]              + -0.125*m[18];
f[18] += -0.055555556*m[9]+ -0.027777778*m[10]                                                     +  0.25*m[14]                           + -0.125*m[17]+  0.125*m[18];
}





void Node::ComputeFeq()
{
    double uxy,delx;
	if(D2Q9i == true && MRT == true){//MRT -- compute moments
		meq[0] = rho;
		meq[3] = u0*rho;
		meq[5] = u1*rho;
		meq[7] = u2*rho;
		meq[1] = -11.0*rho+19.0*(meq[3]*meq[3]+meq[5]*meq[5]+meq[7]*meq[7]);
		meq[2] = -475.0/63.0*(u0*u0+u1*u1+u2*u2); //epsilon_eq (uses rho, Yu)
		meq[4] = -2.0/3.0*meq[3];//qx_eq
		meq[6] = -2.0/3.0*meq[5];//qx_eq
		meq[8] = -2.0/3.0*meq[7];//qx_eq
		meq[9] = (2.0*meq[3]*meq[3]-(meq[5]*meq[5]+meq[7]*meq[7]));//(2.0*u0*u0-(u1*u1+u2*u2));///3.0;//pxx_eq
		meq[10] = 0.0;//-0.5*meq[9];//0.0;//-0.5*meq[9];//0.0;//pixx
		meq[11] = (meq[5]*meq[5]-meq[7]*meq[7]);//pww_eq
		meq[12] = 0.0;//-0.5*meq[11];//0.0;//-0.5*meq[9];//0.0;//piww
		meq[13] = meq[3]*meq[5];//pxy_eq
		meq[14] = meq[5]*meq[7];//pyz_eq
		meq[15] = meq[3]*meq[7];//pxz_eq
		meq[16] = 0.0;//mx_eq
		meq[17] = 0.0;//my_eq
		meq[18] = 0.0;//mz_eq
		m[0] = meq[0];
		m[3] = meq[3];
		m[5] = meq[5];
		m[7] = meq[7];

		m[1]  = -30*f[0]+-11*f[1]+-11*f[2]+-11*f[3]+-11*f[4]+  8*f[5]+  8*f[6]+  8*f[7]+  8*f[8]+-11*f[9]+  8*f[10]+  8*f[11]+  8*f[12]+  8*f[13]+-11*f[14]+  8*f[15]+  8*f[16]+  8*f[17]+  8*f[18];
		m[2]  =  12*f[0]+ -4*f[1]+ -4*f[2]+ -4*f[3]+ -4*f[4]+    f[5]+    f[6]+    f[7]+  1*f[8]+ -4*f[9]+    f[10]+  1*f[11]+    f[12]+    f[13]+ -4*f[14]+    f[15]+    f[16]+    f[17]+    f[18];
		m[4]  =           -4*f[1]         +  4*f[3]         +    f[5]+ -  f[6]+ -  f[7]+    f[8]         +    f[10]          + -  f[12]                    +    f[15]          + -  f[17]          ;
		m[6]  =                    -4*f[2]         +  4*f[4]+    f[5]+    f[6]+ -  f[7]+ -  f[8]                   +    f[11]          + -  f[13]                    +    f[16]          + -  f[18];
		m[8]  =                                                                                 + -4*f[9]+    f[10]+    f[11]+    f[12]+    f[13]+  4*f[14]+ -  f[15]+ -  f[16]+ -  f[17]+ -  f[18];
		m[9]  =            2*f[1]+ -  f[2]+  2*f[3]+ -  f[4]+    f[5]+    f[6]+    f[7]+    f[8]+ -  f[9]+    f[10]+ -2*f[11]+    f[12]+ -2*f[13]+ -  f[14]+    f[15]+ -2*f[16]+    f[17]+ -2*f[18];
		m[10] =           -4*f[1]+  2*f[2]+ -4*f[3]+  2*f[4]+    f[5]+    f[6]+    f[7]+    f[8]+  2*f[9]+    f[10]+ -2*f[11]+    f[12]+ -2*f[13]+  2*f[14]+    f[15]+ -2*f[16]+    f[17]+ -2*f[18];
		m[11] =                       f[2]         +    f[4]+    f[5]+    f[6]+    f[7]+    f[8]+ -  f[9]+ -  f[10]          + -  f[12]          + -  f[14]+ -  f[15]          + -  f[17]          ;
		m[12] =                    -2*f[2]           -2*f[4]+    f[5]+    f[6]+    f[7]+    f[8]+  2*f[9]+ -  f[10]          + -  f[12]          +  2*f[14]+ -  f[15]          + -  f[17]          ;
		m[13] =                                                  f[5]+ -  f[6]+    f[7]+ -  f[8]                                                                                                   ;
		m[14] =                                                                                                         f[11]          + -  f[13]                    + -  f[16]          +    f[18];
		m[15] =                                                                                               f[10]          + -  f[12]                    + -  f[15]          +    f[17]          ;  
		m[16] =                                                  f[5]+ -  f[6]+ -  f[7]+    f[8]           -  f[10]          +    f[12]                    + -  f[15]          +    f[17]          ;  
		m[17] =                                               -  f[5]+ -  f[6]+    f[7]+    f[8]                   +    f[11]          + -  f[13]                    +    f[16]          + -  f[18];  
		m[18] =                                                                                               f[10]+ -  f[11]+    f[12]+ -  f[13]          + -  f[15]+    f[16]+ -  f[17]+    f[18];

	}
	else{
	    for(int i = 0; i<19; i++){  
			uxy = u0*v[i][0]+u1*v[i][1]+u2*v[i][2];
			if(D2Q9i == true){
				feq[i] = t[i]*(rho+3.0*uxy+4.5*uxy*uxy-1.5*usqr);//D2Q9i
			}
			else{
				feq[i] = t[i]*rho*(1.0+3.0*uxy+4.5*uxy*uxy-1.5*usqr);//D2Q9
			}
	    }
	}

	if(MRT == true){
		delx = (pow(2,level));
		//PI is the Strain rate tensor
		double PI11 = -1.0/38.0*(S1[level]*(m[1]-meq[1])+19.0*omega[level]*(m[9]-meq[9]))*delx;
		double PI22 = -1.0/76.0*(2.0*S1[level]*(m[1]-meq[1])-19.0*(omega[level]*(m[9]-meq[9])-3.0*omega[level]*(m[11]-meq[11])))*delx;
		double PI33 = -1.0/76.0*(2.0*S1[level]*(m[1]-meq[1])-19.0*(omega[level]*(m[9]-meq[9])+3.0*omega[level]*(m[11]-meq[11])))*delx;
		double PI12 = -1.5*omega[level]*(m[13]-meq[13])*delx;
		double PI23 = -1.5*omega[level]*(m[14]-meq[14])*delx;
		double PI13 = -1.5*omega[level]*(m[15]-meq[15])*delx;
		nu0 = ((1.0/omega[level])-0.5)/3.0;
		Smag = sqrt(PI11*PI11+PI22*PI22+PI33*PI33+PI12*PI12+PI23*PI23+PI13*PI13);
	}
	if(LES == true){
		Cs = 0.01;//0.0225;//1;
		if(xcoord > (xDim-10)) Cs = 0.0225;
		OmegaStar = 1.0/(3.0*(nu0+delx*delx*Cs*Smag)+0.5);
	}
}

void Node::Collide()
{
    if(MRT == true){
		if (image != BB && image != BB1 && image != BB2 && image != CB1 &&
			image != CB2){// && image != NeumannEast)
double Omega;
if(LES == true) Omega = OmegaStar;
else{
	Omega = omega[level];
}
double s_1,s_2,s_4,s_10,s_16;
s_1  = S1 [level];
s_2  = S2 [level];
s_4  = S4 [level];
s_10 = S10[level];
s_16 = S16[level];

f[ 0] -= - 0.012531328*(m[1]-meq[1])*s_1 +  0.047619048*(m[2]-meq[2])*s_2 ;
f[ 1] -= -0.0045948204*(m[1]-meq[1])*s_1 + -0.015873016*(m[2]-meq[2])*s_2 +   -0.1*(m[4]-meq[4])*s_4                                                       +  0.055555556*(m[9]-meq[9])*Omega + -0.055555556*(m[10]-meq[10])*s_10;
f[ 2] -= -0.0045948204*(m[1]-meq[1])*s_1 + -0.015873016*(m[2]-meq[2])*s_2                            +   -0.1*(m[6]-meq[6])*s_4                            + -0.027777778*(m[9]-meq[9])*Omega +  0.027777778*(m[10]-meq[10])*s_10;
f[ 3] -= -0.0045948204*(m[1]-meq[1])*s_1 + -0.015873016*(m[2]-meq[2])*s_2 +    0.1*(m[4]-meq[4])*s_4                                                       +  0.055555556*(m[9]-meq[9])*Omega + -0.055555556*(m[10]-meq[10])*s_10;
f[ 4] -= -0.0045948204*(m[1]-meq[1])*s_1 + -0.015873016*(m[2]-meq[2])*s_2                            +    0.1*(m[6]-meq[6])*s_4                            + -0.027777778*(m[9]-meq[9])*Omega +  0.027777778*(m[10]-meq[10])*s_10;
f[ 5] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2 +  0.025*(m[4]-meq[4])*s_4 +  0.025*(m[6]-meq[6])*s_4                            +  0.027777778*(m[9]-meq[9])*Omega +  0.013888889*(m[10]-meq[10])*s_10;
f[ 6] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2 + -0.025*(m[4]-meq[4])*s_4 +  0.025*(m[6]-meq[6])*s_4                            +  0.027777778*(m[9]-meq[9])*Omega +  0.013888889*(m[10]-meq[10])*s_10;
f[ 7] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2 + -0.025*(m[4]-meq[4])*s_4 + -0.025*(m[6]-meq[6])*s_4                            +  0.027777778*(m[9]-meq[9])*Omega +  0.013888889*(m[10]-meq[10])*s_10;
f[ 8] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2 +  0.025*(m[4]-meq[4])*s_4 + -0.025*(m[6]-meq[6])*s_4                            +  0.027777778*(m[9]-meq[9])*Omega +  0.013888889*(m[10]-meq[10])*s_10;
f[ 9] -= -0.0045948204*(m[1]-meq[1])*s_1 + -0.015873016*(m[2]-meq[2])*s_2                                                       +   -0.1*(m[8]-meq[8])*s_4 + -0.027777778*(m[9]-meq[9])*Omega +  0.027777778*(m[10]-meq[10])*s_10;
f[10] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2 +  0.025*(m[4]-meq[4])*s_4                            +  0.025*(m[8]-meq[8])*s_4 +  0.027777778*(m[9]-meq[9])*Omega +  0.013888889*(m[10]-meq[10])*s_10;
f[11] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2                            +  0.025*(m[6]-meq[6])*s_4 +  0.025*(m[8]-meq[8])*s_4 + -0.055555556*(m[9]-meq[9])*Omega + -0.027777778*(m[10]-meq[10])*s_10;
f[12] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2 + -0.025*(m[4]-meq[4])*s_4                            +  0.025*(m[8]-meq[8])*s_4 +  0.027777778*(m[9]-meq[9])*Omega +  0.013888889*(m[10]-meq[10])*s_10;
f[13] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2                            + -0.025*(m[6]-meq[6])*s_4 +  0.025*(m[8]-meq[8])*s_4 + -0.055555556*(m[9]-meq[9])*Omega + -0.027777778*(m[10]-meq[10])*s_10;
f[14] -= -0.0045948204*(m[1]-meq[1])*s_1 + -0.015873016*(m[2]-meq[2])*s_2                                                       +    0.1*(m[8]-meq[8])*s_4 + -0.027777778*(m[9]-meq[9])*Omega +  0.027777778*(m[10]-meq[10])*s_10;
f[15] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2 +  0.025*(m[4]-meq[4])*s_4                            + -0.025*(m[8]-meq[8])*s_4 +  0.027777778*(m[9]-meq[9])*Omega +  0.013888889*(m[10]-meq[10])*s_10;
f[16] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2                            +  0.025*(m[6]-meq[6])*s_4 + -0.025*(m[8]-meq[8])*s_4 + -0.055555556*(m[9]-meq[9])*Omega + -0.027777778*(m[10]-meq[10])*s_10;
f[17] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2 + -0.025*(m[4]-meq[4])*s_4                            + -0.025*(m[8]-meq[8])*s_4 +  0.027777778*(m[9]-meq[9])*Omega +  0.013888889*(m[10]-meq[10])*s_10;
f[18] -=  0.0033416876*(m[1]-meq[1])*s_1 +  0.003968254*(m[2]-meq[2])*s_2                            + -0.025*(m[6]-meq[6])*s_4 + -0.025*(m[8]-meq[8])*s_4 + -0.055555556*(m[9]-meq[9])*Omega + -0.027777778*(m[10]-meq[10])*s_10;

//f[ 0] -=            0*(m[11]-meq[11])*Omega +            0*(m[12]-meq[12])*s_10 + (    0*(m[13]-meq[13]) +     0*(m[14]-meq[14]) +     0*(m[15]-meq[15]))*Omega ;
//f[ 1] -=            0*(m[11]-meq[11])*Omega +            0*(m[12]-meq[12])*s_10 + (    0*(m[13]-meq[13]) +     0*(m[14]-meq[14]) +     0*(m[15]-meq[15]))*Omega ;
f[ 2] -=  0.083333333*(m[11]-meq[11])*Omega + -0.083333333*(m[12]-meq[12])*s_10 ;
//f[ 3] -=            0*(m[11]-meq[11])*Omega +            0*(m[12]-meq[12])*s_10 ;
f[ 4] -=  0.083333333*(m[11]-meq[11])*Omega + -0.083333333*(m[12]-meq[12])*s_10 ;
f[ 5] -=  0.083333333*(m[11]-meq[11])*Omega +  0.041666667*(m[12]-meq[12])*s_10 + ( 0.25*(m[13]-meq[13])                                                )*Omega; 
f[ 6] -=  0.083333333*(m[11]-meq[11])*Omega +  0.041666667*(m[12]-meq[12])*s_10 + (-0.25*(m[13]-meq[13])                                                )*Omega; 
f[ 7] -=  0.083333333*(m[11]-meq[11])*Omega +  0.041666667*(m[12]-meq[12])*s_10 + ( 0.25*(m[13]-meq[13])                                                )*Omega; 
f[ 8] -=  0.083333333*(m[11]-meq[11])*Omega +  0.041666667*(m[12]-meq[12])*s_10 + (-0.25*(m[13]-meq[13])                                                )*Omega; 
f[ 9] -= -0.083333333*(m[11]-meq[11])*Omega +  0.083333333*(m[12]-meq[12])*s_10 ;
f[10] -= -0.083333333*(m[11]-meq[11])*Omega + -0.041666667*(m[12]-meq[12])*s_10  +(                                              +  0.25*(m[15]-meq[15]))*Omega ;
f[11] -=                                                                         +(                         0.25*(m[14]-meq[14])                        )*Omega ;
f[12] -= -0.083333333*(m[11]-meq[11])*Omega + -0.041666667*(m[12]-meq[12])*s_10  +(                                              + -0.25*(m[15]-meq[15]))*Omega ;
f[13] -=                                                                         +(                        -0.25*(m[14]-meq[14])                        )*Omega ;
f[14] -= -0.083333333*(m[11]-meq[11])*Omega +  0.083333333*(m[12]-meq[12])*s_10 ;
f[15] -= -0.083333333*(m[11]-meq[11])*Omega + -0.041666667*(m[12]-meq[12])*s_10  +(                                              + -0.25*(m[15]-meq[15]))*Omega ;
f[16] -=                                                                         +(                        -0.25*(m[14]-meq[14])                        )*Omega ;
f[17] -= -0.083333333*(m[11]-meq[11])*Omega + -0.041666667*(m[12]-meq[12])*s_10  +(                                              +  0.25*(m[15]-meq[15]))*Omega ;
f[18] -=                                                                         +(                         0.25*(m[14]-meq[14])                        )*Omega ;

f[ 5] -=  0.125*(m[16]-meq[16])*s_16 + -0.125*(m[17]-meq[17])*s_16 ;                        
f[ 6] -= -0.125*(m[16]-meq[16])*s_16 + -0.125*(m[17]-meq[17])*s_16 ;                        
f[ 7] -= -0.125*(m[16]-meq[16])*s_16 +  0.125*(m[17]-meq[17])*s_16 ;                        
f[ 8] -=  0.125*(m[16]-meq[16])*s_16 +  0.125*(m[17]-meq[17])*s_16 ;                        
f[10] -= -0.125*(m[16]-meq[16])*s_16                               +  0.125*(m[18]-meq[18])*s_16 ;
f[11] -=                             +  0.125*(m[17]-meq[17])*s_16 + -0.125*(m[18]-meq[18])*s_16 ;
f[12] -=  0.125*(m[16]-meq[16])*s_16                               +  0.125*(m[18]-meq[18])*s_16 ;
f[13] -=                             + -0.125*(m[17]-meq[17])*s_16 + -0.125*(m[18]-meq[18])*s_16 ;
f[15] -= -0.125*(m[16]-meq[16])*s_16                               + -0.125*(m[18]-meq[18])*s_16 ;
f[16] -=                             +  0.125*(m[17]-meq[17])*s_16 +  0.125*(m[18]-meq[18])*s_16 ;
f[17] -=  0.125*(m[16]-meq[16])*s_16                               + -0.125*(m[18]-meq[18])*s_16 ;
f[18] -=                             + -0.125*(m[17]-meq[17])*s_16 +  0.125*(m[18]-meq[18])*s_16 ;

feqflag = 0;

		}
	}
	else{
	for(int i = 0; i<19; i++){
		if (image != BB && image != BB1 && image != BB2 && image != CB1 &&
			image != CB2){
		    f[i] = (1.0-omega[level])*f[i]+omega[level]*feq[i];
			f[i] += t[i]*3.0*dPdx*v[i][0];}
//			image != CB2 && image != xsymmetry && image != ysymmetry)
    }
	}
}


void Node::ComputeMacros()
{
    rho = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8]
			+f[9]+f[10]+f[11]+f[12]+f[13]
			+f[14]+f[15]+f[16]+f[17]+f[18];

	if(D2Q9i == true){
    	u0 = (f[1]-f[3]+f[5]-f[6]-f[7]+f[8]
				+f[10]-f[12]+f[15]-f[17]);//D2Q9i
    	u1 = (f[2]-f[4]+f[5]+f[6]-f[7]-f[8]//rho;
				+f[11]-f[13]+f[16]-f[18]);//D2Q9i
		u2 = (f[9]+f[10]+f[11]+f[12]+f[13]
				-f[14]-f[15]-f[16]-f[17]-f[18]);
	}
	else{
    	u0 = (f[1]-f[3]+f[5]-f[6]-f[7]+f[8])/rho;//D2Q9
    	u1 = (f[2]-f[4]+f[5]+f[6]-f[7]-f[8])/rho;
	}
    usqr = u0*u0+u1*u1+u2*u2;

    if (usqr > 10){
        ErrorFlag = 1;
        cout<<"Velocity explosion"<<xcoord<<", "<<ycoord<<", "<<zcoord<<"\n";
    }
    if (rho < 0){
        ErrorFlag = 1;
        cout<<"Negative density at"<<xcoord<<", "<<ycoord<<", "<<zcoord<<"\n";
    }
}

void Node::VelAverage(int tStep)//compute time averaged velocity
{
	u0avg = u0avg*(tStep-StartAvg)/(tStep-StartAvg+1)+u0/(tStep-StartAvg+1);
	u1avg = u1avg*(tStep-StartAvg)/(tStep-StartAvg+1)+u1/(tStep-StartAvg+1);
	u2avg = u2avg*(tStep-StartAvg)/(tStep-StartAvg+1)+u2/(tStep-StartAvg+1);
}

void Node::BoundaryMacro()
{
    if (image == DirichletNorth){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
        u0 = uMax;
        u1 = 0.0;
        
        fInt1 = f[0]+f[1]+f[3];
        fInt2 = f[6]+f[2]+f[5];
		if(D2Q9i == true)
        	rho = -u1+(fInt1+2.0*fInt2); //D2Q9i
		else
			rho = (fInt1+2.0*fInt2)/(1.0+u1);
        

    }
    else if (image == DirichletSouth){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
        u0 = uMax;
        u1 = 0;

        if((xcoord == 0)||(xcoord == xDim-1))
            u0 = 0;
        
        fInt1 = f[0]+f[1]+f[3];
        fInt2 = f[4]+f[7]+f[8];
		if(D2Q9i == true)
        	rho = u1-(fInt1+2.0*fInt2); //D2Q9i
        else
			rho = u1-(fInt1+2.0*fInt2);


    }
    else if (image == DirichletWest || image == DirichletWest2){

		if(ycoord == 0){
			f[2] = f[4];
			f[6] = f[7];
			f[11] = f[13];
			f[16] = f[18];
		}
		if(ycoord == yDim-1){
			f[4] = f[2];
			f[7] = f[6];
			f[13] = f[11];
			f[18] = f[16];
		}
		if(zcoord == 0){
			f[9] = f[14];
			f[10] = f[15];
			f[11] = f[16];
			f[12] = f[17];
			f[13] = f[18];
		}
		if(zcoord == zDim-1){
			f[14] = f[9];
			f[15] = f[10];
			f[16] = f[11];
			f[17] = f[12];
			f[18] = f[13];
		}

		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
        u0 = uMax;//*PoisProf(zcoord)*1.5;
        u1 = 0.0;//uMax;//0.0;
		u2 = 0.0;
        
		fInt1 = f[0]+f[2]+f[4]+f[9]+f[11]+f[13]+f[14]+f[16]+f[18];
        fInt2 = f[3]+f[6]+f[7]+f[12]+f[17];
		if(D2Q9i == true)
        	rho = u0+(fInt1+2.0*fInt2); //D2Q9i
		else
			rho = u0-(fInt1+2.0*fInt2)/(1.0-u0);
    }   
    else if (image == DirichletEast){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
        u0 = uMax;
        u1 = 0;
        
        fInt1 = f[0]+f[2]+f[4];
        fInt2 = f[1]+f[5]+f[8];
		if(D2Q9i == true)
        	rho = -u0+(fInt1+2.0*fInt2); //D2Q9i
		else
			rho = (fInt1+2.0*fInt2)/(1.0+u0);

    }
	else if (image == NeumannWest){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
        rho = 1.0;
        u1 = 0.0;
        
        fInt1 = f[0]+f[2]+f[4];
        fInt2 = f[3]+f[6]+f[7];
		if(D2Q9i == true)
        u0 = rho-(fInt1+2.0*fInt2); //D2Q9i
		else
		u0 = 1.0-(fInt1+2.0*fInt2)/rho;
        
    }   
    else if (image == NeumannEast || image == NeumannEast2){

		if(yPeriodics == 0){
		if(ycoord == 0){
			f[2] = f[4];
			f[5] = f[8];
			f[11] = f[13];
			f[16] = f[18];
		}
		if(ycoord == yDim-1){
			f[4] = f[2];
			f[8] = f[5];
			f[13] = f[11];
			f[18] = f[16];
		}
		}

		if(zPeriodics == 0){
		if(zcoord == 0){
			f[9] = f[14];
			f[10] = f[15];
			f[11] = f[16];
			f[12] = f[17];
			f[13] = f[18];
		}
		if(zcoord == zDim-1){
			f[14] = f[9];
			f[15] = f[10];
			f[16] = f[11];
			f[17] = f[12];
			f[18] = f[13];
		}
		}

		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
        rho = 1.0;
        u1 = 0.0;//nb[3]->u1;//0.0;
		u2 = 0.0;
        
        fInt1 = f[0]+f[2]+f[4]+f[9]+f[11]+f[13]+f[14]+f[16]+f[18];
        fInt2 = f[1]+f[5]+f[8]+f[10]+f[15];
		if(D2Q9i == true)
        u0 = -rho+(fInt1+2.0*fInt2); //D2Q9i
		else
		u0 = -1.0+(fInt1+2.0*fInt2)/rho;
        
    }

}

void Node::InletOutlet()
{
    if (image == DirichletNorth){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
           
        fDiff = 0.5*(f[1]-f[3]);
		if(D2Q9i == true){
        rhoUy = u1/6.0;
        rhoUx = u0/2.0;
		}
		else{
		rhoUy = rho*u1/6.0;
		rhoUx = rho*u0/2.0;
		}

        f[4] = f[2]-4.0*rhoUy;
        f[7] = f[5]+fDiff-rhoUx+rhoUy;
        f[8] = f[6]-fDiff+rhoUx+rhoUy;
    }
    else if (image == DirichletSouth){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
  
        fDiff = 0.5*(f[1]-f[3]);
		if(D2Q9i == true){
        rhoUy = u1/6.0;
        rhoUx = u0/2.0;
		}
		else{
		rhoUy = rho*u1/6.0;
		rhoUx = rho*u0/2.0;
		}

        f[2] = f[4]+4.0*rhoUy;
        f[5] = f[7]-fDiff+rhoUx+rhoUy;
        f[6] = f[8]+fDiff-rhoUx+rhoUy; 
    }
    else if (image == DirichletWest){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
           
        fDiff = 0.5*(f[2]-f[4]);
		if(D2Q9i == true){
        rhoUy = u1/2.0;
        rhoUx = u0/6.0;
		}
		else{
		rhoUy = rho*u1/2.0;
		rhoUx = rho*u0/6.0;
		}

        f[1] = f[3]+4.0*rhoUx;
        f[5] = f[7]-fDiff+rhoUy+rhoUx;
        f[8] = f[6]+fDiff-rhoUy+rhoUx;
		u0 = (f[1]-f[3]+f[5]-f[6]-f[7]+f[8]);
		
    }   
    else if (image == DirichletEast){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
            
        fDiff = 0.5*(f[2]-f[4]);
		if(D2Q9i == true){
        rhoUy = u1/2.0;
        rhoUx = u0/6.0;
		}
		else{
		rhoUy = rho*u1/2.0;
		rhoUx = rho*u0/2.0;
		}

        f[3] = f[1]-4.0*rhoUx;
        f[7] = f[5]+fDiff-rhoUy-rhoUx;
        f[6] = f[8]-fDiff+rhoUy-rhoUx; 
    }
	else if (image == NeumannWest){
		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
             
        fDiff = 0.5*(f[2]-f[4]);
		if(D2Q9i == true){
        rhoUy = u1/2.0;
        rhoUx = u0/6.0;
		}
		else{
		rhoUy = rho*u1/2.0;
		rhoUx = rho*u0/6.0;
		}

        f[1] = f[3]+4.0*rhoUy;
        f[5] = f[7]-fDiff+rhoUy+rhoUx;
        f[8] = f[6]+fDiff-rhoUy+rhoUx; 
    }   
    else if (image == NeumannEast){

		double rhoUx,rhoUy,fInt1,fInt2,fDiff;
              
        fDiff = 0.5*(f[2]-f[4]);
		if(D2Q9i == true){
        rhoUy = u1/2.0;
        rhoUx = u0/6.0;
		}
		else{
		rhoUy = rho*u1/2.0;
		rhoUx = rho*u0/6.0;
		}

        f[3] = f[1]-4.0*rhoUx;
        f[7] = f[5]+fDiff-rhoUy-rhoUx;
        f[6] = f[8]-fDiff+rhoUy-rhoUx; 
    }
}


void Node::RegularizedBoundary()
{
//	if(edge == 1){
	double PI11 = 0;
	double PI12 = 0;
	double PI22 = 0;
	double PI33 = 0;
	double PI13 = 0;
	double PI23 = 0;
	double Q11 = 0;//this Q is diff from Q used for Smagorinski Model
	double Q12 = 0;
	double Q22 = 0;
	double Q33 = 0;
	double Q23 = 0;
	double Q13 = 0;
	int check = 0;
	double uxy;
    if (image == DirichletNorth2){

    }
    else if (image == DirichletSouth2){

    }
    else if (image == DirichletWest2){
		check = 1;
	    for(int i = 0; i<19; i++){  
			uxy = u0*v[i][0]+u1*v[i][1]+u2*v[i][2];
			usqr = u0*u0+u1*u1+u2*u2;
			if(D2Q9i == true){
				feq[i] = t[i]*(rho+3.0*uxy+4.5*uxy*uxy-1.5*usqr);//D2Q9i
			}
			else{
				feq[i] = t[i]*rho*(1.0+3.0*uxy+4.5*uxy*uxy-1.5*usqr);//D2Q9
			}
	    }
		f[1] = feq[1]+f[opposite[1]]-feq[opposite[1]];
		f[5] = feq[5]+f[opposite[5]]-feq[opposite[5]];
		f[8] = feq[8]+f[opposite[8]]-feq[opposite[8]];
		f[10] = feq[10]+f[opposite[10]]-feq[opposite[10]];
		f[15] = feq[15]+f[opposite[15]]-feq[opposite[15]];
    }
    else if (image == DirichletEast2){

    }
	else if (image == NeumannWest2){
    }   
    else if (image == NeumannEast2){
		check = 1;
	    for(int i = 0; i<19; i++){  
			uxy = u0*v[i][0]+u1*v[i][1]+u2*v[i][2];
			usqr = u0*u0+u1*u1+u2*u2;
			if(D2Q9i == true){
				feq[i] = t[i]*(rho+3.0*uxy+4.5*uxy*uxy-1.5*usqr);//D2Q9i
			}
			else{
				feq[i] = t[i]*rho*(1.0+3.0*uxy+4.5*uxy*uxy-1.5*usqr);//D2Q9
			}
	    }
		f[3] = feq[3]+f[opposite[3]]-feq[opposite[3]];
		f[6] = feq[6]+f[opposite[6]]-feq[opposite[6]];
		f[7] = feq[7]+f[opposite[7]]-feq[opposite[7]];
		f[12] = feq[12]+f[opposite[12]]-feq[opposite[12]];
		f[17] = feq[17]+f[opposite[17]]-feq[opposite[17]];
    }

	if(check == 1){
		PI11 = 0;
		PI22 = 0;
		PI33 = 0;
		PI12 = 0;
		PI23 = 0;
		PI13 = 0;
		for(int i = 0; i<19; i++){
			PI11 += v[i][0]*v[i][0]*(f[i]-feq[i]);
			PI22 += v[i][1]*v[i][1]*(f[i]-feq[i]);
			PI33 += v[i][2]*v[i][2]*(f[i]-feq[i]);
			PI12 += v[i][0]*v[i][1]*(f[i]-feq[i]);
			PI23 += v[i][1]*v[i][2]*(f[i]-feq[i]);
			PI13 += v[i][0]*v[i][2]*(f[i]-feq[i]);
		}
		//regenerate f from Q and PI
		for(int i = 0; i<19; i++){
			Q11 = v[i][0]*v[i][0]-1.0/3.0;
			Q12 = v[i][1]*v[i][0];
			Q13 = v[i][2]*v[i][0];
			Q22 = v[i][1]*v[i][1]-1.0/3.0;
			Q23 = v[i][2]*v[i][1];
			Q33 = v[i][2]*v[i][2]-1.0/3.0;
			f[i] = feq[i]+t[i]*4.5*(Q11*PI11+2.0*(Q12*PI12+Q13*PI13+Q23*PI23)
									+Q22*PI22+Q33*PI33);
		}

	}
}





void Node::BounceBack()
{
    if (image >= BB && image <= BB3)
    {
		for(int k = 0; k<9; k++){
			swap(f[halfstream[k]],f[opposite[halfstream[k]]]);
		}
    }
}
void Node::Symmetry()
{
	if (image == zsymmetry){
		if(zcoord == zDim-1){
			f[14] = f[ 9];
			f[15] = f[10];
			f[16] = f[11];
			f[17] = f[12];
			f[18] = f[13];
		}
		else if(zcoord == 0){//== 0){
			f[ 9] = f[14];
			f[10] = f[15];
			f[11] = f[16];
			f[12] = f[17];
			f[13] = f[18];
		}
	}
	if (image == xsymmetry){
		if(xcoord == xDim-1){
			f[1] = f[xsym[1]];
			f[5] = f[xsym[5]];
			f[8] = f[xsym[8]];
			f[10] = f[xsym[10]];
			f[15] = f[xsym[15]];
		}
		else if(xcoord == 0){//== 0){
			f[3] = f[xsym[3]];
			f[6] = f[xsym[6]];
			f[7] = f[xsym[7]];
			f[12] = f[xsym[12]];
			f[17] = f[xsym[17]];
		}
	}
	if (image == ysymmetry){
		//if(ycoord == yDim-1){
		if(ycoord > yDim-2){
			f[4 ] = f[2 ];
			f[7 ] = f[6 ];
			f[8 ] = f[5 ];
			f[13] = f[11];
			f[18] = f[16];
		}
		//else if(ycoord == 0){//== 0){
		else if(ycoord < 1){//== 0){
			f[2 ] = f[4 ];
			f[5 ] = f[8 ];
			f[6 ] = f[7 ];
			f[11] = f[13];
			f[16] = f[18];
		}
	}
}






void Node::Periodic(Node * Pair)
{
	if (image == xperiodic){
		if(xcoord > Pair->xcoord){
			f[3 ] = Pair->f[3 ];
			f[12] = Pair->f[12];
			f[17] = Pair->f[17];
			f[7 ] = Pair->f[7 ];
			f[6 ] = Pair->f[6 ];
		}
		else{
			f[1 ] = Pair->f[1 ];
			f[5 ] = Pair->f[5 ];
			f[8 ] = Pair->f[8 ];
			f[10] = Pair->f[10];
			f[15] = Pair->f[15];
		}
	}
	if (image == yperiodic){
		if(ycoord > Pair->ycoord){
			f[4 ] = Pair->f[4 ];
			f[7 ] = Pair->f[7 ];
			f[8 ] = Pair->f[8 ];
			f[13] = Pair->f[13];
			f[18] = Pair->f[18];
		}
		else{
			f[2 ] = Pair->f[2 ];
			f[5 ] = Pair->f[5 ];
			f[6 ] = Pair->f[6 ];
			f[11] = Pair->f[11];
			f[16] = Pair->f[16];
		}
	}
	if (image == zperiodic){
		if(zcoord > Pair->zcoord){
			f[14] = Pair->f[14];
			f[15] = Pair->f[15];
			f[16] = Pair->f[16];
			f[17] = Pair->f[17];
			f[18] = Pair->f[18];
		}
		else{
			f[9 ] = Pair->f[9 ];
			f[10] = Pair->f[10];
			f[11] = Pair->f[11];
			f[12] = Pair->f[12];
			f[13] = Pair->f[13];
		}
	}

}




void Node::CurvedB()
{
    if (image >= CB1 && image <= CB2 && edge == 0)
    {
		vector<double > temp(19,0);
		for(int i = 1; i<19; i++)
			temp[i] = f[i];
		for(int k = 1; k<19; k++){
			if(nb[k]->image == fluid){
				int m = opposite[k];

//				double uf1 = nb[k]->u0*v[m][0];
//				double uf2 = nb[k]->u1*v[m][1];
//				double uf3 = nb[k]->u2*v[m][2];
//				double uf = uf1+uf2+uf3;
//				double ubf1,ubf3,ubf2,ubf,si;
				if((delta[k] - 0.5)<-Tolerance){
//					f[k] = nb[k]->nb[k]->f[m]*(1.0-2.0*delta[k])+nb[k]->f[m]*2.0*delta[k];
					f[k] = nb[k]->f[m]*(1.0-2.0*delta[k])+temp[m]*2.0*delta[k];
				}
				else{
//					f[k] = nb[k]->f[k]*(2.0*delta[k]-1.0)/(2.0*delta[k])+
//							nb[k]->f[m]/(2.0*delta[k]);
					f[k] = nb[k]->f[m]*(2.0*delta[k]-1.0)/(2.0*delta[k])+
							temp[m]/(2.0*delta[k]);
				}
//				if(D2Q9i == true){
//				f[k] = (1.0-si)*nb[k]->f[m]+si*t[m]*
//						(nb[k]->rho+3.0*ubf+4.5*uf*uf-1.5*nb[k]->usqr);
//				}
//				else{
//				f[k] = (1.0-si)*nb[k]->f[m]+si*t[m]*nb[k]->rho*
//						(1.0+3.0*ubf+4.5*uf*uf-1.5*nb[k]->usqr);
//				}
			}
		}
    }
}

void Node::BBForce(double& ForceX,double& ForceY,double& ForceZ,int BBimage)
{
	if (image == BBimage){
		for(int i = 1; i<19; i++){
			if(nb[i] != NULL){
			//if(nb[i]->image == fluid){
			if(nb[i]->image < 1){
//			ForceX = ForceX+v[opposite[i]][0]*(f[i]+f[i]);
//			ForceY = ForceY+v[opposite[i]][1]*(f[i]+f[i]);
//			ForceZ = ForceZ+v[opposite[i]][2]*(f[i]+f[i]);

			ForceX = ForceX+v[opposite[i]][0]*(f[i]+nb[i]->f[opposite[i]]);
			ForceY = ForceY+v[opposite[i]][1]*(f[i]+nb[i]->f[opposite[i]]);
			ForceZ = ForceZ+v[opposite[i]][2]*(f[i]+nb[i]->f[opposite[i]]);
//			cout<<xcoord<<", "<<ycoord<<", "<<i<<"\n";
			}
			}
		}
	}
}

void Node::DeleteNode()
{
    if (child[0] == NULL){
//		for(int i = 0; i<4; i++)
//	    	delete child[i];
//			//DeleteNode(n->child[i]);
		return;
   		//delete n;
		//return;
	}
    else if (child[0] != NULL)
    {
	for(int i = 0; i<8; i++){
		child[i]->DeleteNode();
		delete child[i];}
	return;
    }
}

void Node::Stream1()
{
	for(int k = 0; k<9; k++){
		if(nb[halfstream[k]] != NULL){
			swap(f[halfstream[k]],nb[halfstream[k]]->f[opposite[halfstream[k]]]);
		}
//		else
//			edge = 1;//flag for edge nodes
	}
}
void Node::Stream2()
{
	for(int k = 0; k<9; k++){
		swap(f[halfstream[k]],f[opposite[halfstream[k]]]);
	}
}

void Node::MeqToFeq()
{

feq[ 0] = 0.052631579*meq[0]+  -0.012531328*meq[1]+  0.047619048*meq[2]                                                                                    ;
feq[ 1] = 0.052631579*meq[0]+ -0.0045948204*meq[1]+ -0.015873016*meq[2]+  0.1*meq[3]+   -0.1*meq[4]                                                        ;
feq[ 2] = 0.052631579*meq[0]+ -0.0045948204*meq[1]+ -0.015873016*meq[2]                            +  0.1*meq[5]+   -0.1*meq[6]                            ;
feq[ 3] = 0.052631579*meq[0]+ -0.0045948204*meq[1]+ -0.015873016*meq[2]+ -0.1*meq[3]+    0.1*meq[4]                                                        ;
feq[ 4] = 0.052631579*meq[0]+ -0.0045948204*meq[1]+ -0.015873016*meq[2]                            + -0.1*meq[5]+    0.1*meq[6]                            ;
feq[ 5] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+  0.1*meq[3]+  0.025*meq[4]+  0.1*meq[5]+  0.025*meq[6]                            ;
feq[ 6] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+ -0.1*meq[3]+ -0.025*meq[4]+  0.1*meq[5]+  0.025*meq[6]                            ;
feq[ 7] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+ -0.1*meq[3]+ -0.025*meq[4]+ -0.1*meq[5]+ -0.025*meq[6]                            ;
feq[ 8] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+  0.1*meq[3]+  0.025*meq[4]+ -0.1*meq[5]+ -0.025*meq[6]                            ;
feq[ 9] = 0.052631579*meq[0]+ -0.0045948204*meq[1]+ -0.015873016*meq[2]                                                        +  0.1*meq[7]+   -0.1*meq[8];
feq[10] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+  0.1*meq[3]+  0.025*meq[4]                            +  0.1*meq[7]+  0.025*meq[8];
feq[11] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]                            +  0.1*meq[5]+  0.025*meq[6]+  0.1*meq[7]+  0.025*meq[8];
feq[12] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+ -0.1*meq[3]+ -0.025*meq[4]                            +  0.1*meq[7]+  0.025*meq[8];
feq[13] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]                            + -0.1*meq[5]+ -0.025*meq[6]+  0.1*meq[7]+  0.025*meq[8];
feq[14] = 0.052631579*meq[0]+ -0.0045948204*meq[1]+ -0.015873016*meq[2]                                                        + -0.1*meq[7]+    0.1*meq[8];
feq[15] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+  0.1*meq[3]+  0.025*meq[4]                            + -0.1*meq[7]+ -0.025*meq[8];
feq[16] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]                            +  0.1*meq[5]+  0.025*meq[6]+ -0.1*meq[7]+ -0.025*meq[8];
feq[17] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+ -0.1*meq[3]+ -0.025*meq[4]                            + -0.1*meq[7]+ -0.025*meq[8];
feq[18] = 0.052631579*meq[0]+  0.0033416876*meq[1]+  0.003968254*meq[2]+    0*meq[3]+      0*meq[4]+ -0.1*meq[5]+ -0.025*meq[6]+ -0.1*meq[7]+ -0.025*meq[8];

//feq[ 0] +=            0*meq[9]+            0*meq[10]+            0*meq[11]+            0*meq[12]+     0*meq[13]+     0*meq[14]+     0*meq[15]+      0*meq[16]+      0*meq[17]+      0*meq[18];
feq[ 1] +=  0.055555556*meq[9]+ -0.055555556*meq[10];//+            0*meq[11]+            0*meq[12]+     0*meq[13]+     0*meq[14]+     0*meq[15]+      0*meq[16]+      0*meq[17]+      0*meq[18];
feq[ 2] += -0.027777778*meq[9]+  0.027777778*meq[10]+  0.083333333*meq[11]+ -0.083333333*meq[12];//+     0*meq[13]+     0*meq[14]+     0*meq[15]+      0*meq[16]+      0*meq[17]+      0*meq[18];
feq[ 3] +=  0.055555556*meq[9]+ -0.055555556*meq[10];//+            0*meq[11]+            0*meq[12]+     0*meq[13]+     0*meq[14]+     0*meq[15]+      0*meq[16]+      0*meq[17]+      0*meq[18];
feq[ 4] += -0.027777778*meq[9]+  0.027777778*meq[10]+  0.083333333*meq[11]+ -0.083333333*meq[12];//+     0*meq[13]+     0*meq[14]+     0*meq[15]+      0*meq[16]+      0*meq[17]+      0*meq[18];
feq[ 5] +=  0.027777778*meq[9]+  0.013888889*meq[10]+  0.083333333*meq[11]+  0.041666667*meq[12]+  0.25*meq[13]                              +  0.125*meq[16]+ -0.125*meq[17]                ;
feq[ 6] +=  0.027777778*meq[9]+  0.013888889*meq[10]+  0.083333333*meq[11]+  0.041666667*meq[12]+ -0.25*meq[13]                              + -0.125*meq[16]+ -0.125*meq[17]                ;
feq[ 7] +=  0.027777778*meq[9]+  0.013888889*meq[10]+  0.083333333*meq[11]+  0.041666667*meq[12]+  0.25*meq[13]                              + -0.125*meq[16]+  0.125*meq[17]                ;
feq[ 8] +=  0.027777778*meq[9]+  0.013888889*meq[10]+  0.083333333*meq[11]+  0.041666667*meq[12]+ -0.25*meq[13]                              +  0.125*meq[16]+  0.125*meq[17]                ;
feq[ 9] += -0.027777778*meq[9]+  0.027777778*meq[10]+ -0.083333333*meq[11]+  0.083333333*meq[12]                                                                                             ;
feq[10] +=  0.027777778*meq[9]+  0.013888889*meq[10]+ -0.083333333*meq[11]+ -0.041666667*meq[12]                              +  0.25*meq[15]+ -0.125*meq[16]                +  0.125*meq[18];
feq[11] += -0.055555556*meq[9]+ -0.027777778*meq[10]                                                           +  0.25*meq[14]                               +  0.125*meq[17]+ -0.125*meq[18];
feq[12] +=  0.027777778*meq[9]+  0.013888889*meq[10]+ -0.083333333*meq[11]+ -0.041666667*meq[12]                              + -0.25*meq[15]+  0.125*meq[16]                +  0.125*meq[18];
feq[13] += -0.055555556*meq[9]+ -0.027777778*meq[10]                                                             -0.25*meq[14]                               + -0.125*meq[17]+ -0.125*meq[18];
feq[14] += -0.027777778*meq[9]+  0.027777778*meq[10]+ -0.083333333*meq[11]+  0.083333333*meq[12]                                                                                             ;
feq[15] +=  0.027777778*meq[9]+  0.013888889*meq[10]+ -0.083333333*meq[11]+ -0.041666667*meq[12]                              + -0.25*meq[15]+ -0.125*meq[16]                + -0.125*meq[18];
feq[16] += -0.055555556*meq[9]+ -0.027777778*meq[10]                                                           + -0.25*meq[14]                               +  0.125*meq[17]+  0.125*meq[18];
feq[17] +=  0.027777778*meq[9]+  0.013888889*meq[10]+ -0.083333333*meq[11]+ -0.041666667*meq[12]                              +  0.25*meq[15]+  0.125*meq[16]                + -0.125*meq[18];
feq[18] += -0.055555556*meq[9]+ -0.027777778*meq[10]                                                           +  0.25*meq[14]                               + -0.125*meq[17]+  0.125*meq[18];

feqflag = 1;



}



void Node::SpatialInterp1()
{
/*
This is a parent level function. Use on parent node to interpolate values for its children nodes
*/
	double SF = omega[level]*(1.0-omega[level+1])/((1.0-omega[level])*omega[level+1]*2.0);
//	double SF = -(3.0*nu+0.5*(1.0/(level+1.0)))/(3.0*nu+0.5*(1.0/(level+2.0)));//omega[level]/(omega[level+1]*2.0);
//	double SF = omega[level]/(omega[level+1]*2.0);
	double fcoarse;
	double feqcoarse;
	if(child[0] == NULL){
		cout<<"error in interpolation\n";
		exit;}
	else if(MRT == true){ //MRT
		if(child[0]->edge == 1){
		for(int k = 0; k<19; k++){
			if(feqflag == 0) MeqToFeq();
			if(nb[ 3]->feqflag == 0) nb[ 3]->MeqToFeq();//should put this outside of the k loop
			if(nb[ 4]->feqflag == 0) nb[ 4]->MeqToFeq();
			if(nb[14]->feqflag == 0) nb[14]->MeqToFeq();
			if(nb[ 7]->feqflag == 0) nb[ 7]->MeqToFeq();
			if(nb[17]->feqflag == 0) nb[17]->MeqToFeq();
			if(nb[18]->feqflag == 0) nb[18]->MeqToFeq();
			if(nb[14]->nb[7]->feqflag == 0) nb[14]->nb[7]->MeqToFeq();

			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[4]->f[k]+nb[14]->f[k])+
								3.0*(nb[7]->f[k]+nb[17]->f[k]+nb[18]->f[k])+
								nb[14]->nb[7]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[4]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[7]->feq[k]+nb[17]->feq[k]+nb[18]->feq[k])+
								nb[14]->nb[7]->feq[k]);
			child[0]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[0]->feq[k] = feqcoarse;
		}
		}
		if(child[1]->edge == 1){
		for(int k = 0; k<19; k++){
			if(feqflag == 0) MeqToFeq();
			if(nb[ 1]->feqflag == 0) nb[ 1]->MeqToFeq();
			if(nb[ 4]->feqflag == 0) nb[ 4]->MeqToFeq();
			if(nb[14]->feqflag == 0) nb[14]->MeqToFeq();
			if(nb[ 8]->feqflag == 0) nb[ 8]->MeqToFeq();
			if(nb[15]->feqflag == 0) nb[15]->MeqToFeq();
			if(nb[18]->feqflag == 0) nb[18]->MeqToFeq();
			if(nb[14]->nb[8]->feqflag == 0) nb[14]->nb[8]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[4]->f[k]+nb[14]->f[k])+
								3.0*(nb[8]->f[k]+nb[18]->f[k]+nb[15]->f[k])+
								nb[14]->nb[8]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[4]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[8]->feq[k]+nb[18]->feq[k]+nb[15]->feq[k])+
								nb[14]->nb[8]->feq[k]);
			child[1]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[1]->feq[k] = feqcoarse;
		}
		}
		if(child[2]->edge == 1){
		for(int k = 0; k<19; k++){
			if(feqflag == 0) MeqToFeq();
			if(nb[ 1]->feqflag == 0) nb[ 1]->MeqToFeq();
			if(nb[ 2]->feqflag == 0) nb[ 2]->MeqToFeq();
			if(nb[14]->feqflag == 0) nb[14]->MeqToFeq();
			if(nb[ 5]->feqflag == 0) nb[ 5]->MeqToFeq();
			if(nb[15]->feqflag == 0) nb[15]->MeqToFeq();
			if(nb[16]->feqflag == 0) nb[16]->MeqToFeq();
			if(nb[14]->nb[5]->feqflag == 0) nb[14]->nb[5]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[2]->f[k]+nb[14]->f[k])+
								3.0*(nb[5]->f[k]+nb[16]->f[k]+nb[15]->f[k])+
								nb[14]->nb[5]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[2]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[5]->feq[k]+nb[16]->feq[k]+nb[15]->feq[k])+
								nb[14]->nb[5]->feq[k]);
			child[2]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[2]->feq[k] = feqcoarse;
		}
		}
		if(child[3]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 3]->feqflag == 0) nb[ 3]->MeqToFeq();
			if(nb[ 2]->feqflag == 0) nb[ 2]->MeqToFeq();
			if(nb[14]->feqflag == 0) nb[14]->MeqToFeq();
			if(nb[ 6]->feqflag == 0) nb[ 6]->MeqToFeq();
			if(nb[17]->feqflag == 0) nb[17]->MeqToFeq();
			if(nb[16]->feqflag == 0) nb[16]->MeqToFeq();
			if(nb[14]->nb[6]->feqflag == 0) nb[14]->nb[6]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[2]->f[k]+nb[14]->f[k])+
								3.0*(nb[6]->f[k]+nb[16]->f[k]+nb[17]->f[k])+
								nb[14]->nb[6]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[2]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[6]->feq[k]+nb[16]->feq[k]+nb[17]->feq[k])+
								nb[14]->nb[6]->feq[k]);
			child[3]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[3]->feq[k] = feqcoarse;
		}
		}


		if(child[4]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 3]->feqflag == 0) nb[ 3]->MeqToFeq();
			if(nb[ 4]->feqflag == 0) nb[ 4]->MeqToFeq();
			if(nb[ 9]->feqflag == 0) nb[ 9]->MeqToFeq();
			if(nb[ 7]->feqflag == 0) nb[ 7]->MeqToFeq();
			if(nb[12]->feqflag == 0) nb[12]->MeqToFeq();
			if(nb[13]->feqflag == 0) nb[13]->MeqToFeq();
			if(nb[9]->nb[7]->feqflag == 0) nb[9]->nb[7]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[4]->f[k]+nb[9]->f[k])+
								3.0*(nb[7]->f[k]+nb[12]->f[k]+nb[13]->f[k])+
								nb[9]->nb[7]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[4]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[7]->feq[k]+nb[12]->feq[k]+nb[13]->feq[k])+
								nb[9]->nb[7]->feq[k]);
			child[4]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[4]->feq[k] = feqcoarse;
		}
		}
		if(child[5]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 1]->feqflag == 0) nb[ 1]->MeqToFeq();
			if(nb[ 4]->feqflag == 0) nb[ 4]->MeqToFeq();
			if(nb[ 9]->feqflag == 0) nb[ 9]->MeqToFeq();
			if(nb[ 8]->feqflag == 0) nb[ 8]->MeqToFeq();
			if(nb[10]->feqflag == 0) nb[10]->MeqToFeq();
			if(nb[13]->feqflag == 0) nb[13]->MeqToFeq();
			if(nb[9]->nb[8]->feqflag == 0) nb[9]->nb[8]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[4]->f[k]+nb[9]->f[k])+
								3.0*(nb[8]->f[k]+nb[13]->f[k]+nb[10]->f[k])+
								nb[9]->nb[8]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[4]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[8]->feq[k]+nb[13]->feq[k]+nb[10]->feq[k])+
								nb[9]->nb[8]->feq[k]);
			child[5]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[5]->feq[k] = feqcoarse;
		}
		}
		if(child[6]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 1]->feqflag == 0) nb[ 1]->MeqToFeq();
			if(nb[ 2]->feqflag == 0) nb[ 2]->MeqToFeq();
			if(nb[ 9]->feqflag == 0) nb[ 9]->MeqToFeq();
			if(nb[ 5]->feqflag == 0) nb[ 5]->MeqToFeq();
			if(nb[10]->feqflag == 0) nb[10]->MeqToFeq();
			if(nb[11]->feqflag == 0) nb[11]->MeqToFeq();
			if(nb[9]->nb[8]->feqflag == 0) nb[9]->nb[8]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[2]->f[k]+nb[9]->f[k])+
								3.0*(nb[5]->f[k]+nb[11]->f[k]+nb[10]->f[k])+
								nb[9]->nb[5]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[2]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[5]->feq[k]+nb[11]->feq[k]+nb[10]->feq[k])+
								nb[9]->nb[5]->feq[k]);
			child[6]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[6]->feq[k] = feqcoarse;
		}
		}
		if(child[7]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 3]->feqflag == 0) nb[ 3]->MeqToFeq();
			if(nb[ 2]->feqflag == 0) nb[ 2]->MeqToFeq();
			if(nb[ 9]->feqflag == 0) nb[ 9]->MeqToFeq();
			if(nb[ 6]->feqflag == 0) nb[ 6]->MeqToFeq();
			if(nb[12]->feqflag == 0) nb[12]->MeqToFeq();
			if(nb[11]->feqflag == 0) nb[11]->MeqToFeq();
			if(nb[9]->nb[6]->feqflag == 0) nb[9]->nb[6]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[2]->f[k]+nb[9]->f[k])+
								3.0*(nb[6]->f[k]+nb[11]->f[k]+nb[12]->f[k])+
								nb[9]->nb[6]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[2]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[6]->feq[k]+nb[11]->feq[k]+nb[12]->feq[k])+
								nb[9]->nb[6]->feq[k]);
			child[7]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[7]->feq[k] = feqcoarse;
		}
		}
	}//end MRT
	else{ //BGK
		if(child[0]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[4]->f[k]+nb[14]->f[k])+
								3.0*(nb[7]->f[k]+nb[17]->f[k]+nb[18]->f[k])+
								nb[14]->nb[7]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[4]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[7]->feq[k]+nb[17]->feq[k]+nb[18]->feq[k])+
								nb[14]->nb[7]->feq[k]);
			child[0]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[0]->feq[k] = feqcoarse;
		}
		}
		if(child[1]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[4]->f[k]+nb[14]->f[k])+
								3.0*(nb[8]->f[k]+nb[18]->f[k]+nb[15]->f[k])+
								nb[14]->nb[8]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[4]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[8]->feq[k]+nb[18]->feq[k]+nb[15]->feq[k])+
								nb[14]->nb[8]->feq[k]);
			child[1]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[1]->feq[k] = feqcoarse;
		}
		}
		if(child[2]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[2]->f[k]+nb[14]->f[k])+
								3.0*(nb[5]->f[k]+nb[16]->f[k]+nb[15]->f[k])+
								nb[14]->nb[5]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[2]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[5]->feq[k]+nb[16]->feq[k]+nb[15]->feq[k])+
								nb[14]->nb[5]->feq[k]);
			child[2]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[2]->feq[k] = feqcoarse;
		}
		}
		if(child[3]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[2]->f[k]+nb[14]->f[k])+
								3.0*(nb[6]->f[k]+nb[16]->f[k]+nb[17]->f[k])+
								nb[14]->nb[6]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[2]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[6]->feq[k]+nb[16]->feq[k]+nb[17]->feq[k])+
								nb[14]->nb[6]->feq[k]);
			child[3]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[3]->feq[k] = feqcoarse;
		}
		}


		if(child[4]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[4]->f[k]+nb[9]->f[k])+
								3.0*(nb[7]->f[k]+nb[12]->f[k]+nb[13]->f[k])+
								nb[9]->nb[7]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[4]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[7]->feq[k]+nb[12]->feq[k]+nb[13]->feq[k])+
								nb[9]->nb[7]->feq[k]);
			child[4]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[4]->feq[k] = feqcoarse;
		}
		}
		if(child[5]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[4]->f[k]+nb[9]->f[k])+
								3.0*(nb[8]->f[k]+nb[13]->f[k]+nb[10]->f[k])+
								nb[9]->nb[8]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[4]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[8]->feq[k]+nb[13]->feq[k]+nb[10]->feq[k])+
								nb[9]->nb[8]->feq[k]);
			child[5]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[5]->feq[k] = feqcoarse;
		}
		}
		if(child[6]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[2]->f[k]+nb[9]->f[k])+
								3.0*(nb[5]->f[k]+nb[11]->f[k]+nb[10]->f[k])+
								nb[9]->nb[5]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[2]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[5]->feq[k]+nb[11]->feq[k]+nb[10]->feq[k])+
								nb[9]->nb[5]->feq[k]);
			child[6]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[6]->feq[k] = feqcoarse;
		}
		}
		if(child[7]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[2]->f[k]+nb[9]->f[k])+
								3.0*(nb[6]->f[k]+nb[11]->f[k]+nb[12]->f[k])+
								nb[9]->nb[6]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[2]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[6]->feq[k]+nb[11]->feq[k]+nb[12]->feq[k])+
								nb[9]->nb[6]->feq[k]);
			child[7]->f[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
			child[7]->feq[k] = feqcoarse;
		}
		}
	}
}
void Node::SpatialInterp2()
{
/*
This function is identical to SpatialInterp1, except the interpolated f values are stored in ftemp
*/
	double SF = omega[level]*(1.0-omega[level+1])/((1.0-omega[level])*omega[level+1]*2.0);
//	double SF = -(3.0*nu+0.5*(1.0/(level+1.0)))/(3.0*nu+0.5*(1.0/(level+2.0)));//omega[level]/(omega[level+1]*2.0);
//	double SF = omega[level]/(omega[level+1]*2.0);
	double fcoarse;
	double feqcoarse;
	if(child[0] == NULL){
		cout<<"error in interpolation\n";
		exit;}
	else if(MRT == true){ //MRT
		if(child[0]->edge == 1){
		for(int k = 0; k<19; k++){
			if(feqflag == 0) MeqToFeq();
			if(nb[ 3]->feqflag == 0) nb[ 3]->MeqToFeq();
			if(nb[ 4]->feqflag == 0) nb[ 4]->MeqToFeq();
			if(nb[14]->feqflag == 0) nb[14]->MeqToFeq();
			if(nb[ 7]->feqflag == 0) nb[ 7]->MeqToFeq();
			if(nb[17]->feqflag == 0) nb[17]->MeqToFeq();
			if(nb[18]->feqflag == 0) nb[18]->MeqToFeq();
			if(nb[14]->nb[7]->feqflag == 0) nb[14]->nb[7]->MeqToFeq();

			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[4]->f[k]+nb[14]->f[k])+
								3.0*(nb[7]->f[k]+nb[17]->f[k]+nb[18]->f[k])+
								nb[14]->nb[7]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[4]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[7]->feq[k]+nb[17]->feq[k]+nb[18]->feq[k])+
								nb[14]->nb[7]->feq[k]);
			child[0]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[1]->edge == 1){
		for(int k = 0; k<19; k++){
			if(feqflag == 0) MeqToFeq();
			if(nb[ 1]->feqflag == 0) nb[ 1]->MeqToFeq();
			if(nb[ 4]->feqflag == 0) nb[ 4]->MeqToFeq();
			if(nb[14]->feqflag == 0) nb[14]->MeqToFeq();
			if(nb[ 8]->feqflag == 0) nb[ 8]->MeqToFeq();
			if(nb[15]->feqflag == 0) nb[15]->MeqToFeq();
			if(nb[18]->feqflag == 0) nb[18]->MeqToFeq();
			if(nb[14]->nb[8]->feqflag == 0) nb[14]->nb[8]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[4]->f[k]+nb[14]->f[k])+
								3.0*(nb[8]->f[k]+nb[18]->f[k]+nb[15]->f[k])+
								nb[14]->nb[8]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[4]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[8]->feq[k]+nb[18]->feq[k]+nb[15]->feq[k])+
								nb[14]->nb[8]->feq[k]);
			child[1]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[2]->edge == 1){
		for(int k = 0; k<19; k++){
			if(feqflag == 0) MeqToFeq();
			if(nb[ 1]->feqflag == 0) nb[ 1]->MeqToFeq();
			if(nb[ 2]->feqflag == 0) nb[ 2]->MeqToFeq();
			if(nb[14]->feqflag == 0) nb[14]->MeqToFeq();
			if(nb[ 5]->feqflag == 0) nb[ 5]->MeqToFeq();
			if(nb[15]->feqflag == 0) nb[15]->MeqToFeq();
			if(nb[16]->feqflag == 0) nb[16]->MeqToFeq();
			if(nb[14]->nb[5]->feqflag == 0) nb[14]->nb[5]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[2]->f[k]+nb[14]->f[k])+
								3.0*(nb[5]->f[k]+nb[16]->f[k]+nb[15]->f[k])+
								nb[14]->nb[5]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[2]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[5]->feq[k]+nb[16]->feq[k]+nb[15]->feq[k])+
								nb[14]->nb[5]->feq[k]);
			child[2]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[3]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 3]->feqflag == 0) nb[ 3]->MeqToFeq();
			if(nb[ 2]->feqflag == 0) nb[ 2]->MeqToFeq();
			if(nb[14]->feqflag == 0) nb[14]->MeqToFeq();
			if(nb[ 6]->feqflag == 0) nb[ 6]->MeqToFeq();
			if(nb[17]->feqflag == 0) nb[17]->MeqToFeq();
			if(nb[16]->feqflag == 0) nb[16]->MeqToFeq();
			if(nb[14]->nb[6]->feqflag == 0) nb[14]->nb[6]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[2]->f[k]+nb[14]->f[k])+
								3.0*(nb[6]->f[k]+nb[16]->f[k]+nb[17]->f[k])+
								nb[14]->nb[6]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[2]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[6]->feq[k]+nb[16]->feq[k]+nb[17]->feq[k])+
								nb[14]->nb[6]->feq[k]);
			child[3]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}


		if(child[4]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 3]->feqflag == 0) nb[ 3]->MeqToFeq();
			if(nb[ 4]->feqflag == 0) nb[ 4]->MeqToFeq();
			if(nb[ 9]->feqflag == 0) nb[ 9]->MeqToFeq();
			if(nb[ 7]->feqflag == 0) nb[ 7]->MeqToFeq();
			if(nb[12]->feqflag == 0) nb[12]->MeqToFeq();
			if(nb[13]->feqflag == 0) nb[13]->MeqToFeq();
			if(nb[9]->nb[7]->feqflag == 0) nb[9]->nb[7]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[4]->f[k]+nb[9]->f[k])+
								3.0*(nb[7]->f[k]+nb[12]->f[k]+nb[13]->f[k])+
								nb[9]->nb[7]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[4]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[7]->feq[k]+nb[12]->feq[k]+nb[13]->feq[k])+
								nb[9]->nb[7]->feq[k]);
			child[4]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[5]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 1]->feqflag == 0) nb[ 1]->MeqToFeq();
			if(nb[ 4]->feqflag == 0) nb[ 4]->MeqToFeq();
			if(nb[ 9]->feqflag == 0) nb[ 9]->MeqToFeq();
			if(nb[ 8]->feqflag == 0) nb[ 8]->MeqToFeq();
			if(nb[10]->feqflag == 0) nb[10]->MeqToFeq();
			if(nb[13]->feqflag == 0) nb[13]->MeqToFeq();
			if(nb[9]->nb[8]->feqflag == 0) nb[9]->nb[8]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[4]->f[k]+nb[9]->f[k])+
								3.0*(nb[8]->f[k]+nb[13]->f[k]+nb[10]->f[k])+
								nb[9]->nb[8]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[4]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[8]->feq[k]+nb[13]->feq[k]+nb[10]->feq[k])+
								nb[9]->nb[8]->feq[k]);
			child[5]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[6]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 1]->feqflag == 0) nb[ 1]->MeqToFeq();
			if(nb[ 2]->feqflag == 0) nb[ 2]->MeqToFeq();
			if(nb[ 9]->feqflag == 0) nb[ 9]->MeqToFeq();
			if(nb[ 5]->feqflag == 0) nb[ 5]->MeqToFeq();
			if(nb[10]->feqflag == 0) nb[10]->MeqToFeq();
			if(nb[11]->feqflag == 0) nb[11]->MeqToFeq();
			if(nb[9]->nb[8]->feqflag == 0) nb[9]->nb[8]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[2]->f[k]+nb[9]->f[k])+
								3.0*(nb[5]->f[k]+nb[11]->f[k]+nb[10]->f[k])+
								nb[9]->nb[5]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[2]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[5]->feq[k]+nb[11]->feq[k]+nb[10]->feq[k])+
								nb[9]->nb[5]->feq[k]);
			child[6]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[7]->edge == 1){
		for(int k = 0; k<19; k++){
			if(nb[ 3]->feqflag == 0) nb[ 3]->MeqToFeq();
			if(nb[ 2]->feqflag == 0) nb[ 2]->MeqToFeq();
			if(nb[ 9]->feqflag == 0) nb[ 9]->MeqToFeq();
			if(nb[ 6]->feqflag == 0) nb[ 6]->MeqToFeq();
			if(nb[12]->feqflag == 0) nb[12]->MeqToFeq();
			if(nb[11]->feqflag == 0) nb[11]->MeqToFeq();
			if(nb[9]->nb[6]->feqflag == 0) nb[9]->nb[6]->MeqToFeq();
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[2]->f[k]+nb[9]->f[k])+
								3.0*(nb[6]->f[k]+nb[11]->f[k]+nb[12]->f[k])+
								nb[9]->nb[6]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[2]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[6]->feq[k]+nb[11]->feq[k]+nb[12]->feq[k])+
								nb[9]->nb[6]->feq[k]);
			child[7]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
	}//end of MRT

	else{ //BGK
		if(child[0]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[4]->f[k]+nb[14]->f[k])+
								3.0*(nb[7]->f[k]+nb[17]->f[k]+nb[18]->f[k])+
								nb[14]->nb[7]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[4]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[7]->feq[k]+nb[17]->feq[k]+nb[18]->feq[k])+
								nb[14]->nb[7]->feq[k]);
			child[0]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[1]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[4]->f[k]+nb[14]->f[k])+
								3.0*(nb[8]->f[k]+nb[18]->f[k]+nb[15]->f[k])+
								nb[14]->nb[8]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[4]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[8]->feq[k]+nb[18]->feq[k]+nb[15]->feq[k])+
								nb[14]->nb[8]->feq[k]);
			child[1]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[2]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[2]->f[k]+nb[14]->f[k])+
								3.0*(nb[5]->f[k]+nb[16]->f[k]+nb[15]->f[k])+
								nb[14]->nb[5]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[2]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[5]->feq[k]+nb[16]->feq[k]+nb[15]->feq[k])+
								nb[14]->nb[5]->feq[k]);
			child[2]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[3]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[2]->f[k]+nb[14]->f[k])+
								3.0*(nb[6]->f[k]+nb[16]->f[k]+nb[17]->f[k])+
								nb[14]->nb[6]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[2]->feq[k]+nb[14]->feq[k])+
								3.0*(nb[6]->feq[k]+nb[16]->feq[k]+nb[17]->feq[k])+
								nb[14]->nb[6]->feq[k]);
			child[3]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}


		if(child[4]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[4]->f[k]+nb[9]->f[k])+
								3.0*(nb[7]->f[k]+nb[12]->f[k]+nb[13]->f[k])+
								nb[9]->nb[7]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[4]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[7]->feq[k]+nb[12]->feq[k]+nb[13]->feq[k])+
								nb[9]->nb[7]->feq[k]);
			child[4]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[5]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[4]->f[k]+nb[9]->f[k])+
								3.0*(nb[8]->f[k]+nb[13]->f[k]+nb[10]->f[k])+
								nb[9]->nb[8]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[4]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[8]->feq[k]+nb[13]->feq[k]+nb[10]->feq[k])+
								nb[9]->nb[8]->feq[k]);
			child[5]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[6]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[1]->f[k]+nb[2]->f[k]+nb[9]->f[k])+
								3.0*(nb[5]->f[k]+nb[11]->f[k]+nb[10]->f[k])+
								nb[9]->nb[5]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[1]->feq[k]+nb[2]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[5]->feq[k]+nb[11]->feq[k]+nb[10]->feq[k])+
								nb[9]->nb[5]->feq[k]);
			child[6]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
		if(child[7]->edge == 1){
		for(int k = 0; k<19; k++){
			fcoarse = 0.015625*(27.0*f[k]+9.0*(nb[3]->f[k]+nb[2]->f[k]+nb[9]->f[k])+
								3.0*(nb[6]->f[k]+nb[11]->f[k]+nb[12]->f[k])+
								nb[9]->nb[6]->f[k]);
			feqcoarse = 0.015625*(27.0*feq[k]+9.0*(nb[3]->feq[k]+nb[2]->feq[k]+nb[9]->feq[k])+
								3.0*(nb[6]->feq[k]+nb[11]->feq[k]+nb[12]->feq[k])+
								nb[9]->nb[6]->feq[k]);
			child[7]->ftemp[k] = feqcoarse*(1.0-SF)+SF*fcoarse;
		}
		}
	}
}
void Node::SpatialAverage()
{
/*
*/
	double SF = (1.0-omega[level])*omega[level+1]*2.0/(omega[level]*(1.0-omega[level+1]));
//	double SF = -(3.0*nu+0.5*(1.0/(level+2.0)))/(3.0*nu+0.5*(1.0/(level+1.0)));//omega[level]/(omega[level+1]*2.0);
//	double SF = omega[level+1]*2.0/(omega[level]);
	double ffine;
	double feqfine;
	if(child[0] == NULL){
		cout<<"error in averaging\n";
		exit;}
	else if(child[0]->edge+child[1]->edge+child[2]->edge+child[3]->edge+
			child[4]->edge+child[5]->edge+child[6]->edge+child[7]->edge == 0 &&
				child[0]->image == fluid){ 
		if(MRT == true){
			for(int i = 0; i<8; i++){
				child[i]->MeqToFeq();
			}
			for(int k = 0; k<19; k++){
				ffine = 0.125*(child[0]->f[k]+child[1]->f[k]+child[2]->f[k]+child[3]->f[k]+
								child[4]->f[k]+child[5]->f[k]+child[6]->f[k]+child[7]->f[k]);
				feqfine = 0.125*(child[0]->feq[k]+child[1]->feq[k]+child[2]->feq[k]+child[3]->feq[k]+
								child[4]->feq[k]+child[5]->feq[k]+child[6]->feq[k]+child[7]->feq[k]);
				f[k] = feqfine*(1.0-SF)+ffine*SF;
			}
		}
		else{//BGK
		for(int k = 0; k<19; k++){
			ffine = 0.125*(child[0]->f[k]+child[1]->f[k]+child[2]->f[k]+child[3]->f[k]+
							child[4]->f[k]+child[5]->f[k]+child[6]->f[k]+child[7]->f[k]);
			feqfine = 0.125*(child[0]->feq[k]+child[1]->feq[k]+child[2]->feq[k]+child[3]->feq[k]+
							child[4]->feq[k]+child[5]->feq[k]+child[6]->feq[k]+child[7]->feq[k]);
			f[k] = feqfine*(1.0-SF)+ffine*SF;
		}
		}

	}
}

void Node::TemporalInterp()
{
	for(int k=0;k<8;k++){
		if(child[k] == NULL){
			cout<<"something wrong in TemporalInterp\n";
			exit;}
		else if(child[k]->edge == 1){
			for(int i=0;i<19;i++){
				child[k]->ftemp[i] = 0.5*(child[k]->f[i]+child[k]->ftemp[i]);
			}
		}
	}
}
void Node::TemporalInterp2()
{
	for(int k=0;k<8;k++){
		if(child[k] == NULL){
			cout<<"something wrong in TemporalInterp\n";
			exit;}
		else if(child[k]->edge == 1){
			for(int i=0;i<19;i++)
				child[k]->f[i] = child[k]->ftemp[i];
		}
	}
}

void Node::ComputeDelta(double obX, double obY, double obX2, double obY2, double R)
/*
this method is only for cylinders in the xz plane. obY refers to z location of cylinder
*/
{
	if(image == CB1){
		for(int i = 1; i<19; i++){
			if(nb[i]->image == fluid){
				if(i == 1 || i == 5 || i == 8)
					delta[i] = (sqrt(R*R-(zcoord-obY)*(zcoord-obY))-(xcoord-obX))/(pow(2.0,-level));
				else if(i == 3 || i == 6 || i == 7)
					delta[i] = (sqrt(R*R-(zcoord-obY)*(zcoord-obY))-(obX-xcoord))/(pow(2.0,-level));
				else if(i == 9 || i == 11 || i == 13)
					delta[i] = (sqrt(R*R-(xcoord-obX)*(xcoord-obX))-(zcoord-obY))/(pow(2.0,-level));
				else if(i == 14 || i == 16 || i == 18)
					delta[i] = (sqrt(R*R-(xcoord-obX)*(xcoord-obX))-(obY-zcoord))/(pow(2.0,-level));
				else{//coord transformation
					double x2 =  (xcoord-obX)*sqrt(2.0)*0.5+(zcoord-obY)*sqrt(2.0)*0.5;
					double y2 = -(xcoord-obX)*sqrt(2.0)*0.5+(zcoord-obY)*sqrt(2.0)*0.5;
					if(i == 10)
						delta[i] = (sqrt(R*R-y2*y2)-x2)/(sqrt(2.0)*pow(2.0,-level));
					else if(i == 17)
						delta[i] = (sqrt(R*R-y2*y2)+x2)/(sqrt(2.0)*pow(2.0,-level));
					else if(i == 12)
						delta[i] = (sqrt(R*R-x2*x2)-y2)/(sqrt(2.0)*pow(2.0,-level));
					else if(i == 15)
						delta[i] = (sqrt(R*R-x2*x2)+y2)/(sqrt(2.0)*pow(2.0,-level));
				}
			}
		}
	}
	else if(image == CB2){
		for(int i = 1; i<19; i++){
			if(nb[i]->image == fluid){
				if(i == 1 || i == 5 || i == 8)
					delta[i] = (sqrt(R*R-(zcoord-obY2)*(zcoord-obY2))-(xcoord-obX2))/(pow(2.0,-level));
				else if(i == 3 || i == 6 || i == 7)
					delta[i] = (sqrt(R*R-(zcoord-obY2)*(zcoord-obY2))-(obX2-xcoord))/(pow(2.0,-level));
				else if(i == 9 || i == 11 || i == 13)
					delta[i] = (sqrt(R*R-(xcoord-obX2)*(xcoord-obX2))-(zcoord-obY2))/(pow(2.0,-level));
				else if(i == 14 || i == 16 || i == 18)
					delta[i] = (sqrt(R*R-(xcoord-obX2)*(xcoord-obX2))-(obY2-zcoord))/(pow(2.0,-level));
				else{//coord transformation
					double x2 =  (xcoord-obX2)*sqrt(2.0)*0.5+(zcoord-obY2)*sqrt(2.0)*0.5;
					double y2 = -(xcoord-obX2)*sqrt(2.0)*0.5+(zcoord-obY2)*sqrt(2.0)*0.5;
					if(i == 10)
						delta[i] = (sqrt(R*R-y2*y2)-x2)/(sqrt(2.0)*pow(2.0,-level));
					else if(i == 17)
						delta[i] = (sqrt(R*R-y2*y2)+x2)/(sqrt(2.0)*pow(2.0,-level));
					else if(i == 12)
						delta[i] = (sqrt(R*R-x2*x2)-y2)/(sqrt(2.0)*pow(2.0,-level));
					else if(i == 15)
						delta[i] = (sqrt(R*R-x2*x2)+y2)/(sqrt(2.0)*pow(2.0,-level));
				}
			}
		}
	}
}


