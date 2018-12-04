#pragma once

/*

	Definitions for other include files

*/

//General

/*!
	Directions
*/

enum Direction{
	Dir_X,
	Dir_Y,
	Dir_Z
};


//Grid variables

/*!
	Specifies a vector position
*/

struct GridPosition{

public:

	long long x, y, z;

private:

	void Set(long long X, long long Y, long long Z){
		x = X;
		y = Y;
		z = Z;
	}

public:

	GridPosition(long long X, long long Y, long long Z){
		Set(X, Y, Z);
	}

	GridPosition(int X, int Y, int Z){
		Set(X, Y, Z);
	}

	GridPosition(long long X, long long Y){
		Set(X, Y, 0);
	}

	GridPosition(int X, int Y){
		Set(X, Y, 0);
	}

	GridPosition(long long X){
		Set(X, 0, 0);
	}

	GridPosition(int X){
		Set(X, 0, 0);
	}

	long long operator[](int i){
	
		if(i==0) return x;
		if(i==1) return y;
		if(i==2) return z;

		return 0;
	}

};

/*!
	Specifies a 3D region
*/



struct GridRegion{

public:

	long long x0, x1;
	long long y0, y1;
	long long z0, z1;

private:

	static char OutStr[256];	//Static char for writing region to string

	void Set(long long XStart, long long XEnd, long long YStart, long long YEnd, long long ZStart, long long ZEnd){

		x0 = XStart;
		x1 = XEnd;
		y0 = YStart;
		y1 = YEnd;
		z0 = ZStart;
		z1 = ZEnd;

	}

public:

	GridRegion(long long XStart, long long XEnd, long long YStart, long long YEnd, long long ZStart, long long ZEnd){
		Set(XStart, XEnd, YStart, YEnd, ZStart, ZEnd);
	}

	GridRegion(int XStart, int XEnd, int YStart, int YEnd, int ZStart, int ZEnd){
		Set(XStart, XEnd, YStart, YEnd, ZStart, ZEnd);
	}

	GridRegion(long long XStart, long long XEnd, long long YStart, long long YEnd){
		Set(XStart, XEnd, YStart, YEnd, 0, 1);
	}

	GridRegion(int XStart, int XEnd, int YStart, int YEnd){
		Set(XStart, XEnd, YStart, YEnd, 0, 1);
	}

	GridRegion(long long XStart, long long XEnd){
		Set(XStart, XEnd, 0, 1, 0, 1);
	}

	GridRegion(int XStart, int XEnd){
		Set(XStart, XEnd, 0, 1, 0, 1);
	}

	GridRegion(){
		Set(0,1,0,1,0,1);
	}

	void SetRangeX(long long XStart, long long XEnd){
		x0 = XStart;
		x1 = XEnd;
	}

	void SetRangeY(long long YStart, long long YEnd){
		y0 = YStart;
		y1 = YEnd;
	}

	void SetRangeZ(long long ZStart, long long ZEnd){
		z0 = ZStart;
		z1 = ZEnd;
	}

	GridRegion Translate(long long dx, long long dy, long long dz){

		return GridRegion(x0+dx, x1+dx, y0+dy, y1+dy, z0+dz, z1+dz);
	
	}

	long long SizeX(){
		return x1 - x0;
	}

	long long SizeY(){
		return y1 - y0;
	}

	long long SizeZ(){
		return z1 - z0;
	}

	long long CountCells(){
		return (z1 - z0)*(y1 - y0)*(x1 - x0);
	}

	bool IsZero(){	//Test if all coordinates are 0

		return (x0 == 0 && x1 == 0 && y0 == 0 && y1 == 0 && z0 == 0 && z1 == 0);

	}

	bool Contains(long long x, long long y, long long z){		//Test if coordinate within region

		return (x >= x0 && x < x1 && y >= y0 && y < y1 && z >= z0 && z < z1);

	}

	bool TestOverlap(GridRegion& Region){

		return TestOverlap(&Region);
	}

	bool TestOverlap(GridRegion* Region){	//Test whether regions overlap

		long long xn_a = (x0 <= x1) ? x0 : x1;		//Ensure works even if x1 < x0 etc.
		long long xp_a = (x1 >= x0) ? x1 : x0;
		long long yn_a = (y0 <= y1) ? y0 : y1;
		long long yp_a = (y1 >= y0) ? y1 : y0;
		long long zn_a = (z0 <= z1) ? z0 : z1;
		long long zp_a = (z1 >= z0) ? z1 : z0;

		long long xn_b = (Region->x0 <= Region->x1) ? Region->x0 : Region->x1;
		long long xp_b = (Region->x1 >= Region->x0) ? Region->x1 : Region->x0;
		long long yn_b = (Region->y0 <= Region->y1) ? Region->y0 : Region->y1;
		long long yp_b = (Region->y1 >= Region->y0) ? Region->y1 : Region->y0;
		long long zn_b = (Region->z0 <= Region->z1) ? Region->z0 : Region->z1;
		long long zp_b = (Region->z1 >= Region->z0) ? Region->z1 : Region->z0;

		long long x0_ = max(xn_a, xn_b);
		long long x1_ = min(xp_a, xp_b);
		long long y0_ = max(yn_a, yn_b);
		long long y1_ = min(yp_a, yp_b);
		long long z0_ = max(zn_a, zn_b);
		long long z1_ = min(zp_a, zp_b);

		if(x1_ <= x0_ || y1_ <= y0_ || z1_ <= z0_)		//No overlap
			return false;

		return true;
	}

	GridRegion Overlap(GridRegion& Region){

		return Overlap(&Region);
	}

	GridRegion Overlap(GridRegion* Region){

		long long xn_a = (x0 <= x1) ? x0 : x1;		//Ensure works even if x1 < x0 etc.
		long long xp_a = (x1 >= x0) ? x1 : x0;
		long long yn_a = (y0 <= y1) ? y0 : y1;
		long long yp_a = (y1 >= y0) ? y1 : y0;
		long long zn_a = (z0 <= z1) ? z0 : z1;
		long long zp_a = (z1 >= z0) ? z1 : z0;

		long long xn_b = (Region->x0 <= Region->x1) ? Region->x0 : Region->x1;
		long long xp_b = (Region->x1 >= Region->x0) ? Region->x1 : Region->x0;
		long long yn_b = (Region->y0 <= Region->y1) ? Region->y0 : Region->y1;
		long long yp_b = (Region->y1 >= Region->y0) ? Region->y1 : Region->y0;
		long long zn_b = (Region->z0 <= Region->z1) ? Region->z0 : Region->z1;
		long long zp_b = (Region->z1 >= Region->z0) ? Region->z1 : Region->z0;

		long long x0_ = max(xn_a, xn_b);
		long long x1_ = min(xp_a, xp_b);
		long long y0_ = max(yn_a, yn_b);
		long long y1_ = min(yp_a, yp_b);
		long long z0_ = max(zn_a, zn_b);
		long long z1_ = min(zp_a, zp_b);

		if(x1_ <= x0_ || y1_ <= y0_ || z1_ <= z0_)		//No overlap, return 0
			return GridRegion(0,0,0,0,0,0);

		return GridRegion(x0_, x1_, y0_, y1_, z0_, z1_);
	}

	bool FitsWithin(GridRegion& Region){

		return FitsWithin(&Region);
	}

	bool FitsWithin(GridRegion* Region){		//Test to see if fits entirely inside other region

		bool CondX = (x0 >= Region->x0 && x0 <= Region->x1) && (x1 >= Region->x0 && x1 <= Region->x1);
		bool CondY = (y0 >= Region->y0 && y0 <= Region->y1) && (y1 >= Region->y0 && y1 <= Region->y1);
		bool CondZ = (z0 >= Region->z0 && z0 <= Region->z1) && (z1 >= Region->z0 && z1 <= Region->z1);

		return (CondX && CondY && CondZ);
	}

	char* ToStr(){

		sprintf(OutStr, "%lli -> %lli, %lli -> %lli, %lli -> %lli", x0, x1, y0, y1, z0, z1);

		return OutStr;
	}
};

char GridRegion::OutStr[256];	//Static char for writing region to string

class VectorDisp{

public:
	char Str[128];

	VectorDisp(GridPosition Pos){

		sprintf(Str, "%lli, %lli, %lli", Pos.x, Pos.y, Pos.z);

	}

	VectorDisp(GridRegion Region){

		sprintf(Str, "%lli -> %lli, %lli -> %lli, %lli -> %lli", Region.x0, Region.x1, Region.y0, Region.y1, Region.z0, Region.z1);

	}

};