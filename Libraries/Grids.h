#pragma once

#include <typeinfo>

#include "Base.h"
#include "DataFiles.h"		//Read and write large arrays in files



/*!
	Specifies new grid size and options for a coarsegraining operation within the Grid class
*/

struct Coarsegrain{

	long long SzX;		//! New grid size in X
	long long SzY;		//! New grid size in X
	long long SzZ;		//! New grid size in Y

	bool Discrete;		//! Discrete or Continuous mode. Discrete fills a coarsened voxel with the mode overlap value; Continuous fills a coarsened voxel with a value average

private:

	void Set(long long SizeX, long long SizeY, long long SizeZ, bool DiscreteMode){

		SzX = SizeX;
		SzY = SizeY;
		SzZ = SizeZ;

		Discrete = DiscreteMode;

	}

public:

	Coarsegrain(long long SizeX, long long SizeY, long long SizeZ, bool DiscreteMode){
		Set(SizeX, SizeY, SizeZ, DiscreteMode);
	}

	Coarsegrain(int SizeX, int SizeY, int SizeZ, bool DiscreteMode){
		Set(SizeX, SizeY, SizeZ, DiscreteMode);
	}

	Coarsegrain(long long SizeX, long long SizeY, long long SizeZ){
		Set(SizeX, SizeY, SizeZ, false);
	}

	Coarsegrain(int SizeX, int SizeY, int SizeZ){
		Set(SizeX, SizeY, SizeZ, false);
	}

};


//Grid Class#

/*!
	A class for manipulating 3D datasets.
	Includes operations for iterating values, extracting subregions, coarsegraining, reading from/writing to files etc.
*/

template <class T>
class Grid{

public:
	T* Data;			//! Grid data


private:

	long long Sz[4];	//Array dimensions

	bool FreeArr;		//Free array at end (can set false if using shallow copy)

	//

	void Init(long long SzX, long long SzY, long long SzZ, long long Count){

		Sz[0] = Count;
		Sz[1] = SzX;
		Sz[2] = SzY;
		Sz[3] = SzZ;
		
		Data = new T[Sz[0]*Sz[1]*Sz[2]*Sz[3]];

		FreeArr = true;		//Array freed at destructor

		memset(Data, 0, sizeof(T)*Sz[0]*Sz[1]*Sz[2]*Sz[3]);
	}

	long long ArrIndex(long long X, long long Y, long long Z, long long n){

		return Sz[0]*(Sz[1]*(Sz[2]*Z + Y) + X) + n;

	}

	T GetEntry(long long X, long long Y, long long Z, long long n){

		return Data[ArrIndex(X, Y, Z, n)];

	}

	T& GetEntryRef(long long X, long long Y, long long Z, long long n){

		return Data[ArrIndex(X, Y, Z, n)];

	}

	T* GetEntryPtr(long long X, long long Y, long long Z, long long n){

		return &(Data[ArrIndex(X, Y, Z, n)]);

	}

	void SetEntry(long long X, long long Y, long long Z, long long n, T Value){

		Data[ArrIndex(X, Y, Z, n)] = Value;

	}

	DataFileBase::DataType ObtainDataType(const type_info &Tp){

		if(Tp==typeid(int))					return DataFileBase::Type_int32;
		if(Tp==typeid(unsigned int))		return DataFileBase::Type_uint32;
		if(Tp==typeid(short))				return DataFileBase::Type_int16;
		if(Tp==typeid(unsigned short))		return DataFileBase::Type_uint16;
		if(Tp==typeid(long long))			return DataFileBase::Type_int64;
		if(Tp==typeid(unsigned long long))	return DataFileBase::Type_uint64;
		if(Tp==typeid(float))				return DataFileBase::Type_float;
		if(Tp==typeid(double))				return DataFileBase::Type_double;
		if(Tp==typeid(char))				return DataFileBase::Type_char;
		if(Tp==typeid(bool))				return DataFileBase::Type_bool;

		return DataFileBase::Type_char;
	}

	//Subregion

	int ValidateRegion(Grid* RefGrid, GridRegion& Region){

		if(Region.SizeX() <= 0 || Region.SizeY() <= 0 || Region.SizeZ() <= 0){

			return 2;	//Unusable region

		}

		//Check bounds

		int retCode = 0;

		if(Region.x0 < 0 || Region.y0 < 0 || Region.z0 < 0 || Region.x1 > RefGrid->SizeX() || Region.y1 > RefGrid->SizeY() || Region.z1 > RefGrid->SizeZ()){
			
			retCode = 1;	//Usable with adjustment of bounds

		}

		if(Region.x0 < 0) Region.x0 = 0;
		if(Region.y0 < 0) Region.y0 = 0;
		if(Region.z0 < 0) Region.z0 = 0;

		if(Region.x1 > RefGrid->SizeX()) Region.x1 = RefGrid->SizeX();
		if(Region.y1 > RefGrid->SizeY()) Region.y1 = RefGrid->SizeY();
		if(Region.z1 > RefGrid->SizeZ()) Region.z1 = RefGrid->SizeZ();

		return retCode;		// 0 = valid; 1 = adjusted bounds; 2 = unusable
	}

	void ReadSubRegion(Grid<T> &GridRef, GridRegion Region, GridPosition Offset){

		bool Flip[3] = {false, false, false};

		ReadSubRegion(GridRef, Region, Offset, Flip);

	}

	void ReadSubRegion(Grid<T> &GridRef, GridRegion Region, GridPosition Offset, bool Flip[3]){

		for(long long z=0; z<Region.SizeZ(); z++)
		for(long long y=0; y<Region.SizeY(); y++)
		for(long long x=0; x<Region.SizeX(); x++)
		for(long long n=0; n<Sz[0]; n++){

			long long _x = Flip[0] ? Region.x1-x-1 : x+Region.x0;
			long long _y = Flip[1] ? Region.y1-y-1 : y+Region.y0;
			long long _z = Flip[2] ? Region.z1-z-1 : z+Region.z0;

			SetEntry(x+Offset.x, y+Offset.y, z+Offset.z, n, GridRef.Get(_x, _y, _z, n));

		}

	}

	//Coarsegraining

	void ReadCoarsegrained(Grid<T> &GridRef, Coarsegrain CG_Op, GridPosition Offset){

		double CoarseFactorX = (double)GridRef.SizeX() / (double)CG_Op.SzX;
		double CoarseFactorY = (double)GridRef.SizeY() / (double)CG_Op.SzY;
		double CoarseFactorZ = (double)GridRef.SizeZ() / (double)CG_Op.SzZ;

		//Sample arrays

		long long nMaxX = 1 + (long long)ceil(CoarseFactorX);
		long long nMaxY = 1 + (long long)ceil(CoarseFactorY);
		long long nMaxZ = 1 + (long long)ceil(CoarseFactorZ);

		long long nArr = nMaxX*nMaxY*nMaxZ;

		T* OverlapGrid = new T[nArr*SizeVector()];		//Overlap values
		T* UniqueValues = new T[nArr];					//Unique values
		double* ValueCount = new double[nArr];			//Counting value instances
		double* Weights = new double[nArr];				//Overlap fractions

		//Coarsegrain

		for(long long zi=0; zi<CG_Op.SzZ; zi++)		
		for(long long yi=0; yi<CG_Op.SzY; yi++)
		for(long long xi=0; xi<CG_Op.SzX; xi++){	//Each cell in new grid

			//Region of fine grid to coarsegrain to cell

			double pX0 = CoarseFactorX * (double)xi;
			double pX1 = pX0 + CoarseFactorX;
			double pY0 = CoarseFactorY * (double)yi;
			double pY1 = pY0 + CoarseFactorY;
			double pZ0 = CoarseFactorZ * (double)zi;
			double pZ1 = pZ0 + CoarseFactorZ;

			long long x0 = (long long)floor(pX0);
			long long x1 = min( (long long)ceil(pX1) , GridRef.SizeX() );
			long long y0 = (long long)floor(pY0);
			long long y1 = min( (long long)ceil(pY1) , GridRef.SizeY() );
			long long z0 = (long long)floor(pZ0);
			long long z1 = min( (long long)ceil(pZ1) , GridRef.SizeZ() );

			long long dx = x1 - x0;
			long long dy = y1 - y0;
			long long dz = z1 - z0;

			for(long long z=z0; z<z1; z++)
			for(long long y=y0; y<y1; y++)
			for(long long x=x0; x<x1; x++){

				long long _x = x-x0;
				long long _y = y-y0;
				long long _z = z-z0;

				double xFill = min( min ( (double)x + 1.0 - pX0 , pX1 - (double)x ) , 1.0 );
				double yFill = min( min ( (double)y + 1.0 - pY0 , pY1 - (double)y ) , 1.0 );
				double zFill = min( min ( (double)z + 1.0 - pZ0 , pZ1 - (double)z ) , 1.0 );

				long long GridIndex = dx*(dy*_z + _y) + _x;

				for(long long n=0; n<SizeVector(); n++){

					OverlapGrid[GridIndex*SizeVector() + n] = GridRef.Get(x, y, z, n);

				}

				Weights[GridIndex] = xFill*yFill*zFill;

			}

			CoarsegrainFillValue(OverlapGrid, UniqueValues, Weights, ValueCount, dx*dy*dz, CG_Op.Discrete, xi+Offset.x, yi+Offset.y, zi+Offset.z);

		}

		delete[] OverlapGrid;
		delete[] UniqueValues;
		delete[] ValueCount;
		delete[] Weights;

	}

	void CoarsegrainFillValue(T* OverlapGrid, T* UniqueValues, double* Weights, double* ValueInstances, long long ArrLength, bool Discrete, long long GridX, long long GridY, long long GridZ){

		for(long long n=0; n<SizeVector(); n++){

			if(Discrete){

				//Fill with mode value

				long long nUniqueValues = 0;

				for(long long i=0; i<ArrLength; i++){

					T Value = OverlapGrid[i*SizeVector() + n];

					//Search values list

					bool Found = false;

					for(long long i2=0; i2<nUniqueValues; i2++){

						if(Value == UniqueValues[i2]){
							Found = true;
							ValueInstances[i2] += Weights[i];
						}
					}

					if(!Found){
						UniqueValues[nUniqueValues] = Value;
						ValueInstances[nUniqueValues] = Weights[i];
						nUniqueValues++;
					}

				}

				//Find most common value

				T ModeValue = UniqueValues[0];
				double ModeWeight = ValueInstances[0];

				for(long long c=1; c<nUniqueValues; c++){
					if(ValueInstances[c] > ModeWeight){
						ModeValue = UniqueValues[c];
						ModeWeight = ValueInstances[c];
					}
				}

				//Fill

				Set(GridX, GridY, GridZ, n, ModeValue);

			}else{

				//Average values

				bool IntAverage = true;			//Use int or float arithmetic for averaging

				if(typeid(T)==typeid(float) || typeid(T)==typeid(double))
					IntAverage = false;

				double AvSum = 0;
				double WeightSum = 0;

				//Averaging

				for(long long i=0; i<ArrLength; i++){

					T Value = OverlapGrid[i*SizeVector() + n];

					if(IntAverage){
						AvSum += (((long long)Value) * Weights[i]);
					}else{
						AvSum += (((double)Value) * Weights[i]);
					}

					WeightSum += Weights[i];

				}

				//Fill

				double Av = AvSum/WeightSum;

				if(IntAverage){
					
					long long IntAv = (long long)round(Av);

#pragma warning(suppress: 4800)			//Disable warning about casting from int to bool for next line

					Set(GridX, GridY, GridZ, n, (T)IntAv);

				}else{

					Set(GridX, GridY, GridZ, n, (T)Av);

				}

			}

		}

	}

public:

	//Constructors

	Grid(int SizeX,  int SizeY, int SizeZ, int VectorSize){								//! Initialise a 3D grid of size SizeX x SizeY x SizeZ with vector entries of dimension VectorSize
		Init(SizeX, SizeY, SizeZ, VectorSize);
	}

	Grid(long long SizeX,  long long SizeY, long long SizeZ, long long VectorSize){		//! Initialise a 3D grid of size SizeX x SizeY x SizeZ with vector entries of dimension VectorSize
		Init(SizeX, SizeY, SizeZ, VectorSize);
	}

	Grid(int SizeX,  int SizeY, int SizeZ){							//! Initialise a 3D grid of size SizeX x SizeY x SizeZ with scalar entries
		Init(SizeX, SizeY, SizeZ, 1);
	}

	Grid(long long SizeX,  long long SizeY, long long SizeZ){		//! Initialise a 3D grid of size SizeX x SizeY x SizeZ with scalar entries
		Init(SizeX, SizeY, SizeZ, 1);
	}

	Grid(int SizeX,  int SizeY){					//! Initialise a 2D grid of size SizeX x SizeY with scalar entries
		Init(SizeX, SizeY, 1, 1);
	}

	Grid(long long SizeX,  long long SizeY){		//! Initialise a 2D grid of size SizeX x SizeY with scalar entries
		Init(SizeX, SizeY, 1, 1);
	}

	Grid(int SizeX){								//! Initialise an array of size SizeX with scalar entries
		Init(SizeX, 1, 1, 1);
	}

	Grid(long long SizeX){							//! Initialise an array of size SizeX with scalar entries
		Init(SizeX, 1, 1, 1);
	}

	Grid(Grid<T> &GridRef){		//! Create a duplicate of GridRef
		Init(GridRef.SizeX(), GridRef.SizeY(), GridRef.SizeZ(), GridRef.SizeVector());
		memcpy(Data, GridRef.Data, sizeof(T)*Sz[0]*Sz[1]*Sz[2]*Sz[3]);
	}

	Grid(Grid<T>* GridRef){		//! Create a duplicate of GridRef
		Init(GridRef->SizeX(), GridRef->SizeY(), GridRef->SizeZ(), GridRef->SizeVector());
		memcpy(Data, GridRef->Data, sizeof(T)*Sz[0]*Sz[1]*Sz[2]*Sz[3]);
	}

	Grid(Grid<T> &GridRef, GridRegion Region){		//! Create a new grid containing a subregion of GridRef
		Init(Region.SizeX(), Region.SizeY(), Region.SizeZ(), GridRef.SizeVector());
		ReadSubRegion(GridRef, Region, GridPosition(0,0,0));
	}

	Grid(Grid<T> &GridRef, Coarsegrain CG_Op){		//! Create a new grid containing a coarsegained version of GridRef
		Init(CG_Op.SzX, CG_Op.SzY, CG_Op.SzZ, GridRef.SizeVector());
		ReadCoarsegrained(GridRef, CG_Op, GridPosition(0,0,0));
	}

	Grid(T* DataArr, bool ShallowCopy, bool FreeShallowCopy, long long SizeX,  long long SizeY, long long SizeZ, long long VectorSize){
	
		Sz[0] = VectorSize;
		Sz[1] = SizeX;
		Sz[2] = SizeY;
		Sz[3] = SizeZ;

		if(ShallowCopy){
		
			Data = DataArr;				//Copy ptr

			FreeArr = FreeShallowCopy;	//Whether to free at end
		
		}else{
		
			Data = new T[Sz[0]*Sz[1]*Sz[2]*Sz[3]];		//New array

			FreeArr = true;								//Array freed at destructor

			memcpy(Data, DataArr, sizeof(T)*Sz[0]*Sz[1]*Sz[2]*Sz[3]);	//Copy data
		
		}

	}


	//Destructor

	~Grid(){
		if(FreeArr)
			delete[] Data;
	}


	//Metrics

	long long SizeX() const {				//! Returns grid size in X
		return Sz[1];
	}

	long long SizeY() const {				//! Returns grid size in Y
		return Sz[2];
	}

	long long SizeZ() const {				//! Returns grid size in Z
		return Sz[3];
	}

	long long SizeVector() const {		//! Returns grid vector size
		return Sz[0];
	}

	long long CountOccurences(T Value[], GridRegion Region){		//! Counts instances of vectors within a given region

		long long c = 0;

		for(long long z=Region.z0; z<Region.z1; z++)
		for(long long y=Region.y0; y<Region.y1; y++)
		for(long long x=Region.x0; x<Region.x1; x++){

			bool Match = true;

			for(long long n=0; n<SizeVector(); n++){
				if(Get(x, y, z, n) != Value[n])
					Match = false;
			}

			if(Match)
				c++;
		}
			
		return c;
	}

	long long CountOccurences(T Value[]){		//! Counts instances of vectors

		return CountOccurences(Value, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	long long CountOccurences(T* ValueList, int nValues, long long* CountOut, GridRegion Region){

		long long c = 0;

		memset(CountOut, 0, sizeof(long long)*nValues);

		for(long long z=Region.z0; z<Region.z1; z++)
		for(long long y=Region.y0; y<Region.y1; y++)
		for(long long x=Region.x0; x<Region.x1; x++)
		for(long long n=0; n<SizeVector(); n++){

			T Val = Get(x, y, z, n);

			for(int i=0; i<nValues; i++)
				if(Val == ValueList[i]){
					CountOut[i]++;
					c++;
				}
		}
			
		return c;
	}

	long long CountOccurences(T* ValueList, int nValues, long long* CountOut){

		return CountOccurences(ValueList, nValues, CountOut, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	long long CountOccurences(T Value, GridRegion Region){		//! Count instances of a value within a given region

		long long c = 0;

		for(long long z=Region.z0; z<Region.z1; z++)
		for(long long y=Region.y0; y<Region.y1; y++)
		for(long long x=Region.x0; x<Region.x1; x++)
		for(long long n=0; n<SizeVector(); n++){

			if(Get(x, y, z, n) == Value)
				c++;
		}
			
		return c;
	}

	long long CountOccurences(T Value){		//! Count instances of a value

		return CountOccurences(Value, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}


	//Get

	T& operator()(long long X, long long Y, long long Z, long long n){		//! Return reference to value at X, Y, Z; with vector component n
		return GetEntryRef(X, Y, Z, n);
	}

	T& operator()(int X, int Y, int Z, int n){								//! Return reference to value at X, Y, Z; with vector component n
		return GetEntryRef(X, Y, Z, n);
	}

	T& operator()(long long X, long long Y, long long Z){		//! Return reference to value at X, Y, Z
		return GetEntryRef(X, Y, Z, 0);
	}

	T& operator()(int X, int Y, int Z){							//! Return reference to value at X, Y, Z
		return GetEntryRef(X, Y, Z, 0);
	}

	T& operator()(long long X, long long Y){		//! Return reference to value at X, Y
		return GetEntryRef(X, Y, 0, 0);
	}

	T& operator()(int X, int Y){					//! Return reference to value at X, Y
		return GetEntryRef(X, Y, 0, 0);
	}

	T& operator()(long long X){			//! Return reference to value at X
		return GetEntryRef(X, 0, 0, 0);
	}

	T& operator()(int X){				//! Return reference to value at X
		return GetEntryRef(X, 0, 0, 0);
	}

	T Get(long long X, long long Y, long long Z, long long n){		//! Return value at X, Y, Z; with vector component n
		return GetEntry(X, Y, Z, n);
	}

	T Get(int X, int Y, int Z, int n){								//! Return value at X, Y, Z; with vector component n
		return GetEntry(X, Y, Z, n);
	}

	T Get(long long X, long long Y, long long Z){		//! Return value at X, Y, Z
		return GetEntry(X, Y, Z, 0);
	}

	T Get(int X, int Y, int Z){							//! Return value at X, Y, Z
		return GetEntry(X, Y, Z, 0);
	}

	T Get(long long X, long long Y){		//! Return value at X, Y
		return GetEntry(X, Y, 0, 0);
	}

	T Get(int X, int Y){					//! Return value at X, Y
		return GetEntry(X, Y, 0, 0);
	}

	T Get(long long X){			//! Return value at X
		return GetEntry(X, 0, 0, 0);
	}

	T Get(int X){				//! Return value at X
		return GetEntry(X, 0, 0, 0);
	}

	T* GetPtr(long long X, long long Y, long long Z, long long n){		//! Return value at X, Y, Z; with vector component n
		return GetEntryPtr(X, Y, Z, n);
	}

	T* GetPtr(int X, int Y, int Z, int n){								//! Return value at X, Y, Z; with vector component n
		return GetEntryPtr(X, Y, Z, n);
	}

	T* GetPtr(long long X, long long Y, long long Z){		//! Return value at X, Y, Z
		return GetEntryPtr(X, Y, Z, 0);
	}

	T* GetPtr(int X, int Y, int Z){							//! Return value at X, Y, Z
		return GetEntryPtr(X, Y, Z, 0);
	}

	T* GetPtr(long long X, long long Y){		//! Return value at X, Y
		return GetEntryPtr(X, Y, 0, 0);
	}

	T* GetPtr(int X, int Y){					//! Return value at X, Y
		return GetEntryPtr(X, Y, 0, 0);
	}

	T* GetPtr(long long X){			//! Return value at X
		return GetEntryPtr(X, 0, 0, 0);
	}

	T* GetPtr(int X){				//! Return value at X
		return GetEntryPtr(X, 0, 0, 0);
	}


	//Set

	void Set(long long X, long long Y, long long Z, long long n, T Value){		//! Set value at X, Y, Z; vector component n
		SetEntry(X, Y, Z, n, Value);
	}

	void Set(int X, int Y, int Z, int n, T Value){								//! Set value at X, Y, Z; vector component n
		SetEntry(X, Y, Z, n, Value);
	}

	void Set(long long X, long long Y, long long Z, T Value){		//! Set value at X, Y, Z
		SetEntry(X, Y, Z, 0, Value);
	}

	void Set(int X, int Y, int Z, T Value){							//! Set value at X, Y, Z
		SetEntry(X, Y, Z, 0, Value);
	}

	void Set(long long X, long long Y, T Value){		//! Set value at X, Y
		SetEntry(X, Y, 0, 0, Value);
	}

	void Set(int X, int Y, T Value){					//! Set value at X, Y
		SetEntry(X, Y, 0, 0, Value);
	}

	void Set(long long X, T Value){		//! Set value at X
		SetEntry(X, 0, 0, 0, Value);
	}

	void Set(int X, T Value){			//! Set value at X
		SetEntry(X, 0, 0, 0, Value);
	}

	void SetAll(T Value, GridRegion Region){				//! Fill region with given value

		for(long long z=Region.z0; z<Region.z1; z++)
		for(long long y=Region.y0; y<Region.y1; y++)
		for(long long x=Region.x0; x<Region.x1; x++)
		for(long long n=0; n<SizeVector(); n++)
			SetEntry(x, y, z, n, Value);

	}

	void SetAll(T Value[], GridRegion Region){				//! Fill region with given value vector

		for(long long z=Region.z0; z<Region.z1; z++)
		for(long long y=Region.y0; y<Region.y1; y++)
		for(long long x=Region.x0; x<Region.x1; x++)
		for(long long n=0; n<SizeVector(); n++)
			SetEntry(x, y, z, n, Value[n]);

	}

	void SetAll(T Value){									//! Fill grid with given value

		SetAll(Value, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	void SetAll(T Value[]){									//! Fill grid with given value

		SetAll(Value, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	void ReplaceValue(T ValueFind, T ValueReplace, GridRegion Region){		//! Replace value within given region

		for(long long z=Region.z0; z<Region.z1; z++)
		for(long long y=Region.y0; y<Region.y1; y++)
		for(long long x=Region.x0; x<Region.x1; x++)
		for(long long n=0; n<SizeVector(); n++){

			if(Get(x, y, z, n) == ValueFind)
				SetEntry(x, y, z, n, ValueReplace);
		}
	}

	void ReplaceValue(T ValueFind[], T ValueReplace[], GridRegion Region){		//! Replace value vector within given region

		for(long long z=Region.z0; z<Region.z1; z++)
		for(long long y=Region.y0; y<Region.y1; y++)
		for(long long x=Region.x0; x<Region.x1; x++){

			bool Match = true;

			for(long long n=0; n<SizeVector(); n++){
				if(Get(x, y, z, n) != ValueFind[n])
					Match = false;
			}

			if(Match){
				for(long long n=0; n<SizeVector(); n++)
					SetEntry(x, y, z, n, ValueReplace[n]);
			}
		}
	}

	void ReplaceValue(T ValueFind, T ValueReplace){				//! Replace ValueFind with ValueReplace

		ReplaceValue(ValueFind, ValueReplace, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	void ReplaceValue(T ValueFind[], T ValueReplace[]){			//! Replace ValueFind vector with ValueReplace vector

		ReplaceValue(ValueFind, ValueReplace, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	void ScaleValues(T Scale, GridRegion Region){	//! Multiply values by Scale within a given region

		for(long long z=Region.z0; z<Region.z1; z++)
		for(long long y=Region.y0; y<Region.y1; y++)
		for(long long x=Region.x0; x<Region.x1; x++)
		for(long long n=0; n<SizeVector(); n++){
			GetEntryRef(x,y,z,n) *= Scale;
		}
	}

	void ScaleValues(T Scale){						//! Multiply values by Scale
		ScaleValues(Scale, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));
	}


	//Grid operations

	void ReadGridCoarsegrained(Grid& GridRef, Coarsegrain CG_Op, GridPosition Offset){	//! Read GridRef coarsegrained and place at Offset in grid

		ReadCoarsegrained(GridRef, CG_Op, Offset);

	}

	void ReadGridCoarsegrained(Grid& GridRef, Coarsegrain CG_Op){	//! Read GridRef coarsegrained and place at (0,0,0) in grid

		ReadCoarsegrained(GridRef, CG_Op, GridPosition(0,0,0));

	}

	void ReadSubGrid(Grid& GridRef, GridRegion Region, GridPosition Offset, bool Flip[3]){

		ReadSubRegion(GridRef, Region, Offset, Flip);
	}

	void ReadSubGrid(Grid& GridRef, GridRegion Region, GridPosition Offset){	//! Read a given Region of GridRef and place at Offset in grid

		ReadSubRegion(GridRef, Region, Offset);

	}

	void ReadSubGrid(Grid& GridRef, GridRegion Region){		//! Read a given Region of GridRef

		ReadSubRegion(GridRef, Region, GridPosition(0LL,0LL,0LL));

	}

	void ReadGrid(Grid& GridRef, GridPosition Offset, bool Flip[3]){

		ReadSubRegion(GridRef, GridRegion(0LL,GridRef.SizeX(),0LL,GridRef.SizeY(),0LL,GridRef.SizeZ()), Offset, Flip);
	}

	void ReadGrid(Grid& GridRef, GridPosition Offset){		//! Read GridRef into grid at a given offset

		ReadSubRegion(GridRef, GridRegion(0LL,GridRef.SizeX(),0LL,GridRef.SizeY(),0LL,GridRef.SizeZ()), Offset);

	}

	void ReadGrid(Grid& GridRef){		//! Read GridRef into grid

		ReadSubRegion(GridRef, GridRegion(0LL,GridRef.SizeX(),0LL,GridRef.SizeY(),0LL,GridRef.SizeZ()), GridPosition(0LL,0LL,0LL));

	}


	//Outputs

	void PrintValues(const char* SeparatorX, const char* SeparatorY, const char* SeparatorZ, const char* SeparatorN){	//! Print values to console with specified separator strings

		const char* PrtStr;

		if(typeid(T) == typeid(int))					 PrtStr = "%i";
		else if(typeid(T) == typeid(long long))			 PrtStr = "%lli";
		else if(typeid(T) == typeid(unsigned long long)) PrtStr = "%llu";
		else if(typeid(T) == typeid(float))				 PrtStr = "%f";
		else if(typeid(T) == typeid(double))			 PrtStr = "%f";
		else if(typeid(T) == typeid(char))				 PrtStr = "%c";
		else if(typeid(T) == typeid(bool))				 PrtStr = "%i";
		else											 return;

		for(int z=0; z<Sz[3]; z++){
		for(int y=0; y<Sz[2]; y++){
		for(int x=0; x<Sz[1]; x++){
		for(int n=0; n<Sz[0]; n++){

			printf(PrtStr, Get(x, y, z, n));

			if(n<Sz[0]-1)
				printf("%s", SeparatorN);
		}
		if(x<Sz[1]-1)
			printf("%s", SeparatorX);
		}
		if(y<Sz[2]-1)
			printf("%s", SeparatorY);
		}
		if(z<Sz[3]-1)
			printf("%s", SeparatorZ);
		}

	}

	void PrintValues(const char* SeparatorX, const char* SeparatorY, const char* SeparatorZ){		//! Print values to console with specified separator strings
		PrintValues(SeparatorX, SeparatorY, SeparatorZ, "");
	}

	void PrintValues(const char* SeparatorX, const char* SeparatorY){	//! Print values to console with specified separator strings
		PrintValues(SeparatorX, SeparatorY, "", "");
	}

	void PrintValues(const char* Separator){		//! Print values to console with specified separator string
		PrintValues(Separator, Separator, Separator, Separator);
	}

	void PrintValues(){		//! Print values to console with default separator strings
		PrintValues(" ","\n","\n\n",",");
	}

	
	//File data

	int ReadFromFile(char* FileName, bool Binary){

		return ReadFromFile(FileName, Binary, ObtainDataType(typeid(T)));
	}

	int ReadFromFile(char* FileName, bool Binary, DataType FileType){		//! Read dataset from file

		DataFileReader ReadData;

		//Open file

		bool ret = ReadData.Open(FileName);

		if(!ret)
			return 1;	//File could not be opened

		//Read operation

		long long GridSize[3];
		
		GridSize[0] = Sz[1];	//x
		GridSize[1] = Sz[2];	//y
		GridSize[2] = Sz[3];	//z

		ArrayReadOp ReadOp(Binary, GridSize, (int)Sz[0], FileType, Data, ObtainDataType(typeid(T)));

		//Read

		ReadData.ReadArrayFromEnd(&ReadOp);

		if(ReadData.LastError() == 0){

			return 0;

		}

		if(ReadData.LastError() == 2) return 2;		//Insufficient data
		if(ReadData.LastError() == 3) return 3;		//Invalid data
		if(ReadData.LastError() == 4) return 4;		//Memory allocation error
		
		return 5;
	}

	int WriteToFile(const char* FileName,		/** File path */
					bool Binary,				/** File in binary format, else ASCII text */
					bool VTKHeader				/** Has VTK header */ ){									//! Write dataset to file

		DataFileWriter WriteData;

		bool ret = WriteData.Open(FileName);

		if(!ret)
			return 1;	//File could not be opened

		//Write operation

		long long GridSize[3];
		
		GridSize[0] = Sz[1];	//x
		GridSize[1] = Sz[2];	//y
		GridSize[2] = Sz[3];	//z

		ArrayWriteOp WriteOp(Binary, GridSize, (int)Sz[0], ObtainDataType(typeid(T)), Data, ObtainDataType(typeid(T)));

		if(VTKHeader){
			WriteOp.AddVTKHeader("Geometry Output", "Geometry");
		}

		//Write
		
		bool s = WriteData.WriteArray(&WriteOp);						//Output

		if(s){
			return 0;
		}

		return 5;
	}

	//Error codes

	static const char* ErrorCodeText(int ECode){		//! Returns an error message based on an error code

		const char* EStrings[] = {
			"The operation completed succesfully",				//0
			"The file could not be opened",						//1
			"The file contained insufficient data",				//2
			"The file contained invalid data",					//3
			"Memory allocation error",							//4
			"Unknown error",									//5
		};

		return EStrings[ECode];
	}

};