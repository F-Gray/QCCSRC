#pragma once

#include <typeinfo>

#include "Base.h"
#include "DataFiles.h"		//Read and write large arrays in files
#include "Grids.h"

#include "mpi.h"


//GridMPI class

namespace GridMPIDefinitions{

	const int Q26Basis[26][3] = {
		{ 1, 0, 0},		//+X			[0]
		{-1, 0, 0},		//-X			[1]
		{ 0, 1, 0},		//   +Y			[2]
		{ 0,-1, 0},		//   -Y			[3]
		{ 0, 0, 1},		//      +Z		[4]
		{ 0, 0,-1},		//      -Z		[5]
		{ 1, 1, 0},		//+X +Y			[6]
		{-1, 1, 0},		//-X +Y			[7]
		{ 1,-1, 0},		//+X -Y			[8]
		{-1,-1, 0},		//-X -Y			[9]
		{ 0, 1, 1},		//   +Y +Z		[10]
		{ 0,-1, 1},		//   -Y +Z		[11]
		{ 0, 1,-1},		//   +Y -Z		[12]
		{ 0,-1,-1},		//   -Y -Z		[13]
		{ 1, 0, 1},		//+X    +Z		[14]
		{-1, 0, 1},		//-X    +Z		[15]
		{ 1, 0,-1},		//+X    -Z		[16]
		{-1, 0,-1},		//-X    -Z		[17]
		{ 1, 1, 1},		//+X +Y +Z		[18]
		{-1, 1, 1},		//-X +Y +Z		[19]
		{ 1,-1, 1},		//+X -Y +Z		[20]
		{-1,-1, 1},		//-X -Y +Z		[21]
		{ 1, 1,-1},		//+X +Y -Z		[22]
		{-1, 1,-1},		//-X +Y -Z		[23]
		{ 1,-1,-1},		//+X -Y -Z		[24]
		{-1,-1,-1}		//-X -Y -Z		[25]
	};

};

template <class T>
class GridMPI{

public:
	T* Data;			//Grid data
	Grid<T>* GridLocal;	//Grid container for local Data

	MPI_Comm Comm;
	int Tag;

	GridRegion* Regions;

private:

	long long BufSizeMax;

	long long Sz[4];		//Array dimensions
	long long SzLocal[4];

	bool FreeArr;		//Free array at end (can set false if using shallow copy)

	int nThreads;		//Number of threads
	int tId;			//Local thread ID

	struct TransferGroup{
		vector<int> MainThreads;			//List of main threads
		vector<vector<int> > ThreadList;		//Lists of sub threads for each main thread
	};

	//

	void ComputeRegions(){

		//3D decomp

		double RatioXY = (double)Sz[1] / (double)Sz[2];
		double RatioXZ = (double)Sz[1] / (double)Sz[3];

		double N = max(1.0, pow( (double)nThreads*RatioXY*RatioXZ , 1.0/3.0));

		int nDivX = min(nThreads, (int)round(N));								//Number of partitions in X

		int nSubDivX = nThreads / nDivX;
		int rmSubDivX = nThreads % nDivX;

		double DivSzX = (double)Sz[1] / (double)nThreads;

		int nXCount = 0;

		for(int i=0; i<nDivX; i++){

			int tCountDivX = nSubDivX + (i<rmSubDivX ? 1 : 0);
			int SubCount = (int)round(tCountDivX*DivSzX);

			for(int c=0; c<tCountDivX; c++){

				Regions[nXCount+c].x0 = i==0 ?			 0 : Regions[nXCount-1].x1;
				Regions[nXCount+c].x1 = i==nDivX-1 ? Sz[1] : Regions[nXCount+c].x0 + SubCount;

			}

			//Divide in Y

			double RatioYZ = (double)Sz[2] / (double)Sz[3];

			double NY = max(1.0, sqrt( (double)tCountDivX*RatioYZ ));

			int nDivY = min(tCountDivX, (int)round(NY));

			int* tCountDivY = new int[nDivY];

			int nSubDivY = tCountDivX / nDivY;
			int rmSubDivY = tCountDivX % nDivY;

			double DivSzY = (double)Sz[2] / (double)tCountDivX;

			int nYCount = 0;

			for(int i2=0; i2<nDivY; i2++){

				int tCountDivY = nSubDivY + (i2<rmSubDivY ? 1 : 0);
				int SubCountY = (int)round(tCountDivY*DivSzY);

				for(int c=0; c<tCountDivY; c++){

					Regions[nXCount+nYCount+c].y0 = i2==0 ?			  0 : Regions[nXCount+nYCount-1].y1;
					Regions[nXCount+nYCount+c].y1 = i2==nDivY-1 ? Sz[2] : Regions[nXCount+nYCount+c].y0 + SubCountY;

				}

				//Divide in Z

				int nZ = (int)Sz[3] / tCountDivY;
				int rmZ = (int)Sz[3] % tCountDivY;

				for(int i3=0; i3<tCountDivY; i3++){

					Regions[nXCount+nYCount+i3].z0 = nZ*i3 + min(i3, rmZ);
					Regions[nXCount+nYCount+i3].z1 = Regions[nXCount+nYCount+i3].z0 + nZ + (i3<rmZ ? 1 : 0);

				}

				nYCount += tCountDivY;
			}
			
			nXCount += tCountDivX;
		}

	/*	if(tId == 0)
			for(int i=0; i<nThreads; i++)
		printf("Thread %i region (%i -> %i, %i -> %i, %i -> %i)\n", i, Regions[i].x0, Regions[i].x1, Regions[i].y0, Regions[i].y1, Regions[i].z0, Regions[i].z1);	//*/

	}

	void Init(MPI_Comm MPIComm, int MPITag, long long SzX, long long SzY, long long SzZ, long long Count){

		const long long BufSizeMaxDefault = 1024*1024*32;		//Default per thread max buffer size

		memset(this, 0, sizeof(GridMPI));

		Comm = MPIComm;
		Tag = MPITag;

		MPI_Comm_size(Comm, &nThreads);	//Get number of threads running
		MPI_Comm_rank(Comm, &tId);		//Get thread ID (index from 0 -> nThreads-1)

		Regions = new GridRegion[nThreads];

		Sz[0] = Count;
		Sz[1] = SzX;
		Sz[2] = SzY;
		Sz[3] = SzZ;

		ComputeRegions();

		GridRegion& Rgn = Regions[tId];

		SzLocal[0] = Sz[0];
		SzLocal[1] = Rgn.SizeX();
		SzLocal[2] = Rgn.SizeY();
		SzLocal[3] = Rgn.SizeZ();

		long long NodeCount = Sz[0]*Rgn.SizeX()*Rgn.SizeY()*Rgn.SizeZ();

		int Success = 0;
		
		try{
		
			Data = new T[NodeCount];		//Local data grid

		}catch(exception){

			Success = 1;
		}

		int SuccessGlobal;
		MPI_Allreduce(&Success, &SuccessGlobal, 1, MPI_INT, MPI_MAX, Comm);

		if(Success!=0){

			this->~GridMPI();

			throw bad_alloc();
		}

		GridLocal = new Grid<T>(Data, true, false, Rgn.SizeX(), Rgn.SizeY(), Rgn.SizeZ(), Sz[0]);	//Grid container

		FreeArr = true;				//Array freed at destructor

		memset(Data, 0, sizeof(T)*NodeCount);

		BufSizeMax = BufSizeMaxDefault;
	}

	long long ArrIndex(long long X, long long Y, long long Z, long long n){

		return SzLocal[0]*(SzLocal[1]*(SzLocal[2]*Z + Y) + X) + n;

	}

	void SetEntry(long long X, long long Y, long long Z, long long n, T Value){

		if(Regions[tId].Contains(X, Y, Z)){

			int _x = X;
			int _y = Y;
			int _z = Z;

			CoordGlobalToLocal(&_x, &_y, &_z);

			Data[ArrIndex(_x, _y, _z, n)] = Value;

		}

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
		if(Tp==typeid(unsigned char))		return DataFileBase::Type_uchar;
		if(Tp==typeid(bool))				return DataFileBase::Type_bool;

		return DataFileBase::Type_char;
	}

	void GroupThreadOps(vector<TransferGroup>& GroupList, GridRegion* SubRegions, int* BoundaryList){

		using namespace GridMPIDefinitions;
		
	/*	for(int i=0; i<nThreads; i++){

			for(int i2=0; i2<nThreads; i2++){

				if(SubRegions[i].TestOverlap(Regions[i2])){	//Thread i2 required by thread i

					if(tId==0)
					printf(" Thread %i uses thread %i\n", i, i2);

				}
			}
		}

		//*/

		vector<bool> ThreadHasRun(nThreads, false);

		int nRun = nThreads;

		for(int i=0; i<nThreads; i++){

			if(SubRegions[i].CountCells()==0){
				ThreadHasRun[i] = true;
				nRun--;
			}

		}

		for(int Count=0; Count<nRun; ){

			//New group

			GroupList.push_back(TransferGroup());

			TransferGroup& Group = GroupList[GroupList.size() - 1];

			vector<int> ThreadInUse(nThreads, -1);		//Thread used by another
			vector<bool> ThreadRun(nThreads, false);	//Run in this set

			vector<int> ThreadList;
		
			for(int i=0; i<nThreads; i++){

				if(ThreadHasRun[i])
					continue;
			
				if(ThreadInUse[i] != -1)		//Thread in use by another
					continue;

				bool Run = true;
				ThreadList.clear();

				for(int iv=0; iv<27; iv++){		//Vector components for boundary layers

					if(iv > 0 && BoundaryList == 0)
						break;

					GridRegion Rgn = SubRegions[i];

					if(BoundaryList != 0){

						//Incorporate boundary layers

						int* BSz = &BoundaryList[i*6];	//Boundary sizes {+X, -X, +Y, -Y, +Z, -Z}

						bool BSzZero = (BSz[0] == 0 && BSz[1] == 0 && BSz[2] == 0 && BSz[3] == 0 && BSz[4] == 0 && BSz[5] == 0);

						if(iv > 0 && BSzZero)
							break;

						if(iv == 0 && !BSzZero){		//Main grid
			
							if(Rgn.x1 < Sz[1]) Rgn.x1 += BSz[0];			//Add boundary layers if no loop required
							if(Rgn.x0 > 0    ) Rgn.x0 -= BSz[1];
							if(Rgn.y1 < Sz[2]) Rgn.y1 += BSz[2];
							if(Rgn.y0 > 0    ) Rgn.y0 -= BSz[3];
							if(Rgn.z1 < Sz[3]) Rgn.z1 += BSz[4];
							if(Rgn.z0 > 0    ) Rgn.z0 -= BSz[5];

						}else if (!BSzZero){

							int v[3];
							memcpy(v, Q26Basis[iv-1], sizeof(int)*3);

							bool RqXp = (v[0] > 0 && Rgn.x1 == Sz[1] && BSz[0] > 0);
							bool RqXn = (v[0] < 0 && Rgn.x0 == 0     && BSz[1] > 0);
							bool RqYp = (v[1] > 0 && Rgn.y1 == Sz[2] && BSz[2] > 0);
							bool RqYn = (v[1] < 0 && Rgn.y0 == 0     && BSz[3] > 0);
							bool RqZp = (v[2] > 0 && Rgn.z1 == Sz[3] && BSz[4] > 0);
							bool RqZn = (v[2] < 0 && Rgn.z0 == 0     && BSz[5] > 0);

							if(!RqXp && !RqXn && !RqYp && !RqYn && !RqZp && !RqZn)	//Layer not over boundary
								continue;
			
							//Boundary region

							if(v[0] > 0) Rgn.SetRangeX(Rgn.x1, Rgn.x1 + BSz[0]*v[0]);
							if(v[1] > 0) Rgn.SetRangeY(Rgn.y1, Rgn.y1 + BSz[2]*v[1]);
							if(v[2] > 0) Rgn.SetRangeZ(Rgn.z1, Rgn.z1 + BSz[4]*v[2]);
						
							if(v[0] < 0) Rgn.SetRangeX(Rgn.x0 + BSz[1]*v[0], Rgn.x0);
							if(v[1] < 0) Rgn.SetRangeY(Rgn.y0 + BSz[3]*v[1], Rgn.y0);
							if(v[2] < 0) Rgn.SetRangeZ(Rgn.z0 + BSz[5]*v[2], Rgn.z0);
						
							//Loop boundaries

							long long Translation[3] = {0, 0, 0};

							if(Rgn.x0 < 0) Translation[0] += Sz[1];		//Loop any coordinates over boundaries
							if(Rgn.y0 < 0) Translation[1] += Sz[2];
							if(Rgn.z0 < 0) Translation[2] += Sz[3];

							if(Rgn.x1 > Sz[1]) Translation[0] -= Sz[1];
							if(Rgn.y1 > Sz[2]) Translation[1] -= Sz[2];
							if(Rgn.z1 > Sz[3]) Translation[2] -= Sz[3];

							Rgn = Rgn.Translate(Translation[0], Translation[1], Translation[2]);

						}

					}

					for(int i2=0; i2<nThreads; i2++){

						if(i==i2)
							continue;

						if(Rgn.TestOverlap(Regions[i2])){	//Thread i2 required by thread i

							if(ThreadRun[i2] || ThreadInUse[i2] != -1){		//Thread i2 runs, or already in use

								Run = false;

							}else{

								ThreadList.push_back(i2);
							}
						
						}

					}

					if(!Run)
						break;
				}

				if(Run){

					ThreadRun[i] = true;
					ThreadHasRun[i] = true;

					for(int c=0; c<ThreadList.size(); c++)
						ThreadInUse[ThreadList[c]] = i;
					
					Group.MainThreads.push_back(i);
					Group.ThreadList.push_back(ThreadList);

					Count++;
				}

			}

		}

	/*	if(tId == 0){

			printf("%i groups:\n", GroupList.size());

			for(int i=0; i<GroupList.size(); i++){

				printf(" group %i:\n", i);

				for(int c=0; c<GroupList[i].MainThreads.size(); c++){

					printf("   main thread = %i\n", GroupList[i].MainThreads[c]);
					printf("   sub threads = ");

					if(GroupList[i].ThreadList[c].size() > 0)
					for(int c2=0; c2<GroupList[i].ThreadList[c].size(); c2++)
						printf("%i%s", GroupList[i].ThreadList[c][c2], (c2==GroupList[i].ThreadList[c].size()-1) ? "\n" : ", ");

					if(GroupList[i].ThreadList[c].size() == 0)
						printf("(none)\n");

				}

			}

			fflush(stdout);
		}	//*/

	}

	int GroupThreadRole(TransferGroup& Group){

		int FoundInd = -1;

		for(int i=0; i<Group.MainThreads.size(); i++){

			if(Group.MainThreads[i] == tId){		//Is main thread

				FoundInd = i;

			}else{

				for(int c=0; c<Group.ThreadList[i].size(); c++){

					if(Group.ThreadList[i][c] == tId){		//Is sub thread

						FoundInd = i;

					}

				}

			}

		}

		return FoundInd;
	}

	int RegionThreadList(GridRegion* Rgn, int* tList){

		int Ind = 0;
	
		for(int i=0; i<nThreads; i++){

			if(Regions[i].TestOverlap(Rgn)){
			
				tList[Ind++] = i;
			}

		}
	
		return Ind;		//Count
	}
	
public:

	//Constructors

	GridMPI(MPI_Comm MPIComm, int MPITag, int SizeX,  int SizeY, int SizeZ, int VectorSize){								//! Initialise a 3D grid of size SizeX x SizeY x SizeZ with vector entries of dimension VectorSize
		Init(MPIComm, MPITag, SizeX, SizeY, SizeZ, VectorSize);
	}

	GridMPI(MPI_Comm MPIComm, int MPITag, long long SizeX,  long long SizeY, long long SizeZ, long long VectorSize){		//! Initialise a 3D grid of size SizeX x SizeY x SizeZ with vector entries of dimension VectorSize
		Init(MPIComm, MPITag, SizeX, SizeY, SizeZ, VectorSize);
	}

	GridMPI(MPI_Comm MPIComm, int MPITag, int SizeX,  int SizeY, int SizeZ){							//! Initialise a 3D grid of size SizeX x SizeY x SizeZ with scalar entries
		Init(MPIComm, MPITag, SizeX, SizeY, SizeZ, 1);
	}

	GridMPI(MPI_Comm MPIComm, int MPITag, long long SizeX,  long long SizeY, long long SizeZ){		//! Initialise a 3D grid of size SizeX x SizeY x SizeZ with scalar entries
		Init(MPIComm, MPITag, SizeX, SizeY, SizeZ, 1);
	}

	GridMPI(MPI_Comm MPIComm, int MPITag, int SizeX,  int SizeY){					//! Initialise a 2D grid of size SizeX x SizeY with scalar entries
		Init(MPIComm, MPITag, SizeX, SizeY, 1, 1);
	}

	GridMPI(MPI_Comm MPIComm, int MPITag, long long SizeX,  long long SizeY){		//! Initialise a 2D grid of size SizeX x SizeY with scalar entries
		Init(MPIComm, MPITag, SizeX, SizeY, 1, 1);
	}

	GridMPI(MPI_Comm MPIComm, int MPITag, int SizeX){								//! Initialise an array of size SizeX with scalar entries
		Init(MPIComm, MPITag, SizeX, 1, 1, 1);
	}

	GridMPI(MPI_Comm MPIComm, int MPITag, long long SizeX){							//! Initialise an array of size SizeX with scalar entries
		Init(MPIComm, MPITag, SizeX, 1, 1, 1);
	}

	GridMPI(MPI_Comm MPIComm, int MPITag, Grid<T>& GridLocal, GridRegion RegionLocal){

		FreeArr = false;

	}

	GridMPI(MPI_Comm MPIComm, int MPITag, T* DataArr, bool ShallowCopy, bool FreeShallowCopy, long long SizeX,  long long SizeY, long long SizeZ, long long VectorSize){
	
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

	//*/

	//Destructor

	~GridMPI(){

		if(FreeArr && Data!=0){
			delete[] Data;
			Data = 0;
		}

		if(Regions!=0){
			delete[] Regions;
			Regions = 0;
		}

		if(GridLocal!=0){
			delete GridLocal;
			GridLocal = 0;
		}
	}

	//Buffer size

	void SetMaxBufferSize(long long BufSizeBytes){
	
		BufSizeMax = min((long long)INT_MAX, BufSizeBytes);
	}

	//Coordinates

	void CoordLocalToGlobal(int Pos[3]){

		Pos[0] += Regions[tId].x0;
		Pos[1] += Regions[tId].y0;
		Pos[2] += Regions[tId].z0;

	}

	void CoordLocalToGlobal(int* x, int* y, int* z){

		*x += Regions[tId].x0;
		*y += Regions[tId].y0;
		*z += Regions[tId].z0;

	}

	void CoordGlobalToLocal(int Pos[3]){

		Pos[0] -= Regions[tId].x0;
		Pos[1] -= Regions[tId].y0;
		Pos[2] -= Regions[tId].z0;

	}

	void CoordGlobalToLocal(int* x, int* y, int* z){

		*x -= Regions[tId].x0;
		*y -= Regions[tId].y0;
		*z -= Regions[tId].z0;

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

		long long CountTot;
		long long Count = 0;

		GridRegion LocalRgn = Regions[tId].Overlap(Region);

		if(!LocalRgn.IsZero()){

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			Count = GridLocal->CountOccurences(Value, LocalRgn);

		}

		MPI_Allreduce(&Count, &CountTot, 1, MPI_LONG_LONG, MPI_SUM, Comm);

		return CountTot;
	}

	long long CountOccurences(T Value[]){		//! Counts instances of vectors

		return CountOccurences(Value, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	long long CountOccurences(T* ValueList, int nValues, long long* CountOut, GridRegion Region){

		long long* Count = new long long[nValues];
		
		GridRegion LocalRgn = Regions[tId].Overlap(Region);

		if(!LocalRgn.IsZero()){

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			GridLocal->CountOccurences(ValueList, nValues, Count, LocalRgn);

		}

		MPI_Allreduce(Count, CountOut, nValues, MPI_LONG_LONG, MPI_SUM, Comm);

		delete[] Count;

		long long CountTot = 0;

		for(int i=0; i<nValues; i++)
			CountTot += CountOut[i];

		return CountTot;
	}

	long long CountOccurences(T* ValueList, int nValues, long long* CountOut){

		return CountOccurences(ValueList, nValues, CountOut, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	long long CountOccurences(T Value, GridRegion Region){		//! Count instances of a value within a given region

		long long CountTot;
		long long Count = 0;

		GridRegion LocalRgn = Regions[tId].Overlap(Region);

		if(!LocalRgn.IsZero()){

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			Count = GridLocal->CountOccurences(Value, LocalRgn);

		}

		MPI_Allreduce(&Count, &CountTot, 1, MPI_LONG_LONG, MPI_SUM, Comm);

		return CountTot;
	}

	long long CountOccurences(T Value){		//! Count instances of a value

		return CountOccurences(Value, GridRegion(0LL,SizeX(),0LL,SizeY(),0LL,SizeZ()));

	}

	void CountOccurences(long long* Count, T Value, GridRegion* RegionSet, int n){

		long long* CountLocal = new long long[n];
	
		for(int i=0; i<n; i++){

			CountLocal[i] = 0;
			
			GridRegion LocalRgn = RegionSet[i].Overlap(Regions[i]);

			if(LocalRgn.IsZero())
				continue;

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			CountLocal[i] = GridLocal->CountOccurences(Value, LocalRgn);

		}

		MPI_Allreduce(CountLocal, Count, n, MPI_LONG_LONG, MPI_SUM, Comm);

		delete[] CountLocal;
	
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

		GridRegion LocalRgn = Regions[tId].Overlap(Region);

		if(!LocalRgn.IsZero()){

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			GridLocal->SetAll(Value, LocalRgn);

		}

	}

	void SetAll(T Value[], GridRegion Region){				//! Fill region with given value vector

		GridRegion LocalRgn = Regions[tId].Overlap(Region);

		if(!LocalRgn.IsZero()){

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			GridLocal->SetAll(Value, LocalRgn);

		}

	}

	void SetAll(T Value){									//! Fill grid with given value

		GridLocal->SetAll(Value);

	}

	void SetAll(T Value[]){									//! Fill grid with given value

		GridLocal->SetAll(Value);

	}

	void ReplaceValue(T ValueFind, T ValueReplace, GridRegion Region){		//! Replace value within given region

		GridRegion LocalRgn = Regions[tId].Overlap(Region);

		if(!LocalRgn.IsZero()){

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			GridLocal->ReplaceValue(ValueFind, ValueReplace, LocalRgn);

		}

	}

	void ReplaceValue(T ValueFind[], T ValueReplace[], GridRegion Region){		//! Replace value vector within given region

		GridRegion LocalRgn = Regions[tId].Overlap(Region);

		if(!LocalRgn.IsZero()){

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			GridLocal->ReplaceValue(ValueFind, ValueReplace, LocalRgn);

		}
	}

	void ReplaceValue(T ValueFind, T ValueReplace){				//! Replace ValueFind with ValueReplace

		GridLocal->ReplaceValue(ValueFind, ValueReplace);

	}

	void ReplaceValue(T ValueFind[], T ValueReplace[]){			//! Replace ValueFind vector with ValueReplace vector

		GridLocal->ReplaceValue(ValueFind, ValueReplace);

	}

	void ScaleValues(T Scale, GridRegion Region){	//! Multiply values by Scale within a given region

		GridRegion LocalRgn = Regions[tId].Overlap(Region);

		if(!LocalRgn.IsZero()){

			LocalRgn = LocalRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

			GridLocal->ScaleValues(Scale, LocalRgn);

		}
	}

	void ScaleValues(T Scale){						//! Multiply values by Scale

		GridLocal->ScaleValues(Scale);

	}


	//Grid operations

	void ObtainSubGrid(Grid<T>& GridOut, GridRegion Region){

		ObtainSubGrid(GridOut, Region, GridPosition(0,0,0));
	}

	void ObtainSubGrid(Grid<T>& GridOut, GridRegion Region, GridPosition OffsetOut){
	
		int Boundary[] = {0, 0, 0, 0, 0, 0};

		ObtainSubGrid(GridOut, Region, OffsetOut, Boundary);
	
	}

	void ObtainSubGrid(Grid<T>& GridOut, GridRegion Region, GridPosition OffsetOut, int Boundary[6]){

		using namespace GridMPIDefinitions;

		GridRegion* SubRegions = new GridRegion[nThreads];
		int* BoundaryList = new int[nThreads*6];

		MPI_Allgather(&Region, sizeof(GridRegion), MPI_BYTE, SubRegions, sizeof(GridRegion), MPI_BYTE, Comm);
		MPI_Allgather(Boundary, 6, MPI_INT, BoundaryList, 6, MPI_INT, Comm);

		vector<TransferGroup> GroupList;

		GroupThreadOps(GroupList, SubRegions, BoundaryList);

		//Buffer

		long long nNodesTot = Sz[1]*Sz[2]*Sz[3];

		long long BufSize = min(BufSizeMax, (long long)(Sz[0]*nNodesTot*sizeof(T)));

		long long nNodesBuf = BufSize/(Sz[0]*sizeof(T));

		char* Buf = new char[BufSize];

		//Distribute

		MPI_Request* Rq = new MPI_Request[nThreads];
		char** tBuf = new char*[nThreads];
		MPI_Status Stat;
		MPI_Request SendRq;

		for(int GroupInd=0; GroupInd < GroupList.size(); GroupInd++){

			int ThreadGroup = GroupThreadRole(GroupList[GroupInd]);

			MPI_Comm GroupComm;
			MPI_Comm_split(Comm, ThreadGroup, 0, &GroupComm);

			if(ThreadGroup == -1){		//Thread not needed in this set
				MPI_Comm_free(&GroupComm);
				continue;
			}

			int n = GroupList[GroupInd].MainThreads[ThreadGroup];

		//	printf("Thread %i participating in group %i, main thread %i\n", tId, GroupInd, n);
		//	fflush(stdout);

		//	for(int n=0;  n<nThreads; n++ )			//Each thread
			for(int iv=0; iv<27     ; iv++){		//Vector components for boundary layers

				GridRegion Rgn = SubRegions[n];

				//Incorporate boundary layers

				int OffsetBoundary[3] = {0, 0, 0};

				int* BSz = &BoundaryList[n*6];	//Boundary sizes {+X, -X, +Y, -Y, +Z, -Z}

				if(iv == 0){		//Main grid
			
					if(Rgn.x1 < Sz[1]) Rgn.x1 += BSz[0];			//Add boundary layers if no loop required
					if(Rgn.x0 > 0    ) Rgn.x0 -= BSz[1];
					if(Rgn.y1 < Sz[2]) Rgn.y1 += BSz[2];
					if(Rgn.y0 > 0    ) Rgn.y0 -= BSz[3];
					if(Rgn.z1 < Sz[3]) Rgn.z1 += BSz[4];
					if(Rgn.z0 > 0    ) Rgn.z0 -= BSz[5];
				
					OffsetBoundary[0] = (Rgn.x0 == 0) ? BSz[1] : 0;
					OffsetBoundary[1] = (Rgn.y0 == 0) ? BSz[3] : 0;
					OffsetBoundary[2] = (Rgn.z0 == 0) ? BSz[5] : 0;

				}else{

					int v[3];
					memcpy(v, Q26Basis[iv-1], sizeof(int)*3);

					bool RqXp = (v[0] > 0 && Rgn.x1 == Sz[1] && BSz[0] > 0);
					bool RqXn = (v[0] < 0 && Rgn.x0 == 0     && BSz[1] > 0);
					bool RqYp = (v[1] > 0 && Rgn.y1 == Sz[2] && BSz[2] > 0);
					bool RqYn = (v[1] < 0 && Rgn.y0 == 0     && BSz[3] > 0);
					bool RqZp = (v[2] > 0 && Rgn.z1 == Sz[3] && BSz[4] > 0);
					bool RqZn = (v[2] < 0 && Rgn.z0 == 0     && BSz[5] > 0);

					if(!RqXp && !RqXn && !RqYp && !RqYn && !RqZp && !RqZn)	//Layer not over boundary
						continue;
			
					//Boundary region

					OffsetBoundary[0] = BSz[1];
					OffsetBoundary[1] = BSz[3];
					OffsetBoundary[2] = BSz[5];

					if(v[0] > 0){
						Rgn.SetRangeX(Rgn.x1, Rgn.x1 + BSz[0]*v[0]);
						OffsetBoundary[0] = BSz[1] + (int)SubRegions[n].SizeX();
					}
					if(v[1] > 0){
						Rgn.SetRangeY(Rgn.y1, Rgn.y1 + BSz[2]*v[1]);
						OffsetBoundary[1] = BSz[3] + (int)SubRegions[n].SizeY();
					}
					if(v[2] > 0){
						Rgn.SetRangeZ(Rgn.z1, Rgn.z1 + BSz[4]*v[2]);
						OffsetBoundary[2] = BSz[5] + (int)SubRegions[n].SizeZ();
					}

					if(v[0] < 0){
						Rgn.SetRangeX(Rgn.x0 + BSz[1]*v[0], Rgn.x0);
						OffsetBoundary[0] = 0;
					}
					if(v[1] < 0){
						Rgn.SetRangeY(Rgn.y0 + BSz[3]*v[1], Rgn.y0);
						OffsetBoundary[1] = 0;
					}
					if(v[2] < 0){
						Rgn.SetRangeZ(Rgn.z0 + BSz[5]*v[2], Rgn.z0);
						OffsetBoundary[2] = 0;
					}

					//Loop boundaries

					long long Translation[3] = {0, 0, 0};

					if(Rgn.x0 < 0) Translation[0] += Sz[1];		//Loop any coordinates over boundaries
					if(Rgn.y0 < 0) Translation[1] += Sz[2];
					if(Rgn.z0 < 0) Translation[2] += Sz[3];

					if(Rgn.x1 > Sz[1]) Translation[0] -= Sz[1];
					if(Rgn.y1 > Sz[2]) Translation[1] -= Sz[2];
					if(Rgn.z1 > Sz[3]) Translation[2] -= Sz[3];

					Rgn = Rgn.Translate(Translation[0], Translation[1], Translation[2]);

				}

				long long nNodesRgn = Rgn.CountCells();

				if(nNodesRgn == 0)
					continue;

				long long SzRgn[4];

				SzRgn[0] = Sz[0];
				SzRgn[1] = Rgn.SizeX();
				SzRgn[2] = Rgn.SizeY();
				SzRgn[3] = Rgn.SizeZ();

				long long BufSz[3];

				BufSz[2] = min(SzRgn[3], nNodesBuf / (SzRgn[1]*SzRgn[2]) );						//Number of full planes in Z
				BufSz[1] = BufSz[2] > 0 ? SzRgn[2] : min(SzRgn[2], nNodesBuf/SzRgn[1] );		//Number of full lines in Y
				BufSz[0] = BufSz[1] > 0 ? SzRgn[1] : min(SzRgn[1], nNodesBuf );					//Size in X

				GridRegion BufRgn(Rgn.x0, Rgn.x0 + BufSz[0], Rgn.y0, Rgn.y0 + max(1LL, BufSz[1]), Rgn.z0, Rgn.z0 + max(1LL, BufSz[2]));	//Initial buffer region

				for(long long Count=0; Count<nNodesRgn; ){

					//Read data to write thread

					long long RecvBufInd = 0;

					for(int i=0; i<nThreads; i++){
			
						if(!Regions[i].TestOverlap(BufRgn))
							continue;

						GridRegion SendRgn = Regions[i].Overlap(BufRgn);

						if(tId == n && i != n){
				
							MPI_Irecv(&Buf[RecvBufInd], (int)(SendRgn.CountCells()*Sz[0]*sizeof(T)), MPI_CHAR, i, Tag, Comm, &Rq[i]);

							tBuf[i] = &Buf[RecvBufInd];

							RecvBufInd += SendRgn.CountCells()*Sz[0]*sizeof(T);

						}
				
						if(tId == i){

							char* SendBuf = i==n ? &Buf[RecvBufInd] : Buf;
				
							Grid<T> GrdBuf((T*)SendBuf, true, false, SendRgn.SizeX(), SendRgn.SizeY(), SendRgn.SizeZ(), Sz[0]);

							GridRegion RgnLocal = SendRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

							GrdBuf.ReadSubGrid(*GridLocal, RgnLocal);

							if(i != n){
						
								MPI_Isend(Buf, (int)(SendRgn.CountCells()*Sz[0]*sizeof(T)), MPI_CHAR, n, Tag, Comm, &SendRq);

							}else{

								tBuf[i] = &Buf[RecvBufInd];
								RecvBufInd += SendRgn.CountCells()*Sz[0]*sizeof(T);

							}
			
						}

					}

					//Wait for data transfer

					for(int i=0; i<nThreads; i++){

						if(i == n || !Regions[i].TestOverlap(BufRgn))
							continue;

						if(tId == n)
							MPI_Wait(&Rq[i], &Stat);		//Wait recv from neighbour
				
						if(tId == i)
							MPI_Wait(&SendRq, &Stat);		//Wait send to main thread

					}

					//Write data to local grid

					if(tId == n){

						for(int i=0; i<nThreads; i++){

							if(!Regions[i].TestOverlap(BufRgn))
								continue;

							GridRegion SendRgn = Regions[i].Overlap(BufRgn);

							Grid<T> GrdBuf((T*)tBuf[i], true, false, SendRgn.SizeX(), SendRgn.SizeY(), SendRgn.SizeZ(), Sz[0]);
						
							GridPosition Pos(SendRgn.x0 + OffsetBoundary[0] + OffsetOut.x - Rgn.x0,
											 SendRgn.y0 + OffsetBoundary[1] + OffsetOut.y - Rgn.y0,
											 SendRgn.z0 + OffsetBoundary[2] + OffsetOut.z - Rgn.z0);

							GridOut.ReadGrid(GrdBuf, Pos);

						}
					
					}

					MPI_Barrier(GroupComm);

					Count += BufRgn.CountCells();

					//Next buffer region

					if(BufRgn.x1-Rgn.x0 < SzRgn[1]){			//Shift in x
						BufRgn.x0 = BufRgn.x1;
						BufRgn.x1 = min(BufRgn.x1 + BufSz[0], Rgn.x0+SzRgn[1]);
					}else if(BufRgn.y1-Rgn.y0 < SzRgn[2]){		//Shift in y
						BufRgn.x0 = Rgn.x0;
						BufRgn.x1 = Rgn.x0 + BufSz[0];
						BufRgn.y0 = BufRgn.y1;
						BufRgn.y1 = min(BufRgn.y1 + max(1LL, BufSz[1]), Rgn.y0+SzRgn[2]);
					}else if(BufRgn.z1-Rgn.z0 < SzRgn[3]){		//Shift in z
						BufRgn.x0 = Rgn.x0;
						BufRgn.x1 = Rgn.x0 + BufSz[0];
						BufRgn.y0 = Rgn.y0;
						BufRgn.y1 = Rgn.y0 + max(1LL, BufSz[1]);
						BufRgn.z0 = BufRgn.z1;
						BufRgn.z1 = min(BufRgn.z1 + max(1LL, BufSz[2]), Rgn.z0+SzRgn[3]);
					}

				}

			}

			MPI_Comm_free(&GroupComm);
		}

		delete[] Rq;
		delete[] tBuf;
		delete[] Buf;
		delete[] SubRegions;
		delete[] BoundaryList;

	}
	
	void ReadSubGrid(Grid<T>& GridRef, GridRegion Region){		//Read a given Region of GridRef

		bool Flip[3] = {false, false, false};

		ReadSubGrid(GridRef, Region, GridPosition(0, 0, 0), Flip);

	}

	void ReadSubGrid(Grid<T>& GridRef, GridRegion Region, GridPosition ReadOffset){

		bool Flip[3] = {false, false, false};

		ReadSubGrid(GridRef, Region, GridPosition(0, 0, 0), Flip);
	}

	void ReadSubGrid(Grid<T>& GridRef, GridRegion Region, GridPosition ReadOffset, bool Flip[3]){		//Read from GridRef into main grid

		GridRegion* SubRegions = new GridRegion[nThreads];

		MPI_Allgather(&Region, sizeof(GridRegion), MPI_BYTE, SubRegions, sizeof(GridRegion), MPI_BYTE, Comm);

		vector<TransferGroup> GroupList;

		GroupThreadOps(GroupList, SubRegions, 0);

		//Buffer

		long long nNodesTot = Sz[1]*Sz[2]*Sz[3];

		long long BufSize = min(BufSizeMax, (long long)(Sz[0]*nNodesTot*sizeof(T)));

		long long nNodesBuf = BufSize/(Sz[0]*sizeof(T));

		char* Buf = new char[BufSize];

		//Distribute

		MPI_Request* Rq = new MPI_Request[nThreads];
		char** tBuf = new char*[nThreads];
		MPI_Status Stat;
		MPI_Request RecvRq;

		for(int GroupInd=0; GroupInd < GroupList.size(); GroupInd++){

			//Group

			int ThreadGroup = GroupThreadRole(GroupList[GroupInd]);

			MPI_Comm GroupComm;
			MPI_Comm_split(Comm, ThreadGroup, 0, &GroupComm);

			if(ThreadGroup == -1){		//Thread not needed in this set
				MPI_Comm_free(&GroupComm);
				continue;
			}

			int n = GroupList[GroupInd].MainThreads[ThreadGroup];

			//Read

			GridRegion& Rgn = SubRegions[n];

			long long nNodesRgn = Rgn.CountCells();

			if(nNodesRgn == 0)
				continue;

			long long SzRgn[4];

			SzRgn[0] = Sz[0];
			SzRgn[1] = Rgn.SizeX();
			SzRgn[2] = Rgn.SizeY();
			SzRgn[3] = Rgn.SizeZ();

			long long BufSz[3];

			BufSz[2] = min(SzRgn[3], nNodesBuf / (SzRgn[1]*SzRgn[2]) );						//Number of full planes in Z
			BufSz[1] = BufSz[2] > 0 ? SzRgn[2] : min(SzRgn[2], nNodesBuf/SzRgn[1] );		//Number of full lines in Y
			BufSz[0] = BufSz[1] > 0 ? SzRgn[1] : min(SzRgn[1], nNodesBuf );					//Size in X

			GridRegion BufRgn(Rgn.x0, Rgn.x0 + BufSz[0], Rgn.y0, Rgn.y0 + max(1LL, BufSz[1]), Rgn.z0, Rgn.z0 + max(1LL, BufSz[2]));	//Initial buffer region

			for(long long Count=0; Count<nNodesRgn; ){

				long long RecvBufInd = 0;

				for(int i=0; i<nThreads; i++){
			
					if(!Regions[i].TestOverlap(BufRgn))
						continue;

					GridRegion SendRgn = Regions[i].Overlap(BufRgn);

					if(tId == n){

						tBuf[i] = &Buf[RecvBufInd];

						RecvBufInd += SendRgn.CountCells()*Sz[0]*sizeof(T);
				
					}
				
					if(tId == i && i != n){

						MPI_Irecv(Buf, (int)(SendRgn.CountCells()*Sz[0]*sizeof(T)), MPI_CHAR, n, Tag, Comm, &RecvRq);
			
					}
			
				}

				//Read from local grid

				if(tId == n){

					//Find

					for(int i=0; i<nThreads; i++){

						if(!Regions[i].TestOverlap(BufRgn))
							continue;

						GridRegion SendRgn = Regions[i].Overlap(BufRgn);

						Grid<T> GrdBuf((T*)tBuf[i], true, false, SendRgn.SizeX(), SendRgn.SizeY(), SendRgn.SizeZ(), Sz[0]);

						GridRegion RgnLocal = SendRgn.Translate(-Region.x0, -Region.y0, -Region.z0);

						if(Flip[0]) RgnLocal.SetRangeX(Region.SizeX() - RgnLocal.x1, Region.SizeX() - RgnLocal.x0);
						if(Flip[1]) RgnLocal.SetRangeY(Region.SizeY() - RgnLocal.y1, Region.SizeY() - RgnLocal.y0);
						if(Flip[2]) RgnLocal.SetRangeZ(Region.SizeZ() - RgnLocal.z1, Region.SizeZ() - RgnLocal.z0);

						RgnLocal = RgnLocal.Translate(ReadOffset.x, ReadOffset.y, ReadOffset.z);
						
						GrdBuf.ReadSubGrid(GridRef, RgnLocal, GridPosition(0,0,0), Flip);

						if(i != n){

							MPI_Isend(tBuf[i], (int)(SendRgn.CountCells()*Sz[0]*sizeof(T)), MPI_CHAR, i, Tag, Comm, &Rq[i]);

						}

					}

				}

				//Wait for data transfer

				for(int i=0; i<nThreads; i++){

					if(!Regions[i].TestOverlap(BufRgn))
						continue;

					if(tId == n && i != n)
						MPI_Wait(&Rq[i], &Stat);		//Wait to complete sends
				
					if(tId == i){

						if(i != n)
							MPI_Wait(&RecvRq, &Stat);	//Wait to recv data

						GridRegion SendRgn = Regions[i].Overlap(BufRgn);

						char* RecvBuf = i==n ? tBuf[i] : Buf;

						Grid<T> GrdBuf((T*)RecvBuf, true, false, SendRgn.SizeX(), SendRgn.SizeY(), SendRgn.SizeZ(), Sz[0]);

						GridLocal->ReadGrid(GrdBuf, GridPosition(SendRgn.x0-Regions[tId].x0, SendRgn.y0-Regions[tId].y0, SendRgn.z0-Regions[tId].z0));

					}

				}

				MPI_Barrier(GroupComm);

				Count += BufRgn.CountCells();

				//Next buffer region

				if(BufRgn.x1-Rgn.x0 < SzRgn[1]){			//Shift in x
					BufRgn.x0 = BufRgn.x1;
					BufRgn.x1 = min(BufRgn.x1 + BufSz[0], Rgn.x0+SzRgn[1]);
				}else if(BufRgn.y1-Rgn.y0 < SzRgn[2]){		//Shift in y
					BufRgn.x0 = Rgn.x0;
					BufRgn.x1 = Rgn.x0 + BufSz[0];
					BufRgn.y0 = BufRgn.y1;
					BufRgn.y1 = min(BufRgn.y1 + max(1LL, BufSz[1]), Rgn.y0+SzRgn[2]);
				}else if(BufRgn.z1-Rgn.z0 < SzRgn[3]){		//Shift in z
					BufRgn.x0 = Rgn.x0;
					BufRgn.x1 = Rgn.x0 + BufSz[0];
					BufRgn.y0 = Rgn.y0;
					BufRgn.y1 = Rgn.y0 + max(1LL, BufSz[1]);
					BufRgn.z0 = BufRgn.z1;
					BufRgn.z1 = min(BufRgn.z1 + max(1LL, BufSz[2]), Rgn.z0+SzRgn[3]);
				}

			}

			MPI_Comm_free(&GroupComm);
		}

		delete[] Rq;
		delete[] tBuf;
		delete[] Buf;
		delete[] SubRegions;

	}

	void CopyRegion(GridRegion RegionFrom, GridPosition PosTo){

		bool Flip[3] = {false, false, false};

		CopyRegion(RegionFrom, PosTo, Flip);

	}

	void CopyRegion(GridRegion RegionFrom, GridPosition PosTo, bool Flip[3]){

		GridRegion RgnLocal = Regions[tId].Overlap(RegionFrom);

		GridRegion RgnSend = RgnLocal.Translate(-RegionFrom.x0, -RegionFrom.y0, -RegionFrom.z0);

		if(Flip[0]) RgnSend.SetRangeX(RegionFrom.SizeX() - RgnSend.x1, RegionFrom.SizeX() - RgnSend.x0);
		if(Flip[1]) RgnSend.SetRangeY(RegionFrom.SizeY() - RgnSend.y1, RegionFrom.SizeY() - RgnSend.y0);
		if(Flip[2]) RgnSend.SetRangeZ(RegionFrom.SizeZ() - RgnSend.z1, RegionFrom.SizeZ() - RgnSend.z0);

		RgnSend = RgnSend.Translate(PosTo.x, PosTo.y, PosTo.z);

		GridPosition Offset(RgnLocal.x0 - Regions[tId].x0, RgnLocal.y0 - Regions[tId].y0, RgnLocal.z0 - Regions[tId].z0);

		ReadSubGrid(*GridLocal, RgnSend, Offset, Flip);

	}
	

	//File data

	int ReadFromFile(char* FileName, bool Binary, int ReadThread){

		return ReadFromFile(FileName, Binary, GridRegion(0LL, Sz[1], 0LL, Sz[2], 0LL, Sz[3]), ReadThread);
	}

	int ReadFromFile(char* FileName, bool Binary, GridRegion Region, int ReadThread){

		return ReadFromFile(FileName, Binary, ObtainDataType(typeid(T)), Region, ReadThread);
	}

	int ReadFromFile(char* FileName, bool Binary, DataType FileDataType, GridRegion Region, int ReadThread){		// Read dataset from file

		DataFileReader ReadData;

		//Open file

		bool ret;

		if(tId == ReadThread){

			ret = ReadData.Open(FileName, true);

		}

		MPI_Bcast(&ret, 1, MPI_CHAR, ReadThread, Comm);

		if(!ret)
			return 1;	//File could not be opened

		long long SzIn[4];

		SzIn[0] = Sz[0];
		SzIn[1] = Region.SizeX();
		SzIn[2] = Region.SizeY();
		SzIn[3] = Region.SizeZ();

		//Buffer

		long long nNodesTot = SzIn[1]*SzIn[2]*SzIn[3];

		long long BufSize = min(BufSizeMax, (long long)(Sz[0]*nNodesTot*sizeof(T)));

		char* Buf = new char[BufSize];
		
		long long nNodesBuf = BufSize/(Sz[0]*sizeof(T));

		long long BufSz[3];

		BufSz[2] = min(SzIn[3], nNodesBuf / (SzIn[1]*SzIn[2]) );						//Number of full planes in Z
		BufSz[1] = BufSz[2] > 0 ? SzIn[2] : min(SzIn[2], nNodesBuf/SzIn[1] );			//Number of full lines in Y
		BufSz[0] = BufSz[1] > 0 ? SzIn[1] : min(SzIn[1], nNodesBuf );					//Size in X

		GridRegion BufRgn(Region.x0+SzIn[1]-BufSz[0]          , Region.x0+SzIn[1],
						  Region.y0+SzIn[2]-max(1LL, BufSz[1]), Region.y0+SzIn[2],
						  Region.z0+SzIn[3]-max(1LL, BufSz[2]), Region.z0+SzIn[3]);		//Initial buffer region

		MPI_Request* Rq = new MPI_Request[nThreads];
		char** tBuf = new char*[nThreads];
		int* BuftList = new int[nThreads];
		MPI_Status Stat;
		MPI_Request RecvRq;

		int nBuftList = RegionThreadList(&BufRgn, BuftList);

		int ReadStatus = 0;

		for(long long Count=0; Count<nNodesTot; ){

			long long RecvBufInd = 0;

			for(int i=0; i<nBuftList; i++){

				int ListId = BuftList[i];

				GridRegion SendRgn = Regions[ListId].Overlap(BufRgn);

				if(tId == ReadThread){

					tBuf[i] = &Buf[RecvBufInd];

					RecvBufInd += SendRgn.CountCells()*Sz[0]*sizeof(T);
				
				}
				
				if(tId == ListId && ListId != ReadThread){

					MPI_Irecv(Buf, (int)(SendRgn.CountCells()*Sz[0]*sizeof(T)), MPI_CHAR, ReadThread, Tag, Comm, &RecvRq);
			
				}
			
			}

			//Read from file

			if(tId == ReadThread){

				for(long long z=BufRgn.z1-1; z>=BufRgn.z0; z--)
				for(long long y=BufRgn.y1-1; y>=BufRgn.y0; y--)
				for(long long x=BufRgn.x1-1; x>=BufRgn.x0;    ){

					//Find

					for(int i=0; i<nBuftList; i++){

						int ListId = BuftList[i];

						if(!Regions[ListId].Contains(x, y, z))
							continue;

						GridRegion SendRgn = Regions[ListId].Overlap(BufRgn);

						Grid<T> GrdBuf((T*)tBuf[i], true, false, SendRgn.SizeX(), SendRgn.SizeY(), SendRgn.SizeZ(), Sz[0]);

						long long dx = GrdBuf.SizeX();
						
						long long _y = y - SendRgn.y0;
						long long _z = z - SendRgn.z0;

						bool r = ReadData.ReadDataBack(GrdBuf.GetPtr(0LL, _y, _z, 0LL), dx*Sz[0], Binary, FileDataType, ObtainDataType(typeid(T)));

						if(!r)
							ReadStatus = 1;
						
						x -= dx;
					}
						
				}

				//Send

				for(int i=0; i<nBuftList; i++){

					int ListId = BuftList[i];
			
					if(ListId == ReadThread)
						continue;

					GridRegion SendRgn = Regions[ListId].Overlap(BufRgn);

					MPI_Isend(tBuf[i], (int)(SendRgn.CountCells()*Sz[0]*sizeof(T)), MPI_CHAR, ListId, Tag, Comm, &Rq[i]);
				
				}

			}

			//Wait for data transfer

			for(int i=0; i<nBuftList; i++){

				int ListId = BuftList[i];

				if(tId == ReadThread && ListId != ReadThread)
					MPI_Wait(&Rq[i], &Stat);		//Wait to complete sends
				
				if(tId == ListId){

					if(ListId != ReadThread)
						MPI_Wait(&RecvRq, &Stat);	//Wait to recv data

					GridRegion SendRgn = Regions[ListId].Overlap(BufRgn);

					char* RecvBuf = ListId==ReadThread ? tBuf[i] : Buf;

					Grid<T> GrdBuf((T*)RecvBuf, true, false, SendRgn.SizeX(), SendRgn.SizeY(), SendRgn.SizeZ(), Sz[0]);

					GridLocal->ReadGrid(GrdBuf, GridPosition(SendRgn.x0-Regions[tId].x0, SendRgn.y0-Regions[tId].y0, SendRgn.z0-Regions[tId].z0));

				}

			}

			MPI_Bcast(&ReadStatus, 1, MPI_INT, ReadThread, Comm);

			if(ReadStatus != 0)
				break;

			MPI_Barrier(Comm);
			
			Count += BufRgn.CountCells();

			//Next buffer region

			if(BufRgn.x0 > Region.x0){				//Shift in x
				BufRgn.x1 = BufRgn.x0;
				BufRgn.x0 = max(BufRgn.x0 - BufSz[0], Region.x0);
			}else if(BufRgn.y0 > Region.y0){		//Shift in y
				BufRgn.x0 = Region.x0+SzIn[1]-BufSz[0];
				BufRgn.x1 = Region.x0+SzIn[1];
				BufRgn.y1 = BufRgn.y0;
				BufRgn.y0 = max(BufRgn.y0 - max(1LL, BufSz[1]), Region.y0);
			}else if(BufRgn.z0 > Region.z0){		//Shift in z
				BufRgn.x0 = Region.x0+SzIn[1]-BufSz[0];
				BufRgn.x1 = Region.x0+SzIn[1];
				BufRgn.y0 = Region.y0+SzIn[2]-max(1LL, BufSz[1]);
				BufRgn.y1 = Region.y0+SzIn[2];
				BufRgn.z1 = BufRgn.z0;
				BufRgn.z0 = max(BufRgn.z0 - max(1LL, BufSz[2]), Region.z0);
			}

			nBuftList = RegionThreadList(&BufRgn, BuftList);	//List of overlapping thread regions

		}

		delete[] Rq;
		delete[] tBuf;
		delete[] BuftList;
		delete[] Buf;

		int retcode = 0;

		if(tId == ReadThread){
		
			if(ReadData.LastError() == 2) retcode = 2;		//Insufficient data
			if(ReadData.LastError() == 3) retcode = 3;		//Invalid data
			if(ReadData.LastError() == 4) retcode = 4;		//Memory allocation error

		}

		MPI_Bcast(&retcode, 1, MPI_INT, ReadThread, Comm);
		
		return retcode;
	}

	int WriteToFile(const char* FileName, bool Binary, bool VTKHeader, int WriteThread){

		return WriteToFile(FileName, Binary, VTKHeader, GridRegion(0LL, Sz[1], 0LL, Sz[2], 0LL, Sz[3]), WriteThread);
	}

	int WriteToFile(const char* FileName, bool Binary, bool VTKHeader, GridRegion Region, int WriteThread){	//Write dataset to file

		DataFileWriter WriteData;

		bool ret;

		if(tId == WriteThread){

			ret = WriteData.Open(FileName);

		}

		MPI_Bcast(&ret, 1, MPI_CHAR, WriteThread, Comm);

		if(!ret)
			return 1;	//File could not be opened

		long long SzOut[4];

		SzOut[0] = Sz[0];
		SzOut[1] = Region.SizeX();
		SzOut[2] = Region.SizeY();
		SzOut[3] = Region.SizeZ();

		//VTK header

		if(tId == WriteThread && VTKHeader){

			long long nValues = SzOut[1]*SzOut[2]*SzOut[3];

			char* HeaderBuf = new char[512];

			int l = _snprintf(HeaderBuf, 512,

				"# vtk DataFile Version 2.0\n"
				"%s\n"							//Header name
				"%s\n"							//ASCII or binary
				"DATASET STRUCTURED_POINTS\n"
				"DIMENSIONS %lli %lli %lli\n"	//Grid size
				"ORIGIN 0 0 0\n"
				"SPACING 1 1 1\n"
				"POINT_DATA %lli\n"				//Number of points
				"%s %s float\n"					//Vectors/Scalars; Dataset name
				"%s"							//Lookup table
				,							

				"DataGrid",
				Binary ? "binary" : "ASCII",
				SzOut[1], SzOut[2], SzOut[3],
				nValues,
				Sz[0]>1 ? "Vectors" : "Scalars",
				"Dataset",
				Sz[0]==1 ? "LOOKUP_TABLE default\n" : "" );

			WriteData.Write(HeaderBuf, l);

			delete[] HeaderBuf;

		}

		//Buffer

		long long nNodesTot = SzOut[1]*SzOut[2]*SzOut[3];

		long long BufSize = min(BufSizeMax, (long long)(Sz[0]*nNodesTot*sizeof(T)));

		char* Buf = new char[BufSize];
		
		long long nNodesBuf = BufSize/(Sz[0]*sizeof(T));

		long long BufSz[3];

		BufSz[2] = min(SzOut[3], nNodesBuf / (SzOut[1]*SzOut[2]) );						//Number of full planes in Z
		BufSz[1] = BufSz[2] > 0 ? SzOut[2] : min(SzOut[2], nNodesBuf/SzOut[1] );		//Number of full lines in Y
		BufSz[0] = BufSz[1] > 0 ? SzOut[1] : min(SzOut[1], nNodesBuf );					//Size in X

		GridRegion BufRgn(Region.x0, Region.x0 + BufSz[0], Region.y0, Region.y0 + max(1LL, BufSz[1]), Region.z0, Region.z0 + max(1LL, BufSz[2]));	//Initial buffer region

		MPI_Request* Rq = new MPI_Request[nThreads];
		char** tBuf = new char*[nThreads];
		int* BuftList = new int[nThreads];
		MPI_Status Stat;
		MPI_Request SendRq;

		int nBuftList = RegionThreadList(&BufRgn, BuftList);

		for(long long Count=0; Count<nNodesTot; ){

			//Read data to write thread

			long long RecvBufInd = 0;

			for(int i=0; i<nBuftList; i++){

				int ListId = BuftList[i];

				GridRegion SendRgn = Regions[ListId].Overlap(BufRgn);

				if(tId == WriteThread && ListId != WriteThread){
				
					MPI_Irecv(&Buf[RecvBufInd], (int)(SendRgn.CountCells()*Sz[0]*sizeof(T)), MPI_CHAR, ListId, Tag, Comm, &Rq[i]);

					tBuf[i] = &Buf[RecvBufInd];

					RecvBufInd += SendRgn.CountCells()*Sz[0]*sizeof(T);
				
				}
				
				if(tId == ListId){

					char* SendBuf = ListId==WriteThread ? &Buf[RecvBufInd] : Buf;
				
					Grid<T> GrdBuf((T*)SendBuf, true, false, SendRgn.SizeX(), SendRgn.SizeY(), SendRgn.SizeZ(), Sz[0]);

					GridRegion RgnLocal = SendRgn.Translate(-Regions[tId].x0, -Regions[tId].y0, -Regions[tId].z0);

					GrdBuf.ReadSubGrid(*GridLocal, RgnLocal);
			
					if(ListId != WriteThread){
						
						MPI_Isend(Buf, (int)(SendRgn.CountCells()*Sz[0]*sizeof(T)), MPI_CHAR, WriteThread, Tag, Comm, &SendRq);

					}else{

						tBuf[i] = &Buf[RecvBufInd];
						RecvBufInd += SendRgn.CountCells()*Sz[0]*sizeof(T);

					}

				}
			
			}

			//Wait for data transfer

			for(int i=0; i<nBuftList; i++){

				int ListId = BuftList[i];

				if(ListId == WriteThread)
					continue;

				if(tId == WriteThread)
					MPI_Wait(&Rq[i], &Stat);		//Wait recv from neighbour
				
				if(tId == ListId)
					MPI_Wait(&SendRq, &Stat);		//Wait send to main thread

			}

			//Write data to file

			if(tId == WriteThread){

				for(long long z=BufRgn.z0; z<BufRgn.z1; z++)
				for(long long y=BufRgn.y0; y<BufRgn.y1; y++)
				for(long long x=BufRgn.x0; x<BufRgn.x1;    ){

					//Find

					for(int i=0; i<nBuftList; i++){

						int ListId = BuftList[i];

						if(!Regions[ListId].Contains(x, y, z))
							continue;

						GridRegion SendRgn = Regions[ListId].Overlap(BufRgn);

						Grid<T> GrdBuf((T*)tBuf[i], true, false, SendRgn.SizeX(), SendRgn.SizeY(), SendRgn.SizeZ(), Sz[0]);

						long long dx = GrdBuf.SizeX();
						
						long long _y = y - SendRgn.y0;
						long long _z = z - SendRgn.z0;

						if(Binary){

							WriteData.Write(GrdBuf.GetPtr(0LL, _y, _z, 0LL), sizeof(T)*Sz[0]*dx);

						}else{
						
							for(long long c=0; c<dx; c++)
							for(long long v=0; v<Sz[0]; v++){
								WriteData.Write(GrdBuf.Get(c, _y, _z, v));
								WriteData.Write(' ');
							}
						
						}
						
						x += dx;
					}
						
				}

			}

			MPI_Barrier(Comm);

			Count += BufRgn.CountCells();

			//Next buffer region

			if(BufRgn.x1-Region.x0 < SzOut[1]){			//Shift in x
				BufRgn.x0 = BufRgn.x1;
				BufRgn.x1 = min(BufRgn.x1 + BufSz[0], Region.x0+SzOut[1]);
			}else if(BufRgn.y1-Region.y0 < SzOut[2]){	//Shift in y
				BufRgn.x0 = Region.x0;
				BufRgn.x1 = Region.x0 + BufSz[0];
				BufRgn.y0 = BufRgn.y1;
				BufRgn.y1 = min(BufRgn.y1 + max(1LL, BufSz[1]), Region.y0+SzOut[2]);
			}else if(BufRgn.z1-Region.z0 < SzOut[3]){	//Shift in z
				BufRgn.x0 = Region.x0;
				BufRgn.x1 = Region.x0 + BufSz[0];
				BufRgn.y0 = Region.y0;
				BufRgn.y1 = Region.y0 + max(1LL, BufSz[1]);
				BufRgn.z0 = BufRgn.z1;
				BufRgn.z1 = min(BufRgn.z1 + max(1LL, BufSz[2]), Region.z0+SzOut[3]);
			}

			nBuftList = RegionThreadList(&BufRgn, BuftList);

		}

		delete[] Rq;
		delete[] tBuf;
		delete[] BuftList;
		delete[] Buf;

		if(tId == WriteThread)
			WriteData.Close();
			
		return 0;
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
