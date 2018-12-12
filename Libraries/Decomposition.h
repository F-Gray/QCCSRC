#pragma once

#include <vector>

#include "Base.h"
#include "Grids.h"

#ifdef MPI_INCLUDED

#include "GridsMPI.h"

#endif

namespace DecompositionDefinitions{

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

enum NodeValueMode{

	List_Active

};

template <class T>
class Decomposition{

public:

	struct TransferRegion{

		GridRegion Region;						//Region of data to send local to partition

		vector<GridRegion> ExcludeRegions;		//Subregions of Region to exclude due to overlap with other TransferRegions

		vector<int> NeighList;					//Neighbour(s) to send/recv this region to
		vector<int> NeighIndexList;				//Index of corresponding region in neighbour(s) send/recv

		vector<vector<int> > VectorIndexList;	//List of vectors which map to this area for each neighbour


		//Functions

		bool IsExcluded(long long x, long long y, long long z){

			for(int i=0; i<ExcludeRegions.size(); i++){
				if(ExcludeRegions[i].Contains(x,y,z))
					return true;
			}

			return false;
		}

		long long CountCells(){					//The number of grid cells not in excluded areas

			long long Count = 0;

			for(long long z=Region.z0; z<Region.z1; z++)
			for(long long y=Region.y0; y<Region.y1; y++)
			for(long long x=Region.x0; x<Region.x1; x++){

				bool Excl = false;

				for(int i=0; i<ExcludeRegions.size(); i++){		//Remove contributions from excluded areas

					if(ExcludeRegions[i].Contains(x,y,z)){
						Excl = true;
						break;
					}
				}

				if(!Excl)
					Count++;
			}

			return Count;
		}

		void Print(){

			printf("Region (%s) with %lli cells\n", Region.ToStr(), CountCells());

			for(int i=0; i<ExcludeRegions.size(); i++){

				printf("  Excluding (%s)\n", ExcludeRegions[i].ToStr());

			}

			for(int i=0; i<NeighList.size(); i++){

				printf("  Neighbour %i [%i] by vectors: ", NeighList[i], NeighIndexList[i]);

				vector<int> &v = VectorIndexList[i];

				for(int c=0; c<v.size(); c++){

					printf(" %i", v[c]);

					if(c<v.size()-1)
						printf(", ");
				}

				if(i<NeighList.size()-1)
					printf("\n");

			}

			printf("\n");

		}
	};

	struct Neighbour{
	
		int NeighbourIndex;					//Neighbour
		int NeighbourListIndex;				//Index of corresponding region in Send/Recv list of neighbour

		GridRegion RegionGlobal;			//Overlap region in lattice coordinates
		GridRegion RegionLocalPartition;	//Overlap region in partition coordinates
		GridRegion RegionLocalNeighbour;	//Overlap region in neighbour partition coordinates

		long long NodeCount;

		vector<int> VectorIndexList;		//List of vectors which map to this area
	};

	struct Partition{		
	
		GridRegion Region;
	
		long long NodeCount;
	
		vector<Neighbour> NeighbourListRecv;	//Neighbours to receive into ghost nodes
		vector<Neighbour> NeighbourListSend;	//Neighbours to send outer nodes of partition to

		vector<TransferRegion> SendTransferRegions;		//Regions within partition to send to neighbours
		vector<TransferRegion> RecvTransferRegions;		//Regions within partition to receive from neighbours
	};

private:

	Grid<T>* _Lattice;		//Lattice
	long long arrSz[3];		//Lattice dimensions

	Partition DomainPartition;	//Partition referring to whole of domain

	int nDiv;				//Number of partitions total

	T* _NodeValue;			//Values of active nodes if WeightMode == true
	int _NodeValueCount;	//Number of active node values

	inline bool IsActiveNode(T Val){

		for(int i=0; i<_NodeValueCount; i++)
			if(Val == _NodeValue[i])
				return true;

		return false;
	}

	bool WeightMode;				// 0 = weight across all nodes; 1 = weight according to nodes with value _NodeValue

	Direction Dirs[3];		//Order of directions to decompose
	int nDirs;				//Number of directions to decompose

#ifdef MPI_INCLUDED

	GridMPI<T>* _LatticeMPI;	

	bool MPIMode;
	bool MPIGridMode;

	MPI_Comm MPIComm;

	int MPIGeoThread;
	int MPIThreadId;

#endif

	/////////////////////////////

	struct LatticeBlock{
		long long c0;			//Lower bound coordinate integer part
		double c0_;				//Lower bound coordinate fractional part
		long long c1;			//Upper bound coordinate integer part
		double c1_;				//Upper bound coordinate fractional part

		bool LowerRound;		//Round lower border down or up
		bool UpperRound;		//Round upper border down or up

		double RoundError[4];	//Error combinations of rounding borders
		bool RoundActive[4];	//Whether rounding combination is allowed
		int nRoundActive;		//Number of combinations allowed

		long long c0Rounded;	//Rounded lower border
		long long c1Rounded;	//Rounded upper border

		double n_d;				//Exact number of slices in plane
	};

	void InitLattice(Grid<T>* Lattice, bool Mode, T* NodeValue, int NodeValueCount){
	
		_Lattice = Lattice;

		arrSz[0] = Lattice->SizeX();
		arrSz[1] = Lattice->SizeY();
		arrSz[2] = Lattice->SizeZ();

		long long* CountBuf = new long long[NodeValueCount];

		DomainPartition.Region = GridRegion(0LL, arrSz[0], 0LL, arrSz[1], 0LL, arrSz[2]);
		DomainPartition.NodeCount = (!Mode) ? arrSz[0]*arrSz[1]*arrSz[2] : Lattice->CountOccurences(NodeValue, NodeValueCount, CountBuf);
	
		delete[] CountBuf;
	}

	void Init(Grid<T>* Lattice, int nPartitions, bool Mode, T* NodeValue, int NodeValueCount, Direction Dir1, Direction Dir2, Direction Dir3, int CountDirs){

		memset(this, 0, sizeof(Decomposition));

		if(CountDirs == 0){
			Success = false;
			return;
		}

		InitLattice(Lattice, Mode, NodeValue, NodeValueCount);

		nDiv = nPartitions;

		if(Mode){

			_NodeValue = new T[NodeValueCount];
		
			for(int i=0; i<NodeValueCount; i++)
				_NodeValue[i] = NodeValue[i];

			_NodeValueCount = NodeValueCount;

		}
		
		WeightMode = Mode;

		Dirs[0] = Dir1;
		Dirs[1] = Dir2;
		Dirs[2] = Dir3;

		nDirs = CountDirs;

		Partitions = new Partition[nDiv];

		Success = Decompose();
	}

#ifdef MPI_INCLUDED

	void InitLatticeMPI(GridMPI<T>* Lattice, bool Mode, T* NodeValue, int NodeValueCount){

		_LatticeMPI = Lattice;

		arrSz[0] = Lattice->SizeX();
		arrSz[1] = Lattice->SizeY();
		arrSz[2] = Lattice->SizeZ();

		long long* CountBuf = new long long[NodeValueCount];

		DomainPartition.Region = GridRegion(0LL, arrSz[0], 0LL, arrSz[1], 0LL, arrSz[2]);
		DomainPartition.NodeCount = (!Mode) ? arrSz[0]*arrSz[1]*arrSz[2] : Lattice->CountOccurences(NodeValue, NodeValueCount, CountBuf);

		delete[] CountBuf;
	}

	void InitMPI(Grid<T>* Lattice, int nPartitions, bool Mode, T* NodeValue, int NodeValueCount, Direction Dir1, Direction Dir2, Direction Dir3, int CountDirs, MPI_Comm Comm, int GeoThread){

		memset(this, 0, sizeof(Decomposition));

		if(CountDirs == 0){
			Success = false;
			return;
		}

		//Init lattice and transfer info with MPI

		int tId;
		MPI_Comm_rank(Comm, &tId);	//Get thread ID in local group

		if(tId == GeoThread){

			InitLattice(Lattice, Mode, NodeValue, NodeValueCount);

		}

		MPI_Bcast(this, sizeof(Decomposition), MPI_CHAR, GeoThread, Comm);

		//

		MPIGridMode = false;
		MPIMode = true;
		MPIComm = Comm;
		MPIGeoThread = GeoThread;
		MPIThreadId = tId;

		nDiv = nPartitions;

		if(Mode){

			_NodeValue = new T[NodeValueCount];
		
			for(int i=0; i<NodeValueCount; i++)
				_NodeValue[i] = NodeValue[i];

			_NodeValueCount = NodeValueCount;

		}

		WeightMode = Mode;

		Dirs[0] = Dir1;
		Dirs[1] = Dir2;
		Dirs[2] = Dir3;

		nDirs = CountDirs;

		Partitions = new Partition[nDiv];

		Success = Decompose();
	}

	void InitMPIGrid(GridMPI<T>* Lattice, int nPartitions, bool Mode, T* NodeValue, int NodeValueCount, Direction Dir1, Direction Dir2, Direction Dir3, int CountDirs, MPI_Comm Comm){
	
		memset(this, 0, sizeof(Decomposition));

		if(CountDirs == 0){
			Success = false;
			return;
		}

		InitLatticeMPI(Lattice, Mode, NodeValue, NodeValueCount);

		int tId;
		MPI_Comm_rank(Comm, &tId);	//Get thread ID in local group

		MPIGridMode = true;
		MPIMode = true;
		MPIComm = Comm;
		MPIThreadId = tId;

		nDiv = nPartitions;

		if(Mode){

			_NodeValue = new T[NodeValueCount];
		
			for(int i=0; i<NodeValueCount; i++)
				_NodeValue[i] = NodeValue[i];

			_NodeValueCount = NodeValueCount;

		}

		WeightMode = Mode;

		Dirs[0] = Dir1;
		Dirs[1] = Dir2;
		Dirs[2] = Dir3;

		nDirs = CountDirs;

		Partitions = new Partition[nDiv];

		Success = DecomposeMPI();
	
	}

#endif

public:

	//Data

	bool Success;		//Whether the decomposition completed succesfully

	Partition* Partitions;

	//Constructors

		//Single thread with Grid

	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1, T* NodeValue, int nValues){
		Init(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir_X, Dir_X, 1);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, T* NodeValue, int nValues){
		Init(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir2, Dir_X, 2);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3, T* NodeValue, int nValues){
		Init(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir2, Dir3, 3);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, const char* DirStr, T* NodeValue, int nValues){

		int nDecomp = 0;
		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		Init(Lattice, nPartitions, true, NodeValue, nValues, DirList[0], DirList[1], DirList[2], nDecomp);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1, T NodeValue){
		Init(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir_X, Dir_X, 1);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, T NodeValue){
		Init(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir2, Dir_X, 2);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3, T NodeValue){
		Init(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir2, Dir3, 3);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, const char* DirStr, T NodeValue){

		int nDecomp = 0;
		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		Init(Lattice, nPartitions, true, &NodeValue, 1, DirList[0], DirList[1], DirList[2], nDecomp);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1){
		Init(Lattice, nPartitions, false, 0, 0, Dir1, Dir_X, Dir_X, 1);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2){
		Init(Lattice, nPartitions, false, 0, 0, Dir1, Dir2, Dir_X, 2);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3){
		Init(Lattice, nPartitions, false, 0, 0, Dir1, Dir2, Dir3, 3);
	}
	Decomposition(Grid<T>* Lattice, int nPartitions, const char* DirStr){

		int nDecomp = 0;

		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		Init(Lattice, nPartitions, false, 0, 0, DirList[0], DirList[1], DirList[2], nDecomp);
	}

#ifdef MPI_INCLUDED

		//MPI with Grid on a single thread

	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1, T* NodeValue, int nValues){
		InitMPI(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir_X, Dir_X, 1, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, T* NodeValue, int nValues){
		InitMPI(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir2, Dir_X, 2, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3, T* NodeValue, int nValues){
		InitMPI(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir2, Dir3, 3, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, const char* DirStr, T* NodeValue, int nValues){

		int nDecomp = 0;
		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		InitMPI(Lattice, nPartitions, true, NodeValue, nValues, DirList[0], DirList[1], DirList[2], nDecomp, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1, T NodeValue){
		InitMPI(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir_X, Dir_X, 1, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, T NodeValue){
		InitMPI(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir2, Dir_X, 2, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3, T NodeValue){
		InitMPI(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir2, Dir3, 3, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, const char* DirStr, T NodeValue){

		int nDecomp = 0;
		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		InitMPI(Lattice, nPartitions, true, &NodeValue, 1, DirList[0], DirList[1], DirList[2], nDecomp, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1){
		InitMPI(Lattice, nPartitions, false, 0, 0, Dir1, Dir_X, Dir_X, 1, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2){
		InitMPI(Lattice, nPartitions, false, 0, 0, Dir1, Dir2, Dir_X, 2, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3){
		InitMPI(Lattice, nPartitions, false, 0, 0, Dir1, Dir2, Dir3, 3, Comm, GeoThread);
	}
	Decomposition(MPI_Comm Comm, int GeoThread, Grid<T>* Lattice, int nPartitions, const char* DirStr){

		int nDecomp = 0;
		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		InitMPI(Lattice, nPartitions, false, 0, 0, DirList[0], DirList[1], DirList[2], nDecomp, Comm, GeoThread);
	}

		//MPI with GridMPI distributed decomposition

	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1, T* NodeValue, int nValues){
		InitMPIGrid(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir_X, Dir_X, 1, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, T* NodeValue, int nValues){
		InitMPIGrid(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir2, Dir_X, 2, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3, T* NodeValue, int nValues){
		InitMPIGrid(Lattice, nPartitions, true, NodeValue, nValues, Dir1, Dir2, Dir3, 3, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, const char* DirStr, T* NodeValue, int nValues){

		int nDecomp = 0;
		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		InitMPIGrid(Lattice, nPartitions, true, NodeValue, nValues, DirList[0], DirList[1], DirList[2], nDecomp, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1, T NodeValue){
		InitMPIGrid(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir_X, Dir_X, 1, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, T NodeValue){
		InitMPIGrid(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir2, Dir_X, 2, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3, T NodeValue){
		InitMPIGrid(Lattice, nPartitions, true, &NodeValue, 1, Dir1, Dir2, Dir3, 3, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, const char* DirStr, T NodeValue){

		int nDecomp = 0;
		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		InitMPIGrid(Lattice, nPartitions, true, &NodeValue, 1, DirList[0], DirList[1], DirList[2], nDecomp, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1){
		InitMPIGrid(Lattice, nPartitions, false, 0, 0, Dir1, Dir_X, Dir_X, 1, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2){
		InitMPIGrid(Lattice, nPartitions, false, 0, 0, Dir1, Dir2, Dir_X, 2, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, Direction Dir1, Direction Dir2, Direction Dir3){
		InitMPIGrid(Lattice, nPartitions, false, 0, 0, Dir1, Dir2, Dir3, 3, Comm);
	}
	Decomposition(MPI_Comm Comm, GridMPI<T>* Lattice, int nPartitions, const char* DirStr){

		int nDecomp = 0;
		Direction DirList[3];

		InterpretDirStr(DirStr, DirList, &nDecomp);

		InitMPIGrid(Lattice, nPartitions, false, 0, 0, DirList[0], DirList[1], DirList[2], nDecomp, Comm);
	}

#endif

	//Static members

	static bool DecompRegion3D(GridRegion Rgn, int n, GridRegion* RegionsOut){

		long long Sz[3];

		Sz[0] = Rgn.SizeX();
		Sz[1] = Rgn.SizeY();
		Sz[2] = Rgn.SizeZ();
	
		//3D decomp

		double RatioXY = (double)Sz[0] / (double)Sz[1];
		double RatioXZ = (double)Sz[0] / (double)Sz[2];

		double N = max(1.0, pow( (double)n*RatioXY*RatioXZ , 1.0/3.0));

		int nDivX = min(n, (int)round(N));								//Number of partitions in X

		int nSubDivX = n / nDivX;
		int rmSubDivX = n % nDivX;

		double DivSzX = (double)Sz[0] / (double)n;

		int nXCount = 0;

		for(int i=0; i<nDivX; i++){

			int tCountDivX = nSubDivX + (i<rmSubDivX ? 1 : 0);
			int SubCount = (int)round(tCountDivX*DivSzX);

			for(int c=0; c<tCountDivX; c++){

				RegionsOut[nXCount+c].x0 = i==0 ?		 0 : RegionsOut[nXCount-1].x1;
				RegionsOut[nXCount+c].x1 = i==nDivX-1 ? Sz[0] : RegionsOut[nXCount+c].x0 + SubCount;

			}

			//Divide in Y

			double RatioYZ = (double)Sz[1] / (double)Sz[2];

			double NY = max(1.0, sqrt( (double)tCountDivX*RatioYZ ));

			int nDivY = min(tCountDivX, (int)round(NY));

			int* tCountDivY = new int[nDivY];

			int nSubDivY = tCountDivX / nDivY;
			int rmSubDivY = tCountDivX % nDivY;

			double DivSzY = (double)Sz[1] / (double)tCountDivX;

			int nYCount = 0;

			for(int i2=0; i2<nDivY; i2++){

				int tCountDivY = nSubDivY + (i2<rmSubDivY ? 1 : 0);
				int SubCountY = (int)round(tCountDivY*DivSzY);

				for(int c=0; c<tCountDivY; c++){

					RegionsOut[nXCount+nYCount+c].y0 = i2==0 ?			 0 : RegionsOut[nXCount+nYCount-1].y1;
					RegionsOut[nXCount+nYCount+c].y1 = i2==nDivY-1 ? Sz[1] : RegionsOut[nXCount+nYCount+c].y0 + SubCountY;

				}

				//Divide in Z

				int nZ = (int)Sz[2] / tCountDivY;
				int rmZ = (int)Sz[2] % tCountDivY;

				for(int i3=0; i3<tCountDivY; i3++){

					RegionsOut[nXCount+nYCount+i3].z0 = nZ*i3 + min(i3, rmZ);
					RegionsOut[nXCount+nYCount+i3].z1 = RegionsOut[nXCount+nYCount+i3].z0 + nZ + (i3<rmZ ? 1 : 0);

				}

				nYCount += tCountDivY;
			}

			nXCount += tCountDivX;
		}


		return true;
	}

	//Members

	~Decomposition(){

		if(_NodeValue != 0)
			delete[] _NodeValue;

		if(Partitions != 0)
			delete[] Partitions;
	}

	bool OutputDecompositionGeometry(char* FileName, bool Binary, bool VTKHeader){	

		Grid<int> Geo(arrSz[0], arrSz[1], arrSz[2]);

		for(int i=0; i<nDiv; i++){

			for(long long z=Partitions[i].Region.z0; z<Partitions[i].Region.z1; z++)
			for(long long y=Partitions[i].Region.y0; y<Partitions[i].Region.y1; y++)
			for(long long x=Partitions[i].Region.x0; x<Partitions[i].Region.x1; x++){
				
				Geo(x,y,z) = (WeightMode && IsActiveNode(_Lattice->Get(x,y,z))) ? 0 : (i+1);

			}

		}

		int r = Geo.WriteToFile(FileName, Binary, VTKHeader);

		return (r==0);
	}

	bool OutputPartitionNeighbourRegions(int n, char* FileName, bool Binary, bool VTKHeader){

		Partition Ptn = Partitions[n];

		//Width of ghost node areas

		int nBoundary = 0;

		for(int i=0; i<Ptn.NeighbourListRecv.size(); i++)
			nBoundary = (int)max(nBoundary, max(Ptn.Region.x0 - Ptn.NeighbourListRecv[i].RegionLocalPartition.x0, Ptn.NeighbourListRecv[i].RegionLocalPartition.x1 - Ptn.Region.x1));
		
		//Grid to output

		Grid<int> PartOut(Ptn.Region.SizeX() + 2*nBoundary, Ptn.Region.SizeY() + 2*nBoundary, Ptn.Region.SizeZ() + 2*nBoundary);

		PartOut.SetAll(n+1);

		//Recv areas

		for(int i=0; i<Ptn.NeighbourListRecv.size(); i++){

			GridRegion RgnSet = Ptn.NeighbourListRecv[i].RegionLocalPartition.Translate(nBoundary, nBoundary, nBoundary);

			PartOut.SetAll(Ptn.NeighbourListRecv[i].NeighbourIndex+1, RgnSet);

		}

		//Send areas

		for(int i=0; i<Ptn.NeighbourListSend.size(); i++){

			GridRegion RgnSet = Ptn.NeighbourListSend[i].RegionLocalPartition.Translate(nBoundary, nBoundary, nBoundary);

			PartOut.SetAll(Ptn.NeighbourListSend[i].NeighbourIndex+1, RgnSet);

		}

		int r = PartOut.WriteToFile(FileName, Binary, VTKHeader);

		return (r==0);
	}

	bool OutputPartitionSendTransferRegions(int n, char* FileName, bool Binary, bool VTKHeader){

		Partition Ptn = Partitions[n];

		//Grid to output

		Grid<int> PartOut(Ptn.Region.SizeX(), Ptn.Region.SizeY(), Ptn.Region.SizeZ());

		PartOut.SetAll(-1);

		//Send areas

		for(int i=0; i<Ptn.SendTransferRegions.size(); i++){

			PartOut.SetAll(i, Ptn.SendTransferRegions[i].Region);

		}

		int r = PartOut.WriteToFile(FileName, Binary, VTKHeader);

		return (r==0);
	}

	bool FindTransferRegionsTest(){

		vector<TransferRegion> SendTransferRegions;

		TransferRegion TrRgn[7];

		TrRgn[0].Region = GridRegion(0, 10, 0, 10, 0, 1);
		TrRgn[0].NeighList.push_back(1);

		int Vect[] = {1};
		TrRgn[0].VectorIndexList.push_back( vector<int>(Vect, Vect + sizeof(Vect)/sizeof(int)) );

		SendTransferRegions.push_back(TrRgn[0]);	//*/

		TrRgn[1].Region = GridRegion(0, 1, 0, 10, 0, 10);
		TrRgn[1].NeighList.push_back(2);

		int Vect2[] = {2};
		TrRgn[1].VectorIndexList.push_back( vector<int>(Vect2, Vect2 + sizeof(Vect2)/sizeof(int)) );

		SendTransferRegions.push_back(TrRgn[1]);	//*/

		TrRgn[2].Region = GridRegion(0, 10, 9, 10, 0, 10);
		TrRgn[2].NeighList.push_back(3);

		int Vect3[] = {3};
		TrRgn[2].VectorIndexList.push_back( vector<int>(Vect3, Vect3 + sizeof(Vect3)/sizeof(int)) );

		SendTransferRegions.push_back(TrRgn[2]);	//*/

		TrRgn[3].Region = GridRegion(0, 1, 0, 10, 0, 1);
		TrRgn[3].NeighList.push_back(4);

		int Vect4[] = {4};
		TrRgn[3].VectorIndexList.push_back( vector<int>(Vect4, Vect4 + sizeof(Vect4)/sizeof(int)) );

		SendTransferRegions.push_back(TrRgn[3]);	//*/

		TrRgn[4].Region = GridRegion(0, 10, 9, 10, 0, 1);
		TrRgn[4].NeighList.push_back(5);

		int Vect5[] = {5};
		TrRgn[4].VectorIndexList.push_back( vector<int>(Vect5, Vect5 + sizeof(Vect5)/sizeof(int)) );

		SendTransferRegions.push_back(TrRgn[4]);	//*/

		TrRgn[5].Region = GridRegion(0, 1, 9, 10, 0, 10);
		TrRgn[5].NeighList.push_back(6);

		int Vect6[] = {6};
		TrRgn[5].VectorIndexList.push_back( vector<int>(Vect6, Vect6 + sizeof(Vect6)/sizeof(int)) );

		SendTransferRegions.push_back(TrRgn[5]);	//*/

		TrRgn[6].Region = GridRegion(0, 1, 9, 10, 0, 1);
		TrRgn[6].NeighList.push_back(7);

		int Vect7[] = {7};
		TrRgn[6].VectorIndexList.push_back( vector<int>(Vect7, Vect7 + sizeof(Vect7)/sizeof(int)) );

		SendTransferRegions.push_back(TrRgn[6]);	//*/

		//Initial list of transfer regions

		printf("Initial list of regions\n");

		for(int c=0; c<SendTransferRegions.size(); c++){

			SendTransferRegions[c].Print();

		}

		//Find areas overlapping and create new transfer region

		for(int c=0; c<SendTransferRegions.size(); c++){

			//Find any overlapping areas recursively

			for(int c2=c+1; c2<SendTransferRegions.size(); c2++){

				GridRegion Region1 = SendTransferRegions[c].Region;
				GridRegion Region2 = SendTransferRegions[c2].Region;

				if(!Region1.TestOverlap(&Region2))
					continue;

				GridRegion Overlap = Region1.Overlap(&Region2);

				//Check to see if overlap lies entirely in excluded region

				bool Excl = false;

				for(int n=0; n<SendTransferRegions[c].ExcludeRegions.size(); n++){

					if(Overlap.FitsWithin( &(SendTransferRegions[c].ExcludeRegions[n]) )){
						Excl = true;
						break;
					}
				}

				if(Excl)
					continue;

				//Exclude overlap area and add new transfer region for overlap region

				SendTransferRegions[c].ExcludeRegions.push_back(Overlap);
				SendTransferRegions[c2].ExcludeRegions.push_back(Overlap);

				TransferRegion TrRgn;

				TrRgn.Region = Overlap;		//Overlap region

				TrRgn.AddNeighbours( &(SendTransferRegions[c]) );		//Add neighbours and vector lists to new transfer region
				TrRgn.AddNeighbours( &(SendTransferRegions[c2]) );	

				SendTransferRegions.push_back(TrRgn);		//Add overlap region to list
			}

		}

		//Remove any transfer regions with 0 cells

		for(int c=0; c<SendTransferRegions.size(); c++){

			if(SendTransferRegions[c].CountCells() == 0){

				SendTransferRegions.erase(SendTransferRegions.begin() + c);

				c--;
			}

		}

		//Final list of regions

		printf("\nFinal list of regions\n");

		for(int c=0; c<SendTransferRegions.size(); c++){

			SendTransferRegions[c].Print();

		}

		return true;
	}

	bool FindTransferRegions(int VectorBasis[][3], int nVectors, bool LoopBoundary[3]){
		
		//Find neighbour transfer areas

		bool Ret = FindNeighbours(VectorBasis, nVectors, LoopBoundary);

		if(!Ret)
			return false;


		//Initial list of transfer regions

		for(int i=0; i<nDiv; i++){

			for(int c=0; c<Partitions[i].NeighbourListSend.size(); c++){	

				TransferRegion TrRgn;	//Send

				TrRgn.Region = Partitions[i].NeighbourListSend[c].RegionLocalPartition;
				TrRgn.NeighList.push_back( Partitions[i].NeighbourListSend[c].NeighbourIndex );
				TrRgn.NeighIndexList.push_back( Partitions[i].NeighbourListSend[c].NeighbourListIndex );
				TrRgn.VectorIndexList.push_back( vector<int>(Partitions[i].NeighbourListSend[c].VectorIndexList) );

				Partitions[i].SendTransferRegions.push_back(TrRgn);
			}

			for(int c=0; c<Partitions[i].NeighbourListRecv.size(); c++){	

				TransferRegion TrRgn;	//Recv

				TrRgn.Region = Partitions[i].NeighbourListRecv[c].RegionLocalPartition;

				TrRgn.NeighList.push_back( Partitions[i].NeighbourListRecv[c].NeighbourIndex );
				TrRgn.NeighIndexList.push_back( Partitions[i].NeighbourListRecv[c].NeighbourListIndex );
				TrRgn.VectorIndexList.push_back( vector<int>(Partitions[i].NeighbourListRecv[c].VectorIndexList) );

				Partitions[i].RecvTransferRegions.push_back(TrRgn);
			}

		}


		//Find transfer regions

		for(int i=0; i<nDiv; i++){

			//Find areas overlapping and create new transfer region

			for(int c=0; c<Partitions[i].SendTransferRegions.size(); c++){

				//Find any overlapping areas recursively

				for(int c2=c+1; c2<Partitions[i].SendTransferRegions.size(); c2++){

					GridRegion Region1 = Partitions[i].SendTransferRegions[c].Region;
					GridRegion Region2 = Partitions[i].SendTransferRegions[c2].Region;

					if(!Region1.TestOverlap(&Region2))
						continue;

					GridRegion Overlap = Region1.Overlap(&Region2);

					//Check to see if overlap lies entirely in excluded region

					bool Excl = false;

					for(int n=0; n<Partitions[i].SendTransferRegions[c].ExcludeRegions.size(); n++){
						if(Overlap.FitsWithin( &(Partitions[i].SendTransferRegions[c].ExcludeRegions[n]) )){
							Excl = true;
							break;
						}
					}

					if(Excl)
						continue;

					//Exclude overlap area and add new transfer region for overlap region

					AddExcludeRegion(i, c, c2, Overlap);

				}

			}

		}

		//Remove any transfer regions with 0 cells

		vector<int>* SendIndexDec = new vector<int>[nDiv];
		vector<int>* RecvIndexDec = new vector<int>[nDiv];

		for(int i=0; i<nDiv; i++){

			int nSendRegions = (int)Partitions[i].SendTransferRegions.size();
			int nRecvRegions = (int)Partitions[i].RecvTransferRegions.size();

			SendIndexDec[i].resize(nSendRegions, 0);
			RecvIndexDec[i].resize(nRecvRegions, 0);

			int c0 = 0;
			for(int c=0; c<Partitions[i].SendTransferRegions.size(); c++){		//Send regions
				if(Partitions[i].SendTransferRegions[c].CountCells() == 0){
					Partitions[i].SendTransferRegions.erase(Partitions[i].SendTransferRegions.begin() + c);	//Remove entry
					for(int n=c0+1; n<nSendRegions; n++) SendIndexDec[i][n]++;								//Decrement indices references above this
					c--;
				}
				c0++;
			}

			c0 = 0;
			for(int c=0; c<Partitions[i].RecvTransferRegions.size(); c++){		//Recv regions
				if(Partitions[i].RecvTransferRegions[c].CountCells() == 0){
					Partitions[i].RecvTransferRegions.erase(Partitions[i].RecvTransferRegions.begin() + c);	//Remove entry
					for(int n=c0+1; n<nRecvRegions; n++) RecvIndexDec[i][n]++;									//Decrement indices references above this
					c--;
				}
				c0++;
			}

		}

		//Adjust reference indices

		for(int i=0; i<nDiv; i++){

			for(int c=0; c<Partitions[i].SendTransferRegions.size(); c++){		//Send regions
				for(int n=0; n<Partitions[i].SendTransferRegions[c].NeighList.size(); n++){
				
					int Dec = RecvIndexDec[Partitions[i].SendTransferRegions[c].NeighList[n]][Partitions[i].SendTransferRegions[c].NeighIndexList[n]];
				
					Partitions[i].SendTransferRegions[c].NeighIndexList[n] -= Dec;
				}
					
			}

			for(int c=0; c<Partitions[i].RecvTransferRegions.size(); c++){		//Recv regions

				int Dec = SendIndexDec[Partitions[i].RecvTransferRegions[c].NeighList[0]][Partitions[i].RecvTransferRegions[c].NeighIndexList[0]];
				
				Partitions[i].RecvTransferRegions[c].NeighIndexList[0] -= Dec;
			}

		}

		delete[] SendIndexDec;
		delete[] RecvIndexDec;

		return true;
	}

	void PrintNeighbours(){

		for(int n=0; n<nDiv; n++){

			printf("Neighbour %i has %i ghost regions:\n", n, Partitions[n].NeighbourListRecv.size());

			for(int i=0; i<Partitions[n].NeighbourListRecv.size(); i++){

				Decomposition<T>::Neighbour Neigh = Partitions[n].NeighbourListRecv[i];

				GridRegion Rgn = Neigh.RegionGlobal;

				printf("\tRecv %i from neighbour %i [%i] in region (%s) (%i nodes); by vectors: ", i, Neigh.NeighbourIndex, Neigh.NeighbourListIndex, VectorDisp(Rgn).Str, Neigh.NodeCount);

				for(int iv=0; iv<Neigh.VectorIndexList.size(); iv++)
					printf("%i%c ", Neigh.VectorIndexList[iv], (iv==Neigh.VectorIndexList.size()-1) ? ' ': ',');

				printf("\n");

			}

			printf("Neighbour %i has %i send regions:\n", n, Partitions[n].NeighbourListSend.size());

			for(int i=0; i<Partitions[n].NeighbourListSend.size(); i++){

				Decomposition<T>::Neighbour Neigh = Partitions[n].NeighbourListSend[i];

				GridRegion Rgn = Neigh.RegionGlobal;

				printf("\tSend %i to neighbour %i [%i] region (%s) (%i nodes); by vectors: ", i, Neigh.NeighbourIndex, Neigh.NeighbourListIndex, VectorDisp(Rgn).Str, Neigh.NodeCount);

				for(int iv=0; iv<Neigh.VectorIndexList.size(); iv++)
					printf("%i%c ", Neigh.VectorIndexList[iv], (iv==Neigh.VectorIndexList.size()-1) ? ' ': ',');

				printf("\n");

			}

			printf("\n");

		}


	}

	void PrintTransferRegions(){

		for(int n=0; n<nDiv; n++){

		/*	printf("Neighbour %i has %i ghost regions:\n", n, Partitions[n].NeighbourListRecv.size());

			for(int i=0; i<Partitions[n].NeighbourListRecv.size(); i++){

				Decomposition<int>::Neighbour Neigh = Partitions[n].NeighbourListRecv[i];

				GridRegion Rgn = Neigh.RegionGlobal;

				printf("\tRecv %i from neighbour %i in region (%i -> %i, %i -> %i, %i -> %i) (%i nodes)\n\t\tBy vectors: ", i, Neigh.NeighbourIndex, Rgn.x0, Rgn.x1, Rgn.y0, Rgn.y1, Rgn.z0, Rgn.z1, Neigh.NodeCount);

				for(int iv=0; iv<Neigh.VectorIndexList.size(); iv++)
					printf("%i%c ", Neigh.VectorIndexList[iv], (iv==Neigh.VectorIndexList.size()-1) ? ' ': ',');

				printf("\n");

			}	*/

			printf("Partition %i has %i send transfer regions:\n", n, Partitions[n].SendTransferRegions.size());

			for(int i=0; i<Partitions[n].SendTransferRegions.size(); i++){

				printf("[%i] ", i);

				Partitions[n].SendTransferRegions[i].Print();

			}

		//	printf("\n");

			printf("Partition %i has %i recv transfer regions:\n", n, Partitions[n].RecvTransferRegions.size());

			for(int i=0; i<Partitions[n].RecvTransferRegions.size(); i++){

				printf("[%i] ", i);

				Partitions[n].RecvTransferRegions[i].Print();

			}

			printf("\n");

		}

	}

	bool CheckTransferRegions(){

		bool Ret = true;
		
		//Check each send and recv region is paired with correct neighbour

		for(int i=0; i<nDiv; i++){
		
			for(int c=0; c<Partitions[i].SendTransferRegions.size(); c++){

				TransferRegion& TrRgn = (Partitions[i].SendTransferRegions[c]);

				for(int n=0; n<TrRgn.NeighList.size(); n++){

					TransferRegion& TrRgnPair = Partitions[TrRgn.NeighList[n]].RecvTransferRegions[TrRgn.NeighIndexList[n]];

					if(TrRgn.CountCells() != TrRgnPair.CountCells()){
						printf("Transfer region mismatch (send region partition %i [%i] with recv region partition %i [%i]\n", i, c, TrRgn.NeighList[n], TrRgn.NeighIndexList[n]);
						Ret = false;
					}

				}


			}
		
		}

		return Ret;
	}

	bool WriteDecompositionLog(const char* FileName){
	
		FILE* Log = fopen(FileName, "wb");

		if(!Log)
			return false;

		//Write out domain ranges

		int clength[7] = {0,0,0,0,0,0,0};

		long long nNodesTot = 0;

		for(int i=0; i<nDiv; i++){
		
			int tlength[7] = {1,1,1,1,1,1,1};

			GridRegion Rgn = Partitions[i].Region;

			for(int c=10; Rgn.x0 > c; c*=10) tlength[0]++;
			for(int c=10; Rgn.x1 > c; c*=10) tlength[1]++;
			for(int c=10; Rgn.y0 > c; c*=10) tlength[2]++;
			for(int c=10; Rgn.y1 > c; c*=10) tlength[3]++;
			for(int c=10; Rgn.z0 > c; c*=10) tlength[4]++;
			for(int c=10; Rgn.z1 > c; c*=10) tlength[5]++;

			for(int c=10; Partitions[i].NodeCount > c; c*=10) tlength[6]++;

			for(int n=0; n<7; n++)
				clength[n] = max( clength[n], tlength[n] );
		
			nNodesTot += Partitions[i].NodeCount;
		}

		fprintf(Log, "Domain Ranges:\n");

		int il = 1;
		for(int c=10; (nDiv-1) > c; c*=10) il++; 

		char HChar[256];
		memset(HChar, ' ', sizeof(HChar));

		HChar[4 + il + clength[0]] = 'x';
		HChar[10 + il + clength[0] + clength[1] + clength[2]] = 'y';
		HChar[16 + il + clength[0] + clength[1] + clength[2] + clength[3] + clength[4]] = 'z';

		strcpy(&HChar[21 + il + clength[0] + clength[1] + clength[2] + clength[3] + clength[4] + clength[5]], "Count");

		HChar[21 + il + clength[0] + clength[1] + clength[2] + clength[3] + clength[4] + clength[5] + 5] = '\0';

		fprintf(Log, "%s\n", HChar);

		for(int i=0; i<nDiv; i++){

			GridRegion Rgn = Partitions[i].Region;

			char CoordChar[256];
			memset(CoordChar, ' ', sizeof(CoordChar));

			int Ind  = 0;

			_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "(%d)         ", i);
			Ind += il + 3;

			_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "%lld           ", Rgn.x0);
			Ind += clength[0];
			_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, " -> %lld       ", Rgn.x1);
			Ind += clength[1] + 6;
			_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "%lld           ", Rgn.y0);
			Ind += clength[2];
			_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, " -> %lld       ", Rgn.y1);
			Ind += clength[3] + 6;
			_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "%lld           ", Rgn.z0);
			Ind += clength[4];
			_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, " -> %lld       ", Rgn.z1);
			Ind += clength[5] + 6;

			_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "%lli (%.2f%%)\0         ", Partitions[i].NodeCount, 100.0*Partitions[i].NodeCount/(double)nNodesTot);
			Ind += clength[6] + 6;


			fprintf(Log, "%s\n", CoordChar);

		}

		fclose(Log);
	
		return true;
	}

#ifdef MPI_INCLUDED

	//Distribute grid

	bool MPIDistributePartitions(Grid<T>* GridOut, int Offset[3]){
		
		int nBoundary[6] = {0, 0, 0, 0, 0, 0};

		return MPIDistributePartitions(GridOut, Offset, nBoundary);
	}

	bool MPIDistributePartitions(Grid<T>* GridOut, int Offset[3], int nBoundary[6]){

		int OffsetIn[3] = {0, 0, 0};

		if(MPIGridMode){

			_LatticeMPI->ObtainSubGrid(*GridOut, Partitions[MPIThreadId].Region, GridPosition(Offset[0], Offset[1], Offset[2]), nBoundary);

			return true;

		}else{

			return MPIScatterGrid(_Lattice, arrSz, OffsetIn, MPIGeoThread, GridOut, Offset, nBoundary, 1);

		}
	}

	template <class GridType>
	bool MPIScatterGrid(Grid<GridType>* GridIn, int SzGridIn[3], int OffsetIn[3], int ThreadIdMain, Grid<GridType>* GridPart, int Offset[3], int nBoundary[6], int nVector){

		long long Sz[3];
		Sz[0] = SzGridIn[0];
		Sz[1] = SzGridIn[1];
		Sz[2] = SzGridIn[2];

		return MPIScatterGrid(GridIn, Sz, OffsetIn, ThreadIdMain, GridPart, Offset, nBoundary, nVector);
	}

	template <class GridType>
	bool MPIScatterGrid(Grid<GridType>* GridIn, long long SzGridIn[3], int OffsetIn[3], int ThreadIdMain, Grid<GridType>* GridPart, int Offset[3], int nBoundary[6], int nVector){

		using namespace DecompositionDefinitions;

		MPI_Status Stat;

		//Boundary info

		int nThreads;
		MPI_Comm_size(MPIComm, &nThreads);
		
		int* nBoundaryList = (MPIThreadId == MPIGeoThread) ? new int[6*nThreads] : 0;

		MPI_Gather(nBoundary, 6, MPI_INT, nBoundaryList, 6, MPI_INT, MPIGeoThread, MPIComm);

		//Data transfer

		if(MPIThreadId == ThreadIdMain){

			//Find max partition size for array

			long long CountMax = 0;

			for(int i=0; i<nDiv; i++){

				GridRegion Rgn = Partitions[i].Region;

				long long Sz[3];

				Sz[0] = Rgn.SizeX() + nBoundaryList[i*6 + 0] + nBoundaryList[i*6 + 1];	//SizeX + Boundary(+X) + Boundary(-X)
				Sz[1] = Rgn.SizeY() + nBoundaryList[i*6 + 2] + nBoundaryList[i*6 + 3];	//SizeY + Boundary(+Y) + Boundary(-Y)
				Sz[2] = Rgn.SizeZ() + nBoundaryList[i*6 + 4] + nBoundaryList[i*6 + 5];	//SizeZ + Boundary(+Z) + Boundary(-Z)

				CountMax = max( Sz[0]*Sz[1]*Sz[2] , CountMax);
			}

			GridType* DataArr = new GridType[CountMax*nVector];

			//Send to each thread

			for(int i=0; i<nDiv; i++){

				GridRegion Rgn = Partitions[i].Region;

				long long Sz[3];

				Sz[0] = Rgn.SizeX() + nBoundaryList[i*6 + 0] + nBoundaryList[i*6 + 1];	//SizeX + Boundary(+X) + Boundary(-X)
				Sz[1] = Rgn.SizeY() + nBoundaryList[i*6 + 2] + nBoundaryList[i*6 + 3];	//SizeY + Boundary(+Y) + Boundary(-Y)
				Sz[2] = Rgn.SizeZ() + nBoundaryList[i*6 + 4] + nBoundaryList[i*6 + 5];	//SizeZ + Boundary(+Z) + Boundary(-Z)

				Grid<GridType> GrdSend(DataArr, true, false, Sz[0], Sz[1], Sz[2], nVector);

				GrdSend.ReadSubGrid(*GridIn, Rgn, GridPosition(nBoundaryList[i*6]+OffsetIn[0], nBoundaryList[i*6+2]+OffsetIn[1], nBoundaryList[i*6+4]+OffsetIn[2]));		//Read partition into array with offset and boundary

				//Add boundary layers

				for(int iv=0; iv<26; iv++){

					int v[3];
					memcpy(v, Q26Basis[iv], sizeof(int)*3);

					int* Bnd = &(nBoundaryList[i*6]);

					if((v[0] > 0 && Bnd[0] == 0) || (v[0] < 0 && Bnd[1] == 0)) continue;	//Determine if region required
					if((v[1] > 0 && Bnd[2] == 0) || (v[1] < 0 && Bnd[3] == 0)) continue;
					if((v[2] > 0 && Bnd[4] == 0) || (v[2] < 0 && Bnd[5] == 0)) continue;

					GridRegion R0 = GridRegion(0LL, Rgn.SizeX(), 0LL, Rgn.SizeY(), 0LL, Rgn.SizeZ());		//Relative to partition

					//Boundary region

					if(v[0] > 0) R0.SetRangeX(R0.x1, R0.x1 + Bnd[0]*v[0]);
					if(v[1] > 0) R0.SetRangeY(R0.y1, R0.y1 + Bnd[2]*v[1]);
					if(v[2] > 0) R0.SetRangeZ(R0.z1, R0.z1 + Bnd[4]*v[2]);

					if(v[0] < 0) R0.SetRangeX(R0.x0 + Bnd[1]*v[0], R0.x0);
					if(v[1] < 0) R0.SetRangeY(R0.y0 + Bnd[3]*v[1], R0.y0);
					if(v[2] < 0) R0.SetRangeZ(R0.z0 + Bnd[5]*v[2], R0.z0);

					//Loop boundaries

					GridRegion RgnGrid = R0.Translate(Rgn.x0+OffsetIn[0], Rgn.y0+OffsetIn[1], Rgn.z0+OffsetIn[2]);	//Region within grid

					long long Translation[3] = {0, 0, 0};

					if(RgnGrid.x0 < 0) Translation[0] += SzGridIn[0];		//Loop any coordinates over boundaries
					if(RgnGrid.y0 < 0) Translation[1] += SzGridIn[1];
					if(RgnGrid.z0 < 0) Translation[2] += SzGridIn[2];

					if(RgnGrid.x1 > SzGridIn[0]) Translation[0] -= SzGridIn[0];
					if(RgnGrid.y1 > SzGridIn[1]) Translation[1] -= SzGridIn[1];
					if(RgnGrid.z1 > SzGridIn[2]) Translation[2] -= SzGridIn[2];

					RgnGrid = RgnGrid.Translate(Translation[0], Translation[1], Translation[2]);

					//Read area into send grid

					GrdSend.ReadSubGrid(*GridIn, RgnGrid, GridPosition(nBoundaryList[i*6] + R0.x0, nBoundaryList[i*6+2] + R0.y0, nBoundaryList[i*6+4] + R0.z0));

				}

				//Send data

				if(i==MPIGeoThread){

					GridPart->ReadGrid(GrdSend, GridPosition(Offset[0], Offset[1], Offset[2]));
					
				}else{

					size_t SendCount = Sz[0]*Sz[1]*Sz[2]*nVector*sizeof(GridType);

					for(size_t SendInd=0; SendInd < SendCount; ){		//Send in chunks in case SendCount > INT_MAX

						size_t nSend = min((size_t)INT_MAX, SendCount-SendInd);

						MPI_Send(&(GrdSend.Data[SendInd/sizeof(GridType)]), (int)nSend, MPI_BYTE, i, 0, MPIComm);

						SendInd += nSend;
					}

				}

			}

			delete[] nBoundaryList;
			delete[] DataArr;

		}else{

			GridRegion Rgn = Partitions[MPIThreadId].Region;

			long long Sz[3];

			Sz[0] = Rgn.SizeX() + nBoundary[0] + nBoundary[1];
			Sz[1] = Rgn.SizeY() + nBoundary[2] + nBoundary[3];
			Sz[2] = Rgn.SizeZ() + nBoundary[4] + nBoundary[5];
		
			Grid<GridType> GrdRecv(Sz[0], Sz[1], Sz[2], (long long)nVector);

			size_t RecvCount = Sz[0]*Sz[1]*Sz[2]*nVector*sizeof(GridType);

			for(size_t RecvInd=0; RecvInd < RecvCount; ){		//Recv in chunks in case SendCount > INT_MAX

				size_t nRecv = min((size_t)INT_MAX, RecvCount-RecvInd);

				MPI_Recv(&(GrdRecv.Data[RecvInd/sizeof(GridType)]), (int)nRecv, MPI_BYTE, ThreadIdMain, 0, MPIComm, &Stat);

				RecvInd += nRecv;
			}

			//Read into grid with offset

			GridPart->ReadGrid(GrdRecv, GridPosition(Offset[0], Offset[1], Offset[2]));
		
		}
	
		return true;
	}

	//Gather grid

	template <class GridType>
	bool MPIGatherGrid(Grid<GridType>* GridOut, int ThreadIdMain, Grid<GridType>* GridPart, int nVector){

		int Offset[3] = {0, 0, 0};

		return MPIGatherGrid(GridOut, Offset, ThreadIdMain, GridPart, Offset, nVector, 0);
	}

	template <class GridType>
	bool MPIGatherGrid(Grid<GridType>* GridOut, int ThreadIdMain, Grid<GridType>* GridPart, int nVector, int MPITag){

		int Offset[3] = {0, 0, 0};

		return MPIGatherGrid<GridType>(GridOut, Offset, ThreadIdMain, GridPart, Offset, nVector, MPITag);
	}

	template <class GridType>
	bool MPIGatherGrid(Grid<GridType>* GridOut, int OffsetOut[3], int ThreadIdMain, Grid<GridType>* GridPart, int Offset[3], int nVector){

		return MPIGatherGrid<GridType>(GridOut, OffsetOut, ThreadIdMain, GridPart, Offset, nVector, 0);
	}

	template <class GridType>
	bool MPIGatherGrid(Grid<GridType>* GridOut, int OffsetOut[3], int ThreadIdMain, Grid<GridType>* GridPart, int Offset[3], int nVector, int MPITag){
	
		MPI_Status Stat;

		int nThreads;
		MPI_Comm_size(MPIComm, &nThreads);

	
		//Data transfer

		if(MPIThreadId == ThreadIdMain){

			//Find max partition size for array

			size_t CountMax = 0;

			for(int i=0; i<nDiv; i++){
				CountMax = max( (size_t)Partitions[i].Region.CountCells(), CountMax);
			}

			GridType* DataArr = new GridType[CountMax*nVector];

			for(int i=0; i<nDiv; i++){

				GridRegion Rgn = Partitions[i].Region;

				Grid<GridType> GrdRecv(DataArr, true, false, Rgn.SizeX(), Rgn.SizeY(), Rgn.SizeZ(), nVector);

				if(i == ThreadIdMain){

					GridRegion SubRegion((long long)Offset[0], Offset[0]+Rgn.SizeX(), (long long)Offset[1], Offset[1]+Rgn.SizeY(), (long long)Offset[2], Offset[2]+Rgn.SizeZ());

					GridOut->ReadSubGrid(*GridPart, SubRegion, GridPosition(Rgn.x0+OffsetOut[0], Rgn.y0+OffsetOut[1], Rgn.z0+OffsetOut[2]));

				}else{

					size_t RecvCount = Rgn.CountCells()*sizeof(GridType)*nVector;

					for(size_t RecvInd=0; RecvInd < RecvCount; ){		//Recv in chunks in case SendCount > INT_MAX

						size_t nRecv = min((size_t)INT_MAX, RecvCount-RecvInd);

						MPI_Recv(&(GrdRecv.Data[RecvInd/sizeof(GridType)]), (int)nRecv, MPI_BYTE, i, MPITag, MPIComm, &Stat);

						RecvInd += nRecv;
					}

					GridOut->ReadGrid(GrdRecv, GridPosition(Rgn.x0+OffsetOut[0], Rgn.y0+OffsetOut[1], Rgn.z0+OffsetOut[2]));

				}

			}

			delete[] DataArr;

		}else{

			GridRegion Rgn = Partitions[MPIThreadId].Region;

			Grid<GridType> GrdSend(Rgn.SizeX(), Rgn.SizeY(), Rgn.SizeZ(), (long long)nVector);

			GrdSend.ReadSubGrid(*GridPart, GridRegion((long long)Offset[0], Offset[0]+Rgn.SizeX(), (long long)Offset[1], Offset[1]+Rgn.SizeY(), (long long)Offset[2], Offset[2]+Rgn.SizeZ()));

			size_t SendCount = Rgn.CountCells()*sizeof(GridType)*nVector;

			for(size_t SendInd=0; SendInd < SendCount; ){		//Send in chunks in case SendCount > INT_MAX

				size_t nSend = min((size_t)INT_MAX, SendCount-SendInd);

				MPI_Send(&(GrdSend.Data[SendInd/sizeof(GridType)]), (int)nSend, MPI_BYTE, ThreadIdMain, MPITag, MPIComm);

				SendInd += nSend;
			}

		}

		return true;
	}

#endif

private:

	bool InterpretDirStr(const char* Str, Direction* DirList, int* nDecomp){

		int sLen = (int)strlen(Str);

		*nDecomp = 0;

		for(int c=0; c<sLen; c++){

			if(Str[c] == ' ' || Str[c] == ',' || Str[c] == ';' || Str[c] == '\t' || Str[c] == '-')
				continue;

			if(Str[c] == 'x' || Str[c] == 'X'){

				DirList[(*nDecomp)++] = Dir_X;
				
			}else if(Str[c] == 'y' || Str[c] == 'Y'){

				DirList[(*nDecomp)++] = Dir_Y;

			}else if(Str[c] == 'z' || Str[c] == 'Z'){

				DirList[(*nDecomp)++] = Dir_Z;

			}

			if((*nDecomp) == 3)
				break;
		}

		return (*nDecomp) > 0;
	}

	void AddExcludeRegion(int PartitionInd, int SendRegion1Ind, int SendRegion2Ind, GridRegion Overlap){

		//Add exclude region and create new send region
	
		Partitions[PartitionInd].SendTransferRegions[SendRegion1Ind].ExcludeRegions.push_back(Overlap);
		Partitions[PartitionInd].SendTransferRegions[SendRegion2Ind].ExcludeRegions.push_back(Overlap);

		TransferRegion TrRgnSendNew;		//New send region

		TrRgnSendNew.Region = Overlap;		//Overlap region

		int SendIndexNew = (int)Partitions[PartitionInd].SendTransferRegions.size();			//Index of new send region

		//Add new recv regions to neighbour

		vector<int> NeighRecvNewInd;
		vector<vector<int> > NeighRecvNewVectors;

		for(int c=0; c<2; c++){		//SendRegion1Ind and SendRegion2Ind

			TransferRegion& TrRgn = Partitions[PartitionInd].SendTransferRegions[((c==0) ? SendRegion1Ind : SendRegion2Ind)];

			for(int n=0; n<TrRgn.NeighList.size(); n++){		//List of send regions

				int NeighInd = TrRgn.NeighList[n];				//Neighbour index
				int NeighRgnInd = TrRgn.NeighIndexList[n];		//Recv transfer region index

				TransferRegion* RecvTrRgn = &(Partitions[NeighInd].RecvTransferRegions[NeighRgnInd]);	//Corresponding recv region

				//Add exclude region

				long long Translation[3];

				Translation[0] = RecvTrRgn->Region.x0 - TrRgn.Region.x0;
				Translation[1] = RecvTrRgn->Region.y0 - TrRgn.Region.y0;
				Translation[2] = RecvTrRgn->Region.z0 - TrRgn.Region.z0;

				GridRegion RecvRgnOverlap = Overlap.Translate(Translation[0], Translation[1], Translation[2]);

				RecvTrRgn->ExcludeRegions.push_back(RecvRgnOverlap);

				//Ignore if duplicate region

				bool Dupl = false;

				for(int i=0; i<NeighRecvNewInd.size(); i++){
					if(NeighInd == NeighRecvNewInd[i] && RecvTrRgn->VectorIndexList[0] == NeighRecvNewVectors[i]){
						Dupl = true;
						break;
					}
				}

				if(Dupl)
					continue;

				//New recv region

				TransferRegion TrRgnNew;

				TrRgnNew.Region = RecvRgnOverlap;													//Overlap region
				TrRgnNew.NeighList.push_back( RecvTrRgn->NeighList[0] );							//Neighbour index
				TrRgnNew.NeighIndexList.push_back( SendIndexNew );									//Set to new send region index
				TrRgnNew.VectorIndexList.push_back( vector<int>(RecvTrRgn->VectorIndexList[0]) );	//Vectors

				int RecvIndNew = (int)Partitions[NeighInd].RecvTransferRegions.size();				//Index of new recv region

				//Add new recv region ref to send region
			
				TrRgnSendNew.NeighList.push_back( NeighInd );											//Neigh index
				TrRgnSendNew.NeighIndexList.push_back( RecvIndNew );									//Index of recv region in neighbour
				TrRgnSendNew.VectorIndexList.push_back( vector<int>(RecvTrRgn->VectorIndexList[0]) );	//Vectors

				//Add new recv region

				Partitions[NeighInd].RecvTransferRegions.push_back(TrRgnNew);	//New recv region (must be after accessing 'RecvTrRgn->VectorIndexList[0]' to maintain ptr RecvTrRgn)
			
				NeighRecvNewInd.push_back( NeighInd );
				NeighRecvNewVectors.push_back( TrRgnNew.VectorIndexList[0] );
			}

		}

		Partitions[PartitionInd].SendTransferRegions.push_back(TrRgnSendNew);				//Add send overlap region to list
	
	}

#ifdef MPI_INCLUDED

	void CountNodesNeighMPI(){

		long long nRgnTot = 0;

		for(int i=0; i<nDiv; i++){

			nRgnTot += Partitions[i].NeighbourListRecv.size();

		}

		long long* Count = new long long[nRgnTot];

		memset(Count, 0, sizeof(long long)*nRgnTot);

		long long** CountList = new long long*[nDiv];

		nRgnTot = 0;

		for(int i=0; i<nDiv; i++){
		
			CountList[i] = &Count[nRgnTot];

			nRgnTot += Partitions[i].NeighbourListRecv.size();
		}

		//Count nodes

		if(!MPIGridMode){

			for(int i=0; i<nDiv; i++){

				//Main thread counts nodes

				if(MPIThreadId == MPIGeoThread){	

					for(int i2=0; i2<Partitions[i].NeighbourListRecv.size(); i2++){

						CountList[i][i2] = _Lattice->CountOccurences(_NodeValue, Partitions[i].NeighbourListRecv[i2].RegionGlobal);	//Number of fluid nodes in transfer area

					}

				}

				//Distribute

				MPI_Bcast(CountList[i], (int)Partitions[i].NeighbourListRecv.size(), MPI_LONG_LONG, MPIGeoThread, MPIComm);

			}

		}else{

			GridRegion GrdRgn = _LatticeMPI->Regions[MPIThreadId];

			for(int i=0; i<nDiv; i++){

				for(int i2=0; i2<Partitions[i].NeighbourListRecv.size(); i2++){

					GridRegion Overlap = GrdRgn.Overlap(Partitions[i].NeighbourListRecv[i2].RegionGlobal);

					if(Overlap.IsZero())
						continue;

					GridRegion RgnLocal = Overlap.Translate(-GrdRgn.x0, -GrdRgn.y0, -GrdRgn.z0);

					CountList[i][i2] = _LatticeMPI->GridLocal->CountOccurences(_NodeValue, RgnLocal);

				}

			}
		
			//Distribute count

			long long* CountGlobal = new long long[nRgnTot];

			MPI_Allreduce(Count, CountGlobal, (int)nRgnTot, MPI_LONG_LONG, MPI_SUM, MPIComm);

			memcpy(Count, CountGlobal, sizeof(long long)*nRgnTot);

			delete[] CountGlobal;
		
		}

		//Partition

		for(int i=0; i<nDiv; i++){

			for(int i2=0; i2<Partitions[i].NeighbourListRecv.size(); i2++){

				Partitions[i].NeighbourListRecv[i2].NodeCount = CountList[i][i2];

				int NeighInd = Partitions[i].NeighbourListRecv[i2].NeighbourIndex;
				int NeighSendInd = Partitions[i].NeighbourListRecv[i2].NeighbourListIndex;

				Partitions[NeighInd].NeighbourListSend[NeighSendInd].NodeCount = CountList[i][i2];

			//	if(MPIThreadId==0)
			//		printf("Partition [%i] neigh [%i] count = %lli\n", i, i2, CountList[i][i2]);

			}

		}

		//Finalise

		delete[] Count;
		delete[] CountList;
	}

#endif

	bool FindNeighbours(int VectorBasis[][3], int nVectors, bool LoopBoundary[3]){
	
		for(int i=0; i<nDiv; i++){

			GridRegion R0 = Partitions[i].Region;
		
			for(int iv=0; iv<nVectors; iv++){

				int v[3];
				memcpy(v, VectorBasis[iv], sizeof(int)*3);
				
				GridRegion Region = R0;

				//Overlap region
			
				if(v[0] < 0) Region.SetRangeX(R0.x0 + v[0], R0.x0);
				if(v[1] < 0) Region.SetRangeY(R0.y0 + v[1], R0.y0);
				if(v[2] < 0) Region.SetRangeZ(R0.z0 + v[2], R0.z0);

				if(v[0] > 0) Region.SetRangeX(R0.x1, R0.x1 + v[0]);
				if(v[1] > 0) Region.SetRangeY(R0.y1, R0.y1 + v[1]);
				if(v[2] > 0) Region.SetRangeZ(R0.z1, R0.z1 + v[2]);

				//Boundaries

				if((Region.x0 < 0 || Region.x1 > arrSz[0]) && !LoopBoundary[0]) continue;
				if((Region.y0 < 0 || Region.y1 > arrSz[1]) && !LoopBoundary[1]) continue;
				if((Region.z0 < 0 || Region.z1 > arrSz[2]) && !LoopBoundary[2]) continue;

				long long Translation[3] = {0, 0, 0};

				if(Region.x0 < 0) Translation[0] = arrSz[0];		//Loop any coordinates over boundaries
				if(Region.y0 < 0) Translation[1] = arrSz[1];
				if(Region.z0 < 0) Translation[2] = arrSz[2];

				if(Region.x1 > arrSz[0]) Translation[0] = -arrSz[0];
				if(Region.y1 > arrSz[1]) Translation[1] = -arrSz[1];
				if(Region.z1 > arrSz[2]) Translation[2] = -arrSz[2];

				Region = Region.Translate(Translation[0], Translation[1], Translation[2]);
			
		//		if(i==nDiv-1)
		//		printf("Overlap range vector (%i, %i, %i) = (%i -> %i, %i -> %i, %i -> %i)\n", v[0], v[1], v[2], Region.x0, Region.x1, Region.y0, Region.y1, Region.z0, Region.z1);
			
				//Find overlapping neigbours

				for(int i2=0; i2<nDiv; i2++){
				
					if(i2 == i)
						continue;
				
					if(!Region.TestOverlap(&(Partitions[i2].Region)))	//Test overlap
						continue;
		
					GridRegion Overlap = Region.Overlap(&(Partitions[i2].Region));	//Overlap region

		/*			if(i==nDiv-1 && iv==1)
						printf("Boundary [1] overlaps with neighbour [%i] (%i -> %i, %i -> %i, %i -> %i); region (%i -> %i, %i -> %i, %i -> %i)\n", i2,
						Partitions[i2].Region.x0, Partitions[i2].Region.x1, Partitions[i2].Region.y0, Partitions[i2].Region.y1, Partitions[i2].Region.z0, Partitions[i2].Region.z1,
						Overlap.x0, Overlap.x1, Overlap.y0, Overlap.y1, Overlap.z0, Overlap.z1);	//*/

					//Add neighbour

					Neighbour Neigh;

					GridRegion RegionLocal = Overlap.Translate(-R0.x0-Translation[0], -R0.y0-Translation[1], -R0.z0-Translation[2]);
					GridRegion RegionNeighLocal = Overlap.Translate(-Partitions[i2].Region.x0, -Partitions[i2].Region.y0, -Partitions[i2].Region.z0);

					Neigh.NeighbourIndex = i2;						//Index of neighbour
					Neigh.RegionGlobal = Overlap;					//Overlap region in lattice coordinates
					Neigh.RegionLocalPartition = RegionLocal;		//Overlap region in local partition coordinates
					Neigh.RegionLocalNeighbour = RegionNeighLocal;	//Overlap region in local neighbour coordinates

					if(!WeightMode){
					
						Neigh.NodeCount = Overlap.CountCells();

					}
#ifdef MPI_INCLUDED
					else if(!MPIMode){
				
						Neigh.NodeCount = _Lattice->CountOccurences(_NodeValue, Overlap);	//Number of fluid nodes in transfer area
					
					}
#endif

					for(int iv2=0; iv2<nVectors; iv2++){	//Find all vectors mapping to this region

						if((v[0] == 0 || VectorBasis[iv2][0] == v[0]) && (v[1] == 0 || VectorBasis[iv2][1] == v[1]) && (v[2] == 0 || VectorBasis[iv2][2] == v[2])){

							Neigh.VectorIndexList.push_back(iv2);	//Add vector to list
						}
					}

					Partitions[i].NeighbourListRecv.push_back(Neigh);		//Local partition recv

					int NeighListRecvInd = (int)Partitions[i].NeighbourListRecv.size() - 1;	//Index of recv neighbour in list

					//Send area within neighbour

					Neighbour NeighSend = Neigh;

					NeighSend.NeighbourIndex = i;
					NeighSend.RegionLocalPartition = RegionNeighLocal;		//Overlap region in neighbour partition coordinates
					NeighSend.RegionLocalNeighbour = RegionLocal;			//Overlap region in local partition coordinates

					Partitions[i2].NeighbourListSend.push_back(NeighSend);		//Neighbour area to send

					int NeighListSendInd = (int)Partitions[i2].NeighbourListSend.size() - 1;	//Index of send neighbour in list

					//Pair send and recv regions

					Partitions[i].NeighbourListRecv[NeighListRecvInd].NeighbourListIndex = NeighListSendInd;
					Partitions[i2].NeighbourListSend[NeighListSendInd].NeighbourListIndex = NeighListRecvInd;
				}
				
			}

		}

#ifdef MPI_INCLUDED

		if(WeightMode && MPIMode){
		
			CountNodesNeighMPI();
		
		}

#endif
	
		return true;
	}

#ifdef MPI_INCLUDED

	void ObtainNodeDistributionsMPI(GridRegion* Regions, int n, long long** Distributions, Direction Dir){

		if(!WeightMode){
		
			for(int i=0; i<n; i++){

				long long nX = Regions[i].SizeX();
				long long nY = Regions[i].SizeY();
				long long nZ = Regions[i].SizeZ();

				long long Sz = (Dir==Dir_X) ? nX : ( (Dir==Dir_Y) ? nY : nZ );

				long long Count = (Dir==Dir_X) ? nY*nZ : ( (Dir==Dir_Y) ? nX*nZ : nX*nY );

				for(long long c=0; c<Sz; c++)
					Distributions[i][c] = Count;

			}
		
			return;
		}

		long long* Dist = 0;
		long long SzLast = 0;

		for(int i=0; i<n; i++){

			long long Sz = (Dir==Dir_X) ? Regions[i].SizeX() : ( (Dir==Dir_Y) ? Regions[i].SizeY() : Regions[i].SizeZ() );

			if(Sz > SzLast){
				delete[] Dist;
				Dist = new long long[Sz];
				SzLast = Sz;
			}

			memset(Dist, 0, sizeof(long long)*Sz);

			GridRegion LocalRgn = _LatticeMPI->Regions[MPIThreadId].Overlap(Regions[i]);

			if(!LocalRgn.IsZero()){

				long long IndStart = (Dir==Dir_X) ? (LocalRgn.x0-Regions[i].x0) : ( (Dir==Dir_Y) ? (LocalRgn.y0-Regions[i].y0) : (LocalRgn.z0-Regions[i].z0) );

		//		printf("Thread %i map region (%lli -> %lli, %lli -> %lli, %lli -> %lli) from index %i\n", MPIThreadId, LocalRgn.x0, LocalRgn.x1, LocalRgn.y0, LocalRgn.y1, LocalRgn.z0, LocalRgn.z1, IndStart);
		//		fflush(stdout);

				LocalRgn = LocalRgn.Translate(-_LatticeMPI->Regions[MPIThreadId].x0, -_LatticeMPI->Regions[MPIThreadId].y0, -_LatticeMPI->Regions[MPIThreadId].z0);

				for(long long z=0; z<LocalRgn.SizeZ(); z++)
				for(long long y=0; y<LocalRgn.SizeY(); y++)
				for(long long x=0; x<LocalRgn.SizeX(); x++){
			
					if(IsActiveNode( _LatticeMPI->GridLocal->Get(LocalRgn.x0 + x, LocalRgn.y0 + y, LocalRgn.z0 + z) )){

						Dist[IndStart + ((Dir == Dir_X) ? x : ((Dir == Dir_Y) ? y : z))]++;

					}
			
				}

			}

			MPI_Allreduce(Dist, Distributions[i], (int)Sz, MPI_LONG_LONG, MPI_SUM, MPIComm);

		}

		if(Dist!=0)
			delete[] Dist;
	
	}

	bool DecomposeMPI(){

		int d[3];	//The three decomposition directions

		d[0]  = (Dirs[0]==Dir_X) ? 0 : ( (Dirs[0]==Dir_Y) ? 1 : 2 );	//Direction of decomposition 1
		d[1]  = (Dirs[1]==Dir_X) ? 0 : ( (Dirs[1]==Dir_Y) ? 1 : 2 );	//Direction of decomposition 2
		d[2]  = (Dirs[2]==Dir_X) ? 0 : ( (Dirs[2]==Dir_Y) ? 1 : 2 );	//Direction of decomposition 3

		double Ratio01 = (double)arrSz[d[0]] / (double)arrSz[d[1]];
		double Ratio02 = (double)arrSz[d[0]] / (double)arrSz[d[2]];

		double N = pow( (double)nDiv*Ratio01*Ratio02 , 1.0/3.0);

		int nDivA = min(nDiv, (int)round(N));	//Number of partitions in direction 1

		//printf("Ratio01 = %.2e\nRatio02 = %.2e\nNz = %.2e\n", Ratio01, Ratio02, N);
		//printf("Partition dir1 into %i\n", nDivA);

		//Partition in direction 1

		Partition* RegionsDir1 = new Partition[nDivA];

		int* nSubDivA = new int[nDivA];				//Number of sub divisions in each partition
		double* nNodesDivA = new double[nDivA];		//Number of fluid nodes in each partition

		long long nNodesTot = DomainPartition.NodeCount;

		for(int i=0; i<nDivA; i++){

			nSubDivA[i] = nDiv/nDivA + ( (i < nDiv%nDivA) ? 1 : 0 );

			nNodesDivA[i] = (double)nSubDivA[i] * ((double)nNodesTot/(double)nDiv);

			//printf("Partition [%i] into %i with %f nodes\n", i, nSubDivA[i], nNodesDivA[i]);
		}

		long long** Dist = new long long*[nDiv];

		Dist[0] = new long long[arrSz[d[0]]];

		ObtainNodeDistributionsMPI(&(DomainPartition.Region), 1, Dist, Dirs[0]);

	/*	if(MPIThreadId==0)
		for(int c=0; c<arrSz[d[0]]; c++){
			printf("%i, ", Dist[0][c]);

			if((c+1)%10 == 0)
				printf("\n");
		}	//*/

		bool r = DividePartition2(nDivA, Dirs[0], DomainPartition, Dist[0], nNodesDivA, RegionsDir1, 0);

		delete[] Dist[0];

		if(!r){
			delete[] nSubDivA;
			delete[] nNodesDivA;
			delete[] RegionsDir1;
			delete[] Dist;
			return false;
		}


		//Partition in direction 2

		GridRegion* RegionList = new GridRegion[nDivA];

		for(int c=0; c<nDivA; c++){
			RegionList[c] = RegionsDir1[c].Region;
			Dist[c] = new long long[arrSz[d[1]]];
		}

		ObtainNodeDistributionsMPI(RegionList, nDivA, Dist, Dirs[1]);

		delete[] RegionList;

		Partition** RegionsDir2 = new Partition*[nDivA];

		int* nDivB = new int[nDivA];						//
		int** nSubDivB = new int*[nDivA];					//Number of sub divisions in each partition

		int nSubDivBTot = 0;

		for(int i=0; i<nDivA; i++){

			//Number of partitions in direction 2

			double Ratio12 = (double)arrSz[d[1]] / (double)arrSz[d[2]];

			double N2 = sqrt( (double)nSubDivA[i]*Ratio12 );

			nDivB[i] = min(nSubDivA[i], (int)round(N2));

			nSubDivB[i] = new int[nDivB[i]];

			nSubDivBTot += nDivB[i];

			//	printf("Partition [%i] dir2 into %i\n", i, nDivB);

			RegionsDir2[i] = new Partition[nDivB[i]];

			double* nNodesDivB = new double[nDivB[i]];			//Number of fluid nodes in each partition

			for(int c=0; c<nDivB[i]; c++){

				nSubDivB[i][c] = nSubDivA[i]/nDivB[i] + ( (c < nSubDivA[i]%nDivB[i]) ? 1 : 0 );

				nNodesDivB[c] = (double)nSubDivB[i][c] * ((double)RegionsDir1[i].NodeCount/(double)nSubDivA[i]);

				//	printf("Partition [%i][%i] into %i with %f nodes\n", i, c, nSubDivB[c], nNodesDivB[c]);

			}

			r = DividePartition2(nDivB[i], Dirs[1], RegionsDir1[i], Dist[i], nNodesDivB, RegionsDir2[i], 0);

			delete[] nNodesDivB;

			if(!r){
				for(int c=0; c<=i; c++){
					delete[] RegionsDir2[c];
					delete[] nSubDivB[c];
				}
				delete[] RegionsDir2;
				delete[] nSubDivB;
				delete[] nSubDivA;
				delete[] nNodesDivA;
				delete[] RegionsDir1;
				for(int c=0; c<nDivA; c++)
					delete[] Dist[c];
				delete[] Dist;
				return false;
			}

		}

		for(int c=0; c<nDivA; c++)
			delete[] Dist[c];

		//Partition in direction 3

		GridRegion* RegionListB = new GridRegion[nSubDivBTot];

		int nDivCount = 0;

		for(int c=0;  c <nDivA;    c++ )
		for(int c2=0; c2<nDivB[c]; c2++){

			RegionListB[nDivCount] = RegionsDir2[c][c2].Region;
			Dist[nDivCount] = new long long[arrSz[d[2]]];

			nDivCount++;
		}

		ObtainNodeDistributionsMPI(RegionListB, nSubDivBTot, Dist, Dirs[2]);

		delete[] RegionListB;

		int OutIndex = 0;
		int nDivCount2 = 0;

		for(int i=0; i<nDivA; i++)
		for(int i2=0; i2<nDivB[i]; i2++){

			double* nNodesDivC = new double[nSubDivB[i][i2]];

			for(int c=0; c<nSubDivB[i][i2]; c++)
				nNodesDivC[c] = (double)RegionsDir2[i][i2].NodeCount / (double)nSubDivB[i][i2];

			r = DividePartition2(nSubDivB[i][i2], Dirs[2], RegionsDir2[i][i2], Dist[nDivCount2], nNodesDivC, Partitions, OutIndex);

			delete[] nNodesDivC;

			if(!r){
				for(int c=0; c<nDivA; c++){
					delete[] RegionsDir2[c];
					delete[] nSubDivB[c];
				}
				delete[] RegionsDir2;
				delete[] nSubDivB;
				delete[] nSubDivA;
				delete[] nNodesDivA;
				delete[] RegionsDir1;
				for(int c=0; c<nDivA; c++)
					delete[] Dist[c];
				delete[] Dist;

				return false;
			}

			OutIndex += nSubDivB[i][i2];
			nDivCount2++;
		}



		//Finalise

		for(int c=0; c<nDivA; c++){
			delete[] RegionsDir2[c];
			delete[] nSubDivB[c];
		}

		delete[] RegionsDir2;
		delete[] nSubDivB;
		delete[] nSubDivA;
		delete[] nNodesDivA;
		delete[] RegionsDir1;

		for(int c=0; c<nDivA; c++)
			delete[] Dist[c];

		delete[] Dist;

		return true;
	}

#endif

	bool Decompose(){

#ifdef MPI_INCLUDED

		if(!MPIMode){
		
			return DecomposeMain();
		
		}

		//Decompose on main thread and distribute result to all threads

		int tId;
		MPI_Comm_rank(MPIComm, &tId);	//Get thread ID in local group

		bool Ret;

		if(tId == MPIGeoThread){

			Ret = DecomposeMain();

		}

		MPI_Bcast(&Ret, sizeof(bool), MPI_CHAR, MPIGeoThread, MPIComm);

		if(!Ret)
			return false;

		MPI_Bcast(Partitions, sizeof(Partition)*nDiv, MPI_CHAR, MPIGeoThread, MPIComm);		//Distribute partitions list

		return true;

#else

		return DecomposeMain();

#endif

	}

	bool DecomposeMain(){

		if(nDirs == 1){				//1D decomp

			return Decompose1D();

		}else if(nDirs == 2){		//2D decomp

			return Decompose2D();

		}else{						//3D decomp

			return Decompose3D();
				
		}
	
	}

	bool Decompose1D(){

		long long* CountBuf = new long long[_NodeValueCount];

		double* nNodesPerDiv = new double[nDiv];

		long long nNodesTot = WeightMode ? _Lattice->CountOccurences(_NodeValue, _NodeValueCount, CountBuf) : arrSz[0]*arrSz[1]*arrSz[2];

		for(int i=0; i<nDiv; i++)
			nNodesPerDiv[i] = (double)nNodesTot / (double)nDiv;

		bool r = DividePartition(nDiv, Dirs[0], DomainPartition, nNodesPerDiv, Partitions, 0);

		delete[] nNodesPerDiv;

		delete[] CountBuf;
		return r;
	}

	bool Decompose2D(){

		long long* CountBuf = new long long[_NodeValueCount];

		int d[2];	//The two decomposition directions

		d[0]  = (Dirs[0]==Dir_X) ? 0 : ( (Dirs[0]==Dir_Y) ? 1 : 2 );	//Direction of decomposition 1
		d[1]  = (Dirs[1]==Dir_X) ? 0 : ( (Dirs[1]==Dir_Y) ? 1 : 2 );	//Direction of decomposition 2
		
		int d0 = 0;														//Remaining direction
		if(d[0] == d0 || d[1] == d0) d0++;
		if(d[0] == d0 || d[1] == d0) d0++;

		double Ratio = (double)arrSz[d[0]] / (double)arrSz[d[1]];
		double N = sqrt( (double)nDiv * Ratio );

//		printf("Ratio = %.2e\nNz = %.2e\n", Ratio, N);

		int nDivA = min(nDiv, (int)round( N ));		//Number of divisions in direction 1

//		printf("nDivA = %i\n", nDivA);

		//Partition in direction 1

	//	printf("Partition dir1 into %i\n", nDivA);

		Partition* RegionsDir1 = new Partition[nDivA];

		int* nDivB = new int[nDivA];				//Number of sub divisions in each partition
		double* nNodesDivA = new double[nDivA];		//Number of fluid nodes in each partition

		long long nNodesTot = WeightMode ? _Lattice->CountOccurences(_NodeValue, _NodeValueCount, CountBuf) : arrSz[0]*arrSz[1]*arrSz[2];

		for(int i=0; i<nDivA; i++){

			nDivB[i] = nDiv/nDivA + ( (i < nDiv%nDivA) ? 1 : 0 );

			nNodesDivA[i] = (double)nDivB[i] * ((double)nNodesTot/(double)nDiv);

		}

		bool r = DividePartition(nDivA, Dirs[0], DomainPartition, nNodesDivA, RegionsDir1, 0);

		if(!r){
			delete[] nNodesDivA;
			delete[] RegionsDir1;
			delete[] CountBuf;
			return false;
		}

		//Partition in direction 2

		int OutIndex = 0;

		for(int i=0; i<nDivA; i++){

			double* nNodesDivB = new double[nDivB[i]];

			long long nNodesDivA_Rounded = WeightMode ? _Lattice->CountOccurences(_NodeValue, _NodeValueCount, CountBuf, RegionsDir1[i].Region) : RegionsDir1[i].Region.CountCells();

			for(int c=0; c<nDivB[i]; c++)
				nNodesDivB[c] = (double)nNodesDivA_Rounded / (double)nDivB[i];

		//	printf("Partition [%i] into %i with %f nodes\n", i, nDivB[i], (double)nNodesDivA_Rounded / (double)nDivB[i]);

			r = DividePartition(nDivB[i], Dirs[1], RegionsDir1[i], nNodesDivB, Partitions, OutIndex);

			delete[] nNodesDivB;

			if(!r){
				delete[] nNodesDivA;
				delete[] nDivB;
				delete[] RegionsDir1;
				delete[] CountBuf;
				return false;
			}

			OutIndex += nDivB[i];

		}

		delete[] nNodesDivA;
		delete[] nDivB;
		delete[] RegionsDir1;
		delete[] CountBuf;

		return true;
	}

	bool Decompose3D(){

		long long* CountBuf = new long long[_NodeValueCount];

		int d[3];	//The three decomposition directions

		d[0]  = (Dirs[0]==Dir_X) ? 0 : ( (Dirs[0]==Dir_Y) ? 1 : 2 );	//Direction of decomposition 1
		d[1]  = (Dirs[1]==Dir_X) ? 0 : ( (Dirs[1]==Dir_Y) ? 1 : 2 );	//Direction of decomposition 2
		d[2]  = (Dirs[2]==Dir_X) ? 0 : ( (Dirs[2]==Dir_Y) ? 1 : 2 );	//Direction of decomposition 3

		double Ratio01 = (double)arrSz[d[0]] / (double)arrSz[d[1]];
		double Ratio02 = (double)arrSz[d[0]] / (double)arrSz[d[2]];

		double N = pow( (double)nDiv*Ratio01*Ratio02 , 1.0/3.0);

//		printf("Ratio01 = %.2e\nRatio02 = %.2e\nNz = %.2e\n", Ratio01, Ratio02, N);

		int nDivA = min(nDiv, (int)round(N));								//Number of partitions in direction 1

//		printf("Partition dir1 into %i\n", nDivA);

		//Partition in direction 1

		Partition* RegionsDir1 = new Partition[nDivA];

		int* nSubDivA = new int[nDivA];				//Number of sub divisions in each partition
		double* nNodesDivA = new double[nDivA];		//Number of fluid nodes in each partition

		long long nNodesTot = WeightMode ? _Lattice->CountOccurences(_NodeValue, _NodeValueCount, CountBuf) : arrSz[0]*arrSz[1]*arrSz[2];

		for(int i=0; i<nDivA; i++){

			nSubDivA[i] = nDiv/nDivA + ( (i < nDiv%nDivA) ? 1 : 0 );

			nNodesDivA[i] = (double)nSubDivA[i] * ((double)nNodesTot/(double)nDiv);

		//	printf("Partition [%i] into %i with %f nodes\n", i, nSubDivA[i], nNodesDivA[i]);
		}

		bool r = DividePartition(nDivA, Dirs[0], DomainPartition, nNodesDivA, RegionsDir1, 0);

		if(!r){
			delete[] nSubDivA;
			delete[] nNodesDivA;
			delete[] RegionsDir1;
			delete[] CountBuf;
			return false;
		}

		//Partition in direction 2

		int OutIndex = 0;

		for(int i=0; i<nDivA; i++){

			//Number of partitions in direction 2

			double Ratio12 = (double)arrSz[d[1]] / (double)arrSz[d[2]];

			double N2 = sqrt( (double)nSubDivA[i]*Ratio12 );

			int nDivB = min(nSubDivA[i], (int)round(N2));

		//	printf("Partition [%i] dir2 into %i\n", i, nDivB);

			Partition* RegionsDir2 = new Partition[nDivB];

			int* nSubDivB = new int[nDivB];					//Number of sub divisions in each partition
			double* nNodesDivB = new double[nDivB];			//Number of fluid nodes in each partition

			long long nNodesDivA_Rounded = WeightMode ? _Lattice->CountOccurences(_NodeValue, _NodeValueCount, CountBuf, RegionsDir1[i].Region) : RegionsDir1[i].Region.CountCells();

			for(int c=0; c<nDivB; c++){

				nSubDivB[c] = nSubDivA[i]/nDivB + ( (c < nSubDivA[i]%nDivB) ? 1 : 0 );

				nNodesDivB[c] = (double)nSubDivB[c] * ((double)nNodesDivA_Rounded/(double)nSubDivA[i]);

			//	printf("Partition [%i][%i] into %i with %f nodes\n", i, c, nSubDivB[c], nNodesDivB[c]);

			}
			
			r = DividePartition(nDivB, Dirs[1], RegionsDir1[i], nNodesDivB, RegionsDir2, 0);

			if(!r){
				delete[] nSubDivB;
				delete[] nNodesDivB;
				delete[] RegionsDir2;
				delete[] nSubDivA;
				delete[] nNodesDivA;
				delete[] RegionsDir1;
				delete[] CountBuf;
				return false;
			}

			//Partition in direction 3

			for(int i2=0; i2<nDivB; i2++){

				double* nNodesDivC = new double[nSubDivB[i2]];

				long long nNodesDivB_Rounded = WeightMode ? _Lattice->CountOccurences(_NodeValue, _NodeValueCount, CountBuf, RegionsDir2[i2].Region) : RegionsDir2[i2].Region.CountCells();

				for(int c=0; c<nSubDivB[i2]; c++)
					nNodesDivC[c] = (double)nNodesDivB_Rounded / (double)nSubDivB[i2];

				r = DividePartition(nSubDivB[i2], Dirs[2], RegionsDir2[i2], nNodesDivC, Partitions, OutIndex);

				delete[] nNodesDivC;

				if(!r){
					delete[] nSubDivB;
					delete[] nNodesDivB;
					delete[] RegionsDir2;
					delete[] nSubDivA;
					delete[] nNodesDivA;
					delete[] RegionsDir1;
					delete[] CountBuf;
					return false;
				}

				OutIndex += nSubDivB[i2];

			}

			delete[] nSubDivB;
			delete[] nNodesDivB;
			delete[] RegionsDir2;

		}

		delete[] nSubDivA;
		delete[] nNodesDivA;
		delete[] RegionsDir1;
		delete[] CountBuf;

		return true;
	}

	bool DividePartition2(int nDecompose, Direction Dir, Partition PartitionIn, long long* NodeDist, double* nNodesPerDivision, Partition* PartitionsOut, int OutStartIndex){

		int d  = (Dir==Dir_X) ? 0 : ( (Dir==Dir_Y) ? 1 : 2 );	//Direction of decomposition

		//Decompose in 1D

		long long* NodeDistCumulative = new long long[arrSz[d]];		//Cumulative number of non solids

		for(int c=0; c<arrSz[d]; c++){

			NodeDistCumulative[c] = (c > 0) ? (NodeDistCumulative[c-1] + NodeDist[c]) : NodeDist[0];
		
		}

		//Divide

		LatticeBlock* Block = new LatticeBlock[nDecompose];

			//Find precise range
	
		long long x0 = 0;
		double x0_ = 0;

		long long x1 = 0;
		double x1_ = 0;

		double NSolidsPreceding = 0;

		for(int i=0; i<nDecompose; i++){

			while(true){

				double N0 = (x1==0) ? 0 : NodeDistCumulative[x1-1] - NSolidsPreceding;
				double N1 = NodeDistCumulative[x1] - NSolidsPreceding;

				if(i==nDecompose-1){

					x1=arrSz[d];
					x1_ = 0;
					break;

				}else if(N0<nNodesPerDivision[i] && N1>=nNodesPerDivision[i]){

					x1_ = (nNodesPerDivision[i] - N0) / (double)NodeDist[x1];
					break;

				}

				x1++;
			}

			NSolidsPreceding = (x1==0 ? 0 : NodeDistCumulative[x1-1]) + x1_*NodeDist[x1];

			Block[i].c0 = x0;
			Block[i].c0_ = x0_;
			Block[i].c1 = x1;
			Block[i].c1_ = x1_;
			Block[i].n_d = ((double)x1 + x1_) - ((double)x0 + x0_);

			x0 = x1;
			x0_ = x1_;

		}

			//Set rounded range

		if(nDecompose>1){

			bool ret = MinimiseError(Block, nDecompose, NodeDistCumulative, arrSz[d], nNodesPerDivision);

			if(!ret){
				delete[] NodeDistCumulative;
				delete[] Block;
				return false;
			}

		}else{

			Block[0].c0Rounded = 0;
			Block[0].c1Rounded = arrSz[d];

		}

		//Output domains

		for(int i=0; i<nDecompose; i++){

			PartitionsOut[OutStartIndex + i].Region.x0 = (Dir==Dir_X) ? Block[i].c0Rounded : PartitionIn.Region.x0;
			PartitionsOut[OutStartIndex + i].Region.x1 = (Dir==Dir_X) ? Block[i].c1Rounded : PartitionIn.Region.x1;
			PartitionsOut[OutStartIndex + i].Region.y0 = (Dir==Dir_Y) ? Block[i].c0Rounded : PartitionIn.Region.y0;
			PartitionsOut[OutStartIndex + i].Region.y1 = (Dir==Dir_Y) ? Block[i].c1Rounded : PartitionIn.Region.y1;
			PartitionsOut[OutStartIndex + i].Region.z0 = (Dir==Dir_Z) ? Block[i].c0Rounded : PartitionIn.Region.z0;
			PartitionsOut[OutStartIndex + i].Region.z1 = (Dir==Dir_Z) ? Block[i].c1Rounded : PartitionIn.Region.z1;

			PartitionsOut[OutStartIndex + i].NodeCount = NodeDistCumulative[Block[i].c1Rounded-1] - ((Block[i].c0Rounded == 0) ? 0 : (NodeDistCumulative[Block[i].c0Rounded-1]));

		//	if(MPIThreadId==0)
		//		printf("Partition in range (%s) with %lli nodes\n", PartitionsOut[OutStartIndex + i].Region.ToStr(), PartitionsOut[OutStartIndex + i].NodeCount);
		
		}

		delete[] NodeDistCumulative;
		delete[] Block;

		return true;
	}

	bool DividePartition(int nDecompose, Direction Dir, Partition PartitionIn, double* nNodesPerDivision, Partition* PartitionsOut, int OutStartIndex){

		Grid<T>& Lattice = *_Lattice;

		int d  = (Dir==Dir_X) ? 0 : ( (Dir==Dir_Y) ? 1 : 2 );	//Direction of decomposition
		int d1 = (Dir==Dir_X) ? 2 : ( (Dir==Dir_Y) ? 2 : 1 );	//Other directions
		int d0 = (Dir==Dir_X) ? 1 : ( (Dir==Dir_Y) ? 0 : 0 );

		long long Region0[3];
		long long RegionSz[3];

		Region0[0] = PartitionIn.Region.x0;
		Region0[1] = PartitionIn.Region.y0;
		Region0[2] = PartitionIn.Region.z0;

		RegionSz[0] = PartitionIn.Region.SizeX();
		RegionSz[1] = PartitionIn.Region.SizeY();
		RegionSz[2] = PartitionIn.Region.SizeZ();

		//Decompose in 1D

		long long* NonSolidsSlices = new long long[arrSz[d]];			//Count number of non solids in each slice
		long long* NonSolidsCumulative = new long long[arrSz[d]];		//Cumulative number of non solids

		int Count[3];

		for(Count[d]=0; Count[d]!=RegionSz[d]; Count[d]++){	//Decomposition direction

			//Count nodes in slice

			NonSolidsSlices[Count[d]] = 0;

			for(Count[d1]=0; Count[d1]!=RegionSz[d1]; Count[d1]++)
			for(Count[d0]=0; Count[d0]!=RegionSz[d0]; Count[d0]++){

				if(!WeightMode || IsActiveNode( Lattice(Region0[0] + Count[0], Region0[1] + Count[1], Region0[2] + Count[2]) ))
					NonSolidsSlices[Count[d]]++;
			}

			//Cumulative

			if(Count[d]!=0){
				NonSolidsCumulative[Count[d]] = NonSolidsCumulative[Count[d]-1] + NonSolidsSlices[Count[d]];
			}else{
				NonSolidsCumulative[Count[d]] = NonSolidsSlices[0];
			}

		}

	/*	if(MPIThreadId==0 && d==2)
		for(int c=0; c<arrSz[d]; c++){
			printf("%i, ", NonSolidsSlices[c]);

			if((c+1)%10 == 0)
				printf("\n");
		}	//*/

		//Divide

		LatticeBlock* Block = new LatticeBlock[nDecompose];

			//Find precise range
	
		long long x0 = 0;
		double x0_ = 0;

		long long x1 = 0;
		double x1_ = 0;

		double NSolidsPreceding = 0;

		for(int i=0; i<nDecompose; i++){

			while(true){

				double N0 = (x1==0) ? 0 : NonSolidsCumulative[x1-1] - NSolidsPreceding;
				double N1 = NonSolidsCumulative[x1] - NSolidsPreceding;

				if(i==nDecompose-1){

					x1=arrSz[d];
					x1_ = 0;
					break;

				}else if(N0<nNodesPerDivision[i] && N1>=nNodesPerDivision[i]){

					x1_ = (nNodesPerDivision[i] - N0) / (double)NonSolidsSlices[x1];
					break;

				}

				x1++;
			}

			NSolidsPreceding = (x1==0 ? 0 : NonSolidsCumulative[x1-1]) + x1_*NonSolidsSlices[x1];

			Block[i].c0 = x0;
			Block[i].c0_ = x0_;
			Block[i].c1 = x1;
			Block[i].c1_ = x1_;
			Block[i].n_d = ((double)x1 + x1_) - ((double)x0 + x0_);

			x0 = x1;
			x0_ = x1_;

		}

			//Set rounded range

		if(nDecompose>1){

			bool ret = MinimiseError(Block, nDecompose, NonSolidsCumulative, arrSz[d], nNodesPerDivision);

			if(!ret){
				delete[] NonSolidsSlices;
				delete[] NonSolidsCumulative;
				delete[] Block;
				return false;
			}

		}else{

			Block[0].c0Rounded = 0;
			Block[0].c1Rounded = arrSz[d];

		}

		//Output domains

		for(int i=0; i<nDecompose; i++){

			PartitionsOut[OutStartIndex + i].Region.x0 = (Dir==Dir_X) ? Block[i].c0Rounded : PartitionIn.Region.x0;
			PartitionsOut[OutStartIndex + i].Region.x1 = (Dir==Dir_X) ? Block[i].c1Rounded : PartitionIn.Region.x1;
			PartitionsOut[OutStartIndex + i].Region.y0 = (Dir==Dir_Y) ? Block[i].c0Rounded : PartitionIn.Region.y0;
			PartitionsOut[OutStartIndex + i].Region.y1 = (Dir==Dir_Y) ? Block[i].c1Rounded : PartitionIn.Region.y1;
			PartitionsOut[OutStartIndex + i].Region.z0 = (Dir==Dir_Z) ? Block[i].c0Rounded : PartitionIn.Region.z0;
			PartitionsOut[OutStartIndex + i].Region.z1 = (Dir==Dir_Z) ? Block[i].c1Rounded : PartitionIn.Region.z1;

			PartitionsOut[OutStartIndex + i].NodeCount = NonSolidsCumulative[Block[i].c1Rounded-1] - ((Block[i].c0Rounded == 0) ? 0 : (NonSolidsCumulative[Block[i].c0Rounded-1]));
		
		}

		delete[] NonSolidsSlices;
		delete[] NonSolidsCumulative;
		delete[] Block;

		return true;
	}

	bool MinimiseError(LatticeBlock* Plane, int nBlocks, long long* nCumulative, long long NLattice, double* nNodesPerDivision){

		//Obtain minimum error combination of rounded borders

		for(int i=0; i<nBlocks; i++){

			//Combinations of rounded borders for block

			Plane[i].nRoundActive = 0;

			for(int BState=0; BState<4; BState++){		//Combination of borders, 4 permutations

				Plane[i].LowerRound = (BState&1)!=0;	//Bit 0
				Plane[i].UpperRound = (BState&2)!=0;	//Bit 1

				if(i==0 && Plane[i].LowerRound){				//Bottom plane and bottom border rounded up

					Plane[i].RoundActive[BState] = false;	//Not an applicable case
					continue;

				}else if(i==nBlocks-1 && Plane[i].UpperRound){	//Top plane and top border rounded up

					Plane[i].RoundActive[BState] = false;	//Not an applicable case
					continue;
				}

				//Compute error

				long long zIndex0 = Plane[i].LowerRound ? Plane[i].c0+(long long)ceil(Plane[i].c0_) : Plane[i].c0+(long long)floor(Plane[i].c0_);
				long long zIndex1 = Plane[i].UpperRound ? Plane[i].c1+(long long)ceil(Plane[i].c1_) : Plane[i].c1+(long long)floor(Plane[i].c1_);

				long long nCumlN = (i==0) ? 0 : nCumulative[zIndex0-1];
				long long nCumlP = nCumulative[zIndex1-1];

				double delta = (nCumlP - nCumlN) - nNodesPerDivision[i];

				Plane[i].RoundError[BState] = (delta*delta)/(nNodesPerDivision[i]*nNodesPerDivision[i]);	//Error

				Plane[i].RoundActive[BState] = true;
				Plane[i].nRoundActive++;

			}

		}

		//Remove highest error combinations until only lowest error combination remains

		while(true){

			double HighestError = -1;
			int HighestErrorAddr[2];	// { BlockId, BState }

			bool End = true;

			for(int i=0; i<nBlocks; i++){

				if(Plane[i].nRoundActive <= 1)
					continue;

				End = false;

				for(int BState=0; BState<4; BState++){

					if(Plane[i].RoundActive[BState] && Plane[i].RoundError[BState] > HighestError){

						HighestError = Plane[i].RoundError[BState];

						HighestErrorAddr[0] = i;
						HighestErrorAddr[1] = BState;
					}

				}
			}

	/*
			for(int i=0; i<nBlocks; i++){
				bool UpperConstrained = ( Plane[i].nRoundActive == 1 || ( Plane[i].nRoundActive == 2 && ((Plane[i].RoundActive[0]&&Plane[i].RoundActive[1]) || (Plane[i].RoundActive[2]&&Plane[i].RoundActive[3])) ) );		//Are upper and lower boundaries now fixed
				bool LowerConstrained = ( Plane[i].nRoundActive == 1 || ( Plane[i].nRoundActive == 2 && ((Plane[i].RoundActive[0]&&Plane[i].RoundActive[2]) || (Plane[i].RoundActive[1]&&Plane[i].RoundActive[3])) ) );
				
				printf("Partition %i (%i active)", i, Plane[i].nRoundActive);

				if(UpperConstrained)
					printf(" (U constrained)");
				if(LowerConstrained)
					printf(" (L constrained)");

				printf("\n");

				for(int BState=0; BState<4; BState++){
					if(Plane[i].RoundActive[BState])
						printf("  BState [%i] error = %e\n", BState, Plane[i].RoundError[BState]);
				}
			}		//	*/

			if(End){
		//		printf("\nReached minimal error\n");
				break;
			}

			//Remove highest error

			int BlockId = HighestErrorAddr[0];
			int BState  = HighestErrorAddr[1];

		//	printf("\nHighest error at block [%i], BState [%i]\n\n", BlockId, BState);

			Plane[BlockId].RoundActive[BState] = false;
			Plane[BlockId].nRoundActive--;

			int nActive = Plane[BlockId].nRoundActive;
			bool BStates[4];
			memcpy(BStates, Plane[BlockId].RoundActive, sizeof(bool)*4);

			bool UpperConstrained = ( nActive == 1 || ( nActive == 2 && ((BStates[0]&&BStates[1]) || (BStates[2]&&BStates[3])) ) );		//Are upper and lower boundaries now fixed
			bool LowerConstrained = ( nActive == 1 || ( nActive == 2 && ((BStates[0]&&BStates[2]) || (BStates[1]&&BStates[3])) ) );

			bool UpperRound = (BStates[2] || BStates[3]);	//If constrained, boundary up or down
			bool LowerRound = (BStates[1] || BStates[3]);

			for(int i=BlockId+1; UpperConstrained && i<nBlocks; i++){			//Cascade constraint upwards

			//	printf("Lower boundary constrained (%i) of partition [%i]\n\n", (int)(UpperRound ? 1 : 0), i);

				bool UpperConstrained0 = ( Plane[i].nRoundActive == 1 || ( Plane[i].nRoundActive == 2 && ((Plane[i].RoundActive[0]&&Plane[i].RoundActive[1]) || (Plane[i].RoundActive[2]&&Plane[i].RoundActive[3])) ) );
				
				if(UpperRound){
					if(Plane[i].RoundActive[0]) Plane[i].nRoundActive--;
					if(Plane[i].RoundActive[2]) Plane[i].nRoundActive--;
					Plane[i].RoundActive[0] = false;
					Plane[i].RoundActive[2] = false;
				}else{
					if(Plane[i].RoundActive[1]) Plane[i].nRoundActive--;
					if(Plane[i].RoundActive[3]) Plane[i].nRoundActive--;
					Plane[i].RoundActive[1] = false;
					Plane[i].RoundActive[3] = false;
				}

				UpperConstrained = ( Plane[i].nRoundActive == 1 || ( Plane[i].nRoundActive == 2 && ((Plane[i].RoundActive[0]&&Plane[i].RoundActive[1]) || (Plane[i].RoundActive[2]&&Plane[i].RoundActive[3])) ) );
				
				if(UpperConstrained == UpperConstrained0)	//Constraint has propagated
					break;

				UpperRound = (Plane[i].RoundActive[2] || Plane[i].RoundActive[3]);	//Upper border is now constrained as well

			//	printf("Upper boundary now constrained too (%i)\n\n", (int)(UpperRound ? 1 : 0));

			}

			for(int i=BlockId-1; LowerConstrained && i>=0; i--){			//Cascade constraint downwards

			//	printf("Upper boundary constrained (%i) of partition [%i]\n\n", (int)(LowerRound ? 1 : 0), i);

				bool LowerConstrained0 = ( Plane[i].nRoundActive == 1 || ( Plane[i].nRoundActive == 2 && ((Plane[i].RoundActive[0]&&Plane[i].RoundActive[2]) || (Plane[i].RoundActive[1]&&Plane[i].RoundActive[3])) ) );
				
				if(LowerRound){
					if(Plane[i].RoundActive[0]) Plane[i].nRoundActive--;
					if(Plane[i].RoundActive[1]) Plane[i].nRoundActive--;
					Plane[i].RoundActive[0] = false;
					Plane[i].RoundActive[1] = false;
				}else{
					if(Plane[i].RoundActive[2]) Plane[i].nRoundActive--;
					if(Plane[i].RoundActive[3]) Plane[i].nRoundActive--;
					Plane[i].RoundActive[2] = false;
					Plane[i].RoundActive[3] = false;
				}

				LowerConstrained = ( Plane[i].nRoundActive == 1 || ( Plane[i].nRoundActive == 2 && ((Plane[i].RoundActive[0]&&Plane[i].RoundActive[2]) || (Plane[i].RoundActive[1]&&Plane[i].RoundActive[3])) ) );
				
				if(LowerConstrained == LowerConstrained0)	//Constraint has propagated
					break;

				LowerRound = (Plane[i].RoundActive[1] || Plane[i].RoundActive[3]);		//Lower border is now constrained as well

			//	printf("Lower boundary now constrained too (%i)\n\n", (int)(LowerRound ? 1 : 0));

			}

		}

		//Set minimum error combination

		for(int i=0; i<nBlocks; i++){

			if(Plane[i].nRoundActive != 1){		//Partition has not ended up with one single border combination
				return false;
			}

			for(int BState=0; BState<4; BState++){

				if(Plane[i].RoundActive[BState]){

					Plane[i].LowerRound = (BState&1)!=0;	//Bit 0
					Plane[i].UpperRound = (BState&2)!=0;	//Bit 1

					long long zIndex0 = Plane[i].LowerRound ? Plane[i].c0+(long long)ceil(Plane[i].c0_) : Plane[i].c0+(long long)floor(Plane[i].c0_);
					long long zIndex1 = Plane[i].UpperRound ? Plane[i].c1+(long long)ceil(Plane[i].c1_) : Plane[i].c1+(long long)floor(Plane[i].c1_);

					Plane[i].c0Rounded = zIndex0;
					Plane[i].c1Rounded = zIndex1;

					long long nCumlN = (i==0) ? 0 : nCumulative[zIndex0-1];
					long long nCumlP = nCumulative[zIndex1-1];

					if(nCumlP-nCumlN == 0){

						return false;			//Combination leads to 0 nodes in partition

					}

			//		printf("Partition [%i] from z = %i -> %i (%.2f -> %.2f)\n", i, Plane[i].c0Rounded, Plane[i].c1Rounded, (float)(Plane[i].c0 + Plane[i].c0_), (float)(Plane[i].c1 + Plane[i].c1_));
				}

			}
		}

		return true;
	}

};
