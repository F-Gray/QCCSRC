#pragma once

#include <mpi.h>
#include <vector>
#include <stdio.h>

#ifdef _WIN32
#include <Windows.h>
#endif

#ifndef _WIN32
#include <unistd.h>
#include <sys/sysinfo.h>
#endif

using namespace std;

//Memory

void GetSystemMemory(size_t* MemTot, size_t* MemFree){

#ifdef _WIN32

	//Windows

	MEMORYSTATUSEX statex;

	statex.dwLength = sizeof(statex);

	GlobalMemoryStatusEx(&statex);

	if(MemTot != 0)
		*MemTot = statex.ullTotalPhys;

	if(MemFree != 0)
		*MemFree = statex.ullAvailPhys;

#else

	struct sysinfo si;

	sysinfo(&si);

	if(MemTot != 0)
		*MemTot = si.totalram;

	if(MemFree != 0)
		*MemFree = si.freeram;

#endif

}

//System topology

struct GPUInfo;
struct HostStruct;
struct HostThreadInfo;

struct HostThreadInfo{

	int ThreadId;			//MPI thread index

	char HostName[256];		//Host name
	size_t HostMemTot;		//Total memory on host

	HostStruct* Host;		//Host of this thread
	
	int  nGPUs;				//Number of GPUs on this host

	bool HasGPU;			//Associated with GPU
	int GPUIndex;			//Index of thread GPU
	GPUInfo* GPU;			//GPU thread associated with this thread

	HostThreadInfo(){
		memset(this, 0, sizeof(HostThreadInfo));
	}

};

struct GPUInfo{

	char DeviceName[256];		//Device name
	
	int DeviceIndex;			//Index of device on host

	cudaDeviceProp Properties;

	char HostName[256];			//Host name
	HostStruct* Host;			//Ptr to host

	bool HasMPIThread;			//GPU associated with MPI thread
	int MPIThreadId;			//MPI thread associated with device
	HostThreadInfo* MPIThread;	//Ptr to MPI thread

	GPUInfo(){
		memset(this, 0, sizeof(GPUInfo));
	}

};

struct HostStruct{

	char HostName[256];		//Host name

	size_t TotalMem;		//Total system memory of host

	int nThreads;			//Number of MPI threads running on host
	int nGPUs;				//Number of GPUs on this host

	int* ThreadIds;			//Id of threads running on this host
	HostThreadInfo** Threads;	//Array of ptrs to MPI threads on host

	GPUInfo** GPUs;			//Ptr to GPUs on this host

	HostStruct(){
		memset(this, 0, sizeof(HostStruct));
	}

	~HostStruct(){
		if(ThreadIds!=0) delete[] ThreadIds;
		if(Threads!=0)   delete[] Threads;
		if(GPUs!=0)      delete[] GPUs;
	}

};

enum GPUAssignment{
	GPUAssignment_None,				//Does not associate GPUs with MPI threads
	GPUAssignment_Sequential,		//Assign sequential MPI threads sequential GPUs
	GPUAssignment_Default,			//Obtain default GPU using cudaGetDevice for each thread
	GPUAssignment_Auto				//Obtains default GPUs using cudaGetDevice, but reassigns any duplicates sequentially
};

/*!
	Class maps the topology of an MPI system, finding servers, MPI threads and GPUs.
*/

class TopologyMPI{

public:

	HostThreadInfo* Threads;
	HostStruct* Hosts;
	GPUInfo* GPUs;

	int ThreadCount;
	int HostCount;
	int GPUCount;

	MPI_Comm Communicator;

private:

	GPUAssignment _AssignmentMode;

	void Init(){
		memset(this, 0, sizeof(TopologyMPI));
	}

	void CountHostsGPUs(int* nHosts, int* nGPUs, bool FillLists){

		int CountHosts = 0;
		int CountGPUs = 0;

		int HostInd = 0;
		int GPUInd = 0;

		for(int i=0; i<ThreadCount; i++){	//For each thread

			bool Found = false;

			for(int i2=0; i2<i; i2++){		//Check preceding threads

				if(strcmp(Threads[i].HostName, Threads[i2].HostName) == 0){		//Matches
					Found = true;
					break;
				}
			}

			if(!Found){		//New host

				CountHosts++;
				CountGPUs += Threads[i].nGPUs;

				if(FillLists){

					//Add host

					strcpy(Hosts[HostInd].HostName, Threads[i].HostName);	//Name of host
					Hosts[HostInd].TotalMem = Threads[i].HostMemTot;		//Total memory on host
					Hosts[HostInd].nGPUs = Threads[i].nGPUs;				//Number of GPUs on host
					HostInd++;

					//Add host GPUs

					for(int c=0; c<Threads[i].nGPUs; c++){

						strcpy(GPUs[GPUInd].HostName, Threads[i].HostName);	//Name of host
						GPUs[GPUInd].DeviceIndex = c;						//Local device index
						GPUInd++;

					}

				}
			}
		}

		*nHosts = CountHosts;
		*nGPUs = CountGPUs;

	}

	void FillLists(){		//Collect association lists of threads and GPUs on each host

		for(int c=0; c<HostCount; c++){

			vector<HostThreadInfo*> ThreadList;
			vector<GPUInfo*> GPUList;

			HostStruct* Host = &(Hosts[c]);

			//Find threads and GPUs

			for(int i=0; i<ThreadCount; i++){
				if(strcmp(Threads[i].HostName, Host->HostName) == 0)	//Matches
					ThreadList.push_back(&(Threads[i]));
			}

			for(int i=0; i<GPUCount; i++){
				if(strcmp(GPUs[i].HostName, Host->HostName) == 0)		//Matches
					GPUList.push_back(&(GPUs[i]));
			}

			//Create arrays

			Host->Threads = new HostThreadInfo*[ThreadList.size()];
			Host->ThreadIds = new int[ThreadList.size()];
			
			Host->GPUs = new GPUInfo*[GPUList.size()];

			Host->nThreads = (int)ThreadList.size();
			Host->nGPUs = (int)GPUList.size();
			
			//Fill arrays

			for(int i=0; i<ThreadList.size(); i++){

				Host->Threads[i] = ThreadList[i];
				Host->ThreadIds[i] = (ThreadList[i])->ThreadId;

				(ThreadList[i])->Host = Host;	//Thread ptr to host
			}

			for(int i=0; i<GPUList.size(); i++){

				Host->GPUs[i] = GPUList[i];

				(GPUList[i])->Host = Host;		//GPU ptr to host
			}

		}

	}

	void AssociateGPUs(int tId, GPUAssignment AssignmentMode){	//Associate GPUs with local MPI thread

		if(AssignmentMode == GPUAssignment_None)
			return;

		for(int c=0; c<HostCount; c++){

			if(AssignmentMode == GPUAssignment_Sequential){

				//Associate GPUs with MPI threads sequentially

				for(int i=0; i<Hosts[c].nGPUs; i++){

					if(i<Hosts[c].nThreads){	//Associate thread with GPU

						(Hosts[c].GPUs[i])->HasMPIThread = true;
						(Hosts[c].GPUs[i])->MPIThreadId = Hosts[c].ThreadIds[i];
						(Hosts[c].GPUs[i])->MPIThread = Hosts[c].Threads[i];

						(Hosts[c].Threads[i])->HasGPU = true;
						(Hosts[c].Threads[i])->GPUIndex = Hosts[c].GPUs[i]->DeviceIndex;
						(Hosts[c].Threads[i])->GPU = Hosts[c].GPUs[i];
					}
				}

			}else if(AssignmentMode == GPUAssignment_Default || AssignmentMode == GPUAssignment_Auto){
			
				//Use default assigned GPUs

				for(int i=0; i<Hosts[c].nThreads; i++){

					int n = (Hosts[c].Threads[i])->GPUIndex;	//Default MPI assigned GPU

					for(int i2=0; i2<Hosts[c].nGPUs; i2++){

						if(Hosts[c].GPUs[i2]->DeviceIndex == n && !Hosts[c].GPUs[i2]->HasMPIThread){	//Associate
						
							(Hosts[c].GPUs[i2])->HasMPIThread = true;
							(Hosts[c].GPUs[i2])->MPIThreadId = Hosts[c].ThreadIds[i];
							(Hosts[c].GPUs[i2])->MPIThread = Hosts[c].Threads[i];
						
							(Hosts[c].Threads[i])->HasGPU = true;
							(Hosts[c].Threads[i])->GPU = Hosts[c].GPUs[i2];

							break;
						}
					}

				}

				if(AssignmentMode != GPUAssignment_Auto)
					continue;

				//Associate any remaining GPUs sequentially

				for(int i=0; i<Hosts[c].nGPUs; i++){

					if(Hosts[c].GPUs[i]->HasMPIThread)
						continue;

					//Find free thread

					for(int i2=0; i2<Hosts[c].nThreads; i2++){

						if(!Hosts[c].Threads[i2]->HasGPU){		//Associate
							
							(Hosts[c].GPUs[i])->HasMPIThread = true;
							(Hosts[c].GPUs[i])->MPIThreadId = Hosts[c].ThreadIds[i2];
							(Hosts[c].GPUs[i])->MPIThread = Hosts[c].Threads[i2];

							(Hosts[c].Threads[i2])->HasGPU = true;
							(Hosts[c].Threads[i2])->GPUIndex = Hosts[c].GPUs[i]->DeviceIndex;
							(Hosts[c].Threads[i2])->GPU = Hosts[c].GPUs[i];

							break;
						}
					}

				}
			
			}

		}

	}

	bool ObtainGPUProperties(MPI_Comm Comm){

		int ThreadId;
		MPI_Comm_rank(Comm, &ThreadId);		//Get thread ID (index from 0 -> nThreads-1)

		for(int i=0; i<GPUCount; i++){

			cudaError_t cudaStatus;

			int GPUThread = GPUs[i].Host->ThreadIds[0];		//First thread on GPU's host distributes properties

			if(ThreadId == GPUThread){
				cudaStatus = cudaGetDeviceProperties(&(GPUs[i].Properties), GPUs[i].DeviceIndex);
			}

			MPI_Bcast(&cudaStatus, sizeof(cudaError_t), MPI_CHAR, GPUThread, Comm);
			MPI_Bcast(&(GPUs[i].Properties), sizeof(cudaDeviceProp), MPI_CHAR, GPUThread, Comm);

			if(cudaStatus!=cudaSuccess && !MPITestSetup){
				if(ThreadId==0)
					printf("Unable to obtain properties of GPU [%i] on '%s'\n", GPUs[i].DeviceIndex, GPUs[i].HostName);
				return false;
			}
			
		}

		return true;
	}

public:

	TopologyMPI(){

		Init();

	}

	bool ObtainTopology(MPI_Comm Comm){
	
		return ObtainTopology(Comm, GPUAssignment_Auto);
	}

	bool ObtainTopology(MPI_Comm Comm, GPUAssignment AssignmentMode){

		Communicator = Comm;

		MPI_Comm_size(Comm, &ThreadCount);	//Get number of threads running

		int ThreadId;
		MPI_Comm_rank(Comm, &ThreadId);		//Get thread ID (index from 0 -> nThreads-1)

		int Error = 0;


		//Local thread's data

		HostThreadInfo tInfo;		

		tInfo.ThreadId = ThreadId;	//Thread Id

		int ret = gethostname(tInfo.HostName, 256);	//Obtain host name

		if(ret!=0)
			Error = 1;

		if(MPITestSetup){

			const int nHosts = MPITestCount;

			int HostInd = ThreadId/(ThreadCount/nHosts);

			tInfo.HostName[0] += HostInd;
		}

		//Host memory

		GetSystemMemory(&tInfo.HostMemTot, 0);

		//Obtain number of GPUs on this board

		cudaError_t cudaStatus = cudaGetDeviceCount(&(tInfo.nGPUs));	//Obtain number of local GPUs

		if(cudaStatus!=cudaSuccess || cudaStatus==cudaErrorNoDevice || tInfo.nGPUs < 1){	//Unable to find number of GPUs or no GPUs on this host

			tInfo.nGPUs = 0;

			if(cudaStatus!=cudaSuccess && !MPITestSetup)
				Error = 2;

			if(MPITestSetup)
				tInfo.nGPUs = 8;
		}else{

			float* HMem;
			cudaMallocHost(&HMem, 128);		//First have to run a cuda routine to ensure default assignment (bug in cuda presumably)
			cudaFreeHost(HMem);

			if(cudaGetDevice(&(tInfo.GPUIndex)) != cudaSuccess)		//Default GPU assigned to this MPI thread
				Error = 3;
		}


		//Check for errors

		int ErrorMax;
		MPI_Allreduce(&Error, &ErrorMax, 1, MPI_INT, MPI_MAX, Comm);

		if(ErrorMax > 0){
			if(ThreadId==0)
				printf("Unable to obtain number of available GPUs\n");
			return false;
		}


		//Thread data

		Threads = new HostThreadInfo[ThreadCount];

		MPI_Allgather(&tInfo, sizeof(HostThreadInfo), MPI_BYTE, Threads, sizeof(HostThreadInfo), MPI_BYTE, Comm);

		CountHostsGPUs(&HostCount, &GPUCount, false);		//Count total number of unique hosts and GPUs

		Hosts = new HostStruct[HostCount];
		GPUs = new GPUInfo[GPUCount];

		memset(Hosts, 0, sizeof(HostStruct)*HostCount);
		memset(GPUs, 0, sizeof(GPUInfo)*GPUCount);
		
		CountHostsGPUs(&HostCount, &GPUCount, true);		//Fill new arrays

		FillLists();										//Fills sublists of threads, GPUs and hosts

		AssociateGPUs(ThreadId, AssignmentMode);						//Pair GPUs with MPI threads

		bool r = ObtainGPUProperties(Comm);

		return r;
	}

	void PrintTopology(){

		//Write out hosts

		for(int i=0; i<HostCount; i++){

			printf("Host '%s' (Total Mem: %.1fGB) with %i GPUs and %i MPI threads\n", Hosts[i].HostName, Hosts[i].TotalMem/(1024.0*1024.0*1024.0), Hosts[i].nGPUs, Hosts[i].nThreads);
			printf("   Threads: ");

			for(int c=0; c<Hosts[i].nThreads; c++){
				if(c==Hosts[i].nThreads-1) printf("%i", Hosts[i].ThreadIds[c]);
				else					   printf("%i, ", Hosts[i].ThreadIds[c]);
			}

			printf("\n\n");
		}

		printf("Total number of GPUs: %i\n\n", GPUCount);

		//GPUs

		for(int i=0; i<GPUCount; i++){

			if(GPUs[i].HasMPIThread)
				printf("GPU [%i] on '%s' thread %i:\n", GPUs[i].DeviceIndex, GPUs[i].Host->HostName, GPUs[i].MPIThreadId);
			else
				printf("GPU [%i] on '%s' (no thread assigned):\n", GPUs[i].DeviceIndex, GPUs[i].Host->HostName);

			printf("   Name: %s\n   Mem:  %.1fGB\n", GPUs[i].Properties.name, GPUs[i].Properties.totalGlobalMem/(1024.0*1024.0*1024.0));

		}

		/*//Threads

		for(int i=0; i<ThreadCount; i++){

			printf("Thread [%i] on board '%s' (%s GPU)\n", i, Threads[i].HostName, (Threads[i].HasGPU ? "has" : "no"));

		}	*/

	}

	~TopologyMPI(){

		if(Threads!=0) delete[] Threads;
		if(Hosts!=0)   delete[] Hosts;
		if(GPUs!=0)    delete[] GPUs;

	}

};



