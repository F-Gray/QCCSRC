#pragma once

#include <mpi.h>
#include <vector>

#ifdef _WIN32
#include <Windows.h>
#endif

#ifndef _WIN32
#include <unistd.h>
#endif

using namespace std;

struct HostStruct;
struct HostThreadInfo;

struct HostThreadInfo{

	int ThreadId;			//MPI thread index

	char HostName[256];		//Host name

	HostStruct* Host;		//Host of this thread

	HostThreadInfo(){
		memset(this, 0, sizeof(HostThreadInfo));
	}

};

struct HostStruct{

	char HostName[256];		//Host name

	int nThreads;			//Number of MPI threads running on host
	
	int* ThreadIds;			//Id of threads running on this host
	HostThreadInfo** Threads;	//Array of ptrs to MPI threads on host

	HostStruct(){
		memset(this, 0, sizeof(HostStruct));
	}

	~HostStruct(){
		if(ThreadIds!=0) delete[] ThreadIds;
		if(Threads!=0)   delete[] Threads;
	}

};


/*!
	Class maps the topology of an MPI system, finding servers and MPI threads
*/

class TopologyMPI{

public:

	HostThreadInfo* Threads;
	HostStruct* Hosts;

	int ThreadCount;
	int HostCount;

	MPI_Comm Communicator;

private:

	void Init(){
		memset(this, 0, sizeof(TopologyMPI));
	}

	void CountHosts(int* nHosts, bool FillLists){

		int CountHosts = 0;

		int HostInd = 0;

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

				if(FillLists){

					//Add host

					strcpy(Hosts[HostInd].HostName, Threads[i].HostName);	//Name of host
					HostInd++;

				}
			}
		}

		*nHosts = CountHosts;

	}

	void FillLists(){		//Collect association lists of threads and GPUs on each host

		for(int c=0; c<HostCount; c++){

			vector<HostThreadInfo*> ThreadList;

			HostStruct* Host = &(Hosts[c]);

			//Find threads

			for(int i=0; i<ThreadCount; i++){
				if(strcmp(Threads[i].HostName, Host->HostName) == 0)	//Matches
					ThreadList.push_back(&(Threads[i]));
			}

			//Create arrays

			Host->Threads = new HostThreadInfo*[ThreadList.size()];
			Host->ThreadIds = new int[ThreadList.size()];
			
			Host->nThreads = (int)ThreadList.size();
			
			//Fill arrays

			for(int i=0; i<ThreadList.size(); i++){

				Host->Threads[i] = ThreadList[i];
				Host->ThreadIds[i] = (ThreadList[i])->ThreadId;

				(ThreadList[i])->Host = Host;	//Thread ptr to host
			}

		}

	}

public:

	TopologyMPI(){

		Init();

	}

	bool ObtainTopology(MPI_Comm Comm){

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

#ifdef MPITestSetup

		if(MPITestSetup){

			int nHosts = min(ThreadCount, MPITestCount);

			int HostInd = ThreadId/(ThreadCount/nHosts);

			tInfo.HostName[0] += HostInd;
		}
#endif

		//Check for errors

		int ErrorMax;
		MPI_Allreduce(&Error, &ErrorMax, 1, MPI_INT, MPI_MAX, Comm);

		if(ErrorMax > 0){
			if(ThreadId==0)
				printf("Unable to obtain topology\n");
			return false;
		}


		//Thread data

		Threads = new HostThreadInfo[ThreadCount];

		MPI_Allgather(&tInfo, sizeof(HostThreadInfo), MPI_BYTE, Threads, sizeof(HostThreadInfo), MPI_BYTE, Comm);
		
		CountHosts(&HostCount, false);		//Count total number of unique hosts and GPUs

		Hosts = new HostStruct[HostCount];
		
		memset(Hosts, 0, sizeof(HostStruct)*HostCount);

		CountHosts(&HostCount, true);		//Fill new arrays

		FillLists();										//Fills sublists of threads, GPUs and hosts

		return true;
	}

	void PrintTopology(){

		//Write out hosts

		for(int i=0; i<HostCount; i++){

			printf("Host '%s' with %i MPI threads\n", Hosts[i].HostName, Hosts[i].nThreads);
			printf("   Threads: ");

			for(int c=0; c<Hosts[i].nThreads; c++){
				if(c==Hosts[i].nThreads-1) printf("%i", Hosts[i].ThreadIds[c]);
				else					   printf("%i, ", Hosts[i].ThreadIds[c]);
			}

			printf("\n\n");
		}

		/*//Threads

		for(int i=0; i<ThreadCount; i++){

			printf("Thread [%i] on board '%s' (%s GPU)\n", i, Threads[i].HostName, (Threads[i].HasGPU ? "has" : "no"));

		}	*/

	}

	~TopologyMPI(){

		if(Threads!=0) delete[] Threads;
		if(Hosts!=0)   delete[] Hosts;

	}

};



