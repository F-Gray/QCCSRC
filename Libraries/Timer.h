#pragma once

//SimulationTimer class

#include <cmath>
#if _WIN32				//Windows headers
#include <Windows.h>
#else
#include <sys/time.h>
#define _snprintf snprintf
#endif

class SimulationTimer{			//Timer class to record elapsed time

public:
	enum TimeFormat{

		TimeFormat_ms,
		TimeFormat_Sec,
		TimeFormat_Min,
		TimeFormat_Hour,
		TimeFormat_Day,

		TimeFormat_Sec_ms,
		TimeFormat_Min_Sec,
		TimeFormat_Hour_Min_Sec,
		TimeFormat_Day_Hour_Min_Sec,

	};

private:

#if _WIN32	//Windows
	LARGE_INTEGER frequency;        
	LARGE_INTEGER SimulationStart;
	LARGE_INTEGER LastStepTime;
#else
	timeval SimulationStart;
	timeval LastStepTime;
#endif

public:

	SimulationTimer(){

#if _WIN32	//Windows
		QueryPerformanceFrequency(&frequency);
		QueryPerformanceCounter(&SimulationStart);

		LastStepTime = SimulationStart;
#else
		gettimeofday(&SimulationStart, 0);

		LastStepTime.tv_sec = SimulationStart.tv_sec;
		LastStepTime.tv_usec = SimulationStart.tv_usec;
#endif
		
	}

	double GetTimeSinceLastStep(){		//Time since this function last called in seconds

#if _WIN32		//Windows
		LARGE_INTEGER StepTime;
		QueryPerformanceCounter(&StepTime);

		double dt = (double)(StepTime.QuadPart - LastStepTime.QuadPart) / (double)frequency.QuadPart;

		LastStepTime = StepTime;
#else
		timeval StepTime;
		gettimeofday(&StepTime, 0);

		double dt = (double)(StepTime.tv_sec - LastStepTime.tv_sec) + (double)(StepTime.tv_usec - LastStepTime.tv_usec)/1000000.0;

		LastStepTime.tv_sec = StepTime.tv_sec;
		LastStepTime.tv_usec = StepTime.tv_usec;
#endif

		return dt;		//dt in seconds
	}

	double GetTimeSinceLastStepms(){	//Time since this function last called in milliseconds

		return (GetTimeSinceLastStep() * 1000);
	}

	double GetSimulationTime(){			//Time since this SimulationTimer created

#if _WIN32		//Windows
		LARGE_INTEGER Time;
		QueryPerformanceCounter(&Time);

		double dt = (double)(Time.QuadPart - SimulationStart.QuadPart) / (double)frequency.QuadPart;
#else
		timeval Time;
		gettimeofday(&Time, 0);

		double dt = (double)(Time.tv_sec - SimulationStart.tv_sec) + (double)(Time.tv_usec - SimulationStart.tv_usec)/1000000.0;
#endif

		return dt;		//dt in seconds

	}

	void Reset(){						//Set timer start time to now

#if _WIN32	//Windows
		QueryPerformanceFrequency(&frequency);
		QueryPerformanceCounter(&SimulationStart);

		LastStepTime = SimulationStart;
#else
		gettimeofday(&SimulationStart, 0);

		LastStepTime.tv_sec = SimulationStart.tv_sec;
		LastStepTime.tv_usec = SimulationStart.tv_usec;
#endif

	}

	void PrintTime(double Seconds){

		PrintTime(Seconds, TimeFormat_Day_Hour_Min_Sec);

	}

	void PrintTime(double Seconds, TimeFormat Format){

		long long Days  = (long long)floor( Seconds / 86400.0 );
		long long Hours = (long long)floor( Seconds / 3600.0 );
		long long Mins  = (long long)floor( Seconds / 60.0 );
		long long Sec   = (long long)floor( Seconds );
		long long ms    = (long long)floor( Seconds * 1000 );

		if(Format == TimeFormat_Day_Hour_Min_Sec){

			long long Hour_r = (long long)floor( (Seconds - 86400.0*(double)Days) / 3600.0 );
			long long Min_r  = (long long)floor( (Seconds - 86400.0*(double)Days - 3600.0*(double)Hour_r) / 60.0 );
			long long Sec_r  = (long long)floor( (Seconds - 86400.0*(double)Days - 3600.0*(double)Hour_r - 60.0*(double)Min_r) );

			if(Days==0){

				Format = TimeFormat_Hour_Min_Sec;	//Reduce if days = 0

			}else{

				printf("%llid %llih %llim %llis", Days, Hour_r, Min_r, Sec_r);

			}
		}

		if(Format == TimeFormat_Hour_Min_Sec){

			long long Min_r = (long long)floor( (Seconds - 3600.0*(double)Hours) / 60.0 );
			long long Sec_r = (long long)floor( (Seconds - 3600.0*(double)Hours - 60.0*(double)Min_r) );

			if(Hours==0){

				Format = TimeFormat_Min_Sec;	//Reduce if hours = 0

			}else{

				printf("%llih %llim %llis", Hours, Min_r, Sec_r);

			}
		}

		if(Format == TimeFormat_Min_Sec){

			long long Sec_r = (long long)floor(Seconds - 60.0*(double)Mins);

			if(Mins==0){

				Format = TimeFormat_Sec;	//Reduce if minutes = 0

			}else{

				printf("%llim %llis", Mins, Sec_r);

			}
		}

		if(Format == TimeFormat_Sec_ms){
			
			long long ms_r = (long long)floor((Seconds - (double)Sec) * 1000);

			printf("%llis %llims", Sec, ms_r);
		}

		if(Format == TimeFormat_Day)  printf("%llid", Days);
		if(Format == TimeFormat_Hour) printf("%llih", Hours);
		if(Format == TimeFormat_Min)  printf("%llim", Mins);
		if(Format == TimeFormat_Sec)  printf("%llis", Sec);
		if(Format == TimeFormat_ms)   printf("%llims", ms);
		
	}
	
};