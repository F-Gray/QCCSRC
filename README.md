![](ReadMe/QCCSRC.png)

# QCCSRC

Libraries for high performance scientific computing, developed in the Qatar Carbonates and Carbon Storage Research Centre at Imperial College London by Farrel Gray.

![](ReadMe/Graphic.png)

Figure: Fluid flow velocity distribution (middle) and reactant concentration distribution (right) during injection into a rock core sample (left). Computed using a reactive transport model run on GPU cluster. For more details, see [Gray et al. (2018)](https://doi.org/10.1016/j.advwatres.2018.09.007)

## Table of Contents

* [About](#About)
* [Getting Started](#Getting~~Started)
* [Libraries](#Libraries)
* [Contact](#Contact)

## About

The C++ header files included here are used in codes developed within the Qatar Carbonates and Carbon Storage Research Centre (QCCSRC) at Imperial College London. They provide useful library functions for:

* Reading input parameter files
* Buffered file reading/writing with array ordering operations
* 3D data grid class both on local memory systems and distributed over multi-node MPI systems
* Domain decomposition class for multi-node MPI systems
* Multithreading classes
* Accurate timer
* Topology mapping class for multi-node MPI systems, including GPUs

All libraries are compatible with both Linux and Windows systems.

## Getting Started

These static library files can easily be used in a project by including the relevant files, as shown

```C++

#include <InputFile.h>          //Reads input files
#include <DataFiles.h>          //Buffered file reading and writing with array transformation
#include <Grids.h>              //Manipulate 3D grids
#include <Decomposition.h>      //Lattice decomposition
#include <Threads.h>            //Local multi-threading
#include <Timer.h>              //Accurate timing
#include <FilePaths.h>		//Combines folder and file strings

#include <GridsMPI.h>		//Manipulate 3D grids distributed over MPI systems
#include <TopologyMPI.h>        //Map topology of MPI system (no GPUs)

#include <TopologyMPI.cuh>      //Map topology of MPI GPU system

```

The libraries <b>GridsMPI.h</b>, <b>TopologyMPI.h</b> and <b>TopologyMPI.cuh</b> are only usable if MPI libraries are included, and the library <b>TopologyMPI.cuh</b> also requires CUDA to be used. 

## Libraries

### InputFile.h

This class provides a way to read input parameters from a file. It allows a free form structure to be used, and supports both simple data types, arrays and complex nested data-structures.

#### 1. Example with simple data

Here is an example with simple data to be read in from Input.txt.

1.1 <b>Input.txt</b>

  Input values are assigned to named variables using the syntax below, and comments can also be added using C style '//'
    
```Text
//Simple data
    
InputStr = "String for reading"   //C style comments are allowed
		
InputInt = 30                     //Integer to be read in
		
InputFlt = 0.923
    
//Arrays of data
    
ListFlt = 0.035, 0.022, 5.243     //Lists can be separated with commas, spaces, semi-colons...
    
ListFlt2 = { 0.035                //...or spread over multiple lines using curly brackets
             0.022
             5.243 }
```

1.2 <b>InputFile.cpp</b>

In our cpp file we can define the variables to be read in using the <i>InputValueStructEx</i> struct and use the <i>InputFileReader</i> class to perform the read.

```C++

#include <InputFile.h>

//Input file values will be read into these variables

char  InputStr[64];    
int   InputInt;
float InputFlt;

float ListFlt[5];
float ListFlt2[5];

int CountListFlt;      //These will receive a count of how many values were read
int CountListFlt2;

//Define the data struct
    
InputValueStructEx Parameters[] = {

//{ [Name], [Ptr], [Data Length], [Count Max], [Data Type (eg DataType_Int)], [InputValueStructEx], [InputValueStructEx values count], [Ptr Count Read], [Required], [Description] },

  { "InputStr", InputStr,  sizeof(InputStr), 1, DataType_String, 0, 0, 0, false, "An input string" },
  { "InputInt", &InputInt, sizeof(int),      1, DataType_Int,    0, 0, 0, true,  "An input int" },
  { "InputFlt", &InputFlt, sizeof(float),    1, DataType_Float,  0, 0, 0, false, "An input float" },
  
  { "ListFlt",  &ListFlt[0],  sizeof(float), 5, DataType_Float, 0, 0, &CountListFlt,  false, "A list of floats" },
  { "ListFlt2", &ListFlt2[0], sizeof(float), 5, DataType_Float, 0, 0, &CountListFlt2, false, "A list of floats" },
			
};

const int nParams = sizeof(Parameters)/sizeof(InputValueStructEx);

//Perform the file read

void ReadInputFile(){

  InputFileReader InputFile;    //Define instance of InputFileReader class

  char FileName[] = "Input.txt";
  
  bool r = InputFile.ReadInputFile(FileName, Parameters, nParams);  //Read in input file
  
  if(r){    //Success
  
    printf("Values were read successfully\n");   
    
  }else{    //Print error message
  
    printf("Values were not read in: %s\n", InputFile.GetErrorStr(InputFile.GetErrorCode()));   
    
  }
  
}
```

The struct <i>InputValueStructEx</i> is defined

```C++
struct InputValueStructEx{
	const char* DataName;             //Name in input file
	void* Ptr;                        //Pointer to variable
	int   DataSz;                     //Max size of one data entry (bytes)
	int   Count;                      //Max number of values for lists
	int   DataType;                   //Data type enum eg DataType_Bool
	InputValueStructEx* DataStruct;   //If type is DataType_InputValueStructEx, pointer to another InputValueStructEx defining the data within the struct at Ptr (may be null)
	int   DataStructCount;            //If type is DataType_InputValueStructEx, the number of entries in DataStruct (may be null)
	int*  CountRead;                  //Returns the number of values read in (may be null)
	bool  Required;                   //Throws error if value not found in input file
	const char* Description;          //Description of the field (may be null)
};
```

Note that when reading in arrays of strings, the <i>DataSz</i> parameter should be set to the maximum length of each array entry (which includes the ending null character '\0'). Values are read into the char array at intervals of <i>DataSz</i>. For example, if we reading a maximum of 4 strings, each of 31 bytes, we would define our input string variable of 128 bytes, set <i>DataSz</i> = 32, and each array element would written to index 32 * n, for n = {0, 1, 2, 3}.

The DataType enumeration lists supported data types and is defined

```C++
enum InputFileDataType{
	DataType_Bool,
	DataType_Short,
	DataType_Int,
	DataType_LongLongInt,
	DataType_Float,
	DataType_Double,
	DataType_String,
	DataType_InputValueStructEx
};
```

#### 2. Example with complex data

The <i>InputFileReader</i> class also supports reading nested data structures, as shown in the following example

2.1 <b>Input.txt</b>

  Multiple levels of structures can be defined using curly brackets
    
```Text
ComplexData = {   { Name = "Item1"                     //Each element contains a structure of values
                    Value = 80                         //Int
                    FltList = 0.243, 0.334, 0.221 },   //List of floats
                    
                  { Name = "Item2"
                    Value = 27
                    FltList = 0.264, 0.970, 0.396, 0.285, 0.857 },
                    
                  { Name = "Item3"
                    Value = 935
                    FltList = 0.735, 0.497, 0.104 }
              }
```

Datasets such as this can be read in by setting the <i>DataType</i> = DataType_InputValueStructEx in <i>InputValueStructEx</i>, as shown in the following code
  
  2.2 <b>InputFile.cpp</b>

```C++
...

//Input file values will be read into these variables

struct DataStruct{    //This defines the data within each complex data element
  char Name[256];
  int Value;
  float Flt[6];
  int FltCount;
};

DataStruct SubData[4];   //Receives data from complex list
int CountSubData;        //Count how many struct elements were read in

//Define the data structs

InputValueStructEx DataStructParams[] = {     //This InputValueStructEx defines the variables in the sub entries

//{ [Name], [Ptr], [Data Length], [Count Max], [Data Type (eg DataType_Int)], [InputValueStructEx], [InputValueStructEx values count], [Ptr Count Read], [Required], [Description] },

  { "Name",    (void*)offsetof(DataStruct, Name),  sizeof(char)*256, 1, DataType_String, 0, 0, 0, false, 0 },
  { "Value",   (void*)offsetof(DataStruct, Value), sizeof(int),      1, DataType_Int,    0, 0, 0, true, 0 },
  { "FltList", (void*)offsetof(DataStruct, Flt),   sizeof(float),    6, DataType_Float,  0, 0, (int*)offsetof(DataStruct, FltCount), true, 0 },

};
    
InputValueStructEx Parameters[] = {

//{ [Name], [Ptr], [Data Length], [Count Max], [Data Type (eg DataType_Int)], [InputValueStructEx], [InputValueStructEx values count], [Ptr Count Read], [Required], [Description] },

  { "ComplexData", StructRead, sizeof(DataStruct), 4, DataType_InputValueStructEx, DataStructParams, sizeof(DataStructParams)/sizeof(InputValueStructEx), &CountReadStruct, false, 0 },
			
};

...
```

Here we define a second list of parameters using <i>InputValueStructEx</i> struct, which define entries for the sub-data. Now, the "ComplexData" entry in the list of the main parameters has <i>InputValueStructEx.Ptr</i> set to an array of <i>DataStruct</i>, and the <i>InputValueStructEx.DataStruct</i> field is set to the sub parameters list DataStructParams.

Notice how pointers to data locations in the sub data <i>InputValueStructEx</i> are now relative to the struct pointer. This is achieved by using the <i>offsetof</i> macro. Note also that <i>offsetof</i> is only certain to behave correctly if the struct is POD ("plain old data").

#### 3. Example with MPI threads

If MPI is used, and the preprocessor definition MPI_INCLUDED is defined, then a second constructor for the InputFileReader class becomes available

```C++
InputFileReader(MPI_Comm Comm, int SrcThread);
```

All MPI threads initialise the class using the group communicator <i>Comm</i> and the index of the thread which will perform the file read. All threads then call the member function ReadInputFile(...) as in example 1.2, and this now automatically distributes the parameters read from the input file to all threads, and makes the result of GetErrorCode() available to all threads.

### DataFiles.h

### Grids.h

The <i>Grid</i> class provides a container for 3D (as well as 2D and 1D) datasets. It also provides useful capabilities for distributed computing applications, such as the ability to extract subvolumes with extra 'ghost' layers which may be over loop boundaries.

#### 1. Example of 3D grid

The following example creates a 3D grid of type int and reads and writes data from a file.

```C++

int main(){

    Grid<int> DataSet(100, 100, 100);      //Creates a 3D dataset of dimensions 100 x 100 x 100, with a single component

    DataSet.ReadFromFile("C:/Data.raw", true, DataFileBase::Type_char);   //Reads 8-bit data from file in binary format and automatically converts it to int
    
    for(int z=0; z<100; z++)
    for(int y=0; y<100; y++)
    for(int x=0; x<100; x++){
    
	    DataSet(x, y, z) += 2;         //The () operator can be used to return a reference to a particular value
    
    }

    DataSet.WriteToFile("C:/DataOut.raw", true, false);  //Writes data in 32-bit int binary format

    return 0;
}
```
#### 2. Grid Members

A number of different constructors are available
```C++
Grid<class T>(int SizeX,  int SizeY, int SizeZ, int VectorSize);                          //Create grid of dimensions SizeX x SizeY x SizeZ with VectorSize components
Grid<class T>(long long SizeX,  long long SizeY, long long SizeZ, long long VectorSize);
Grid<class T>(int SizeX,  int SizeY, int SizeZ);                                          //Create grid of dimensions SizeX x SizeY x SizeZ with a single component
Grid<class T>(long long SizeX,  long long SizeY, long long SizeZ);
Grid<class T>(int SizeX,  int SizeY);                                                     //Create grid of dimensions SizeX x SizeY with a single component
Grid<class T>(long long SizeX,  long long SizeY);
Grid<class T>(int SizeX);                                                                 //Create 1D array of SizeX elements
Grid<class T>(long long SizeX);

Grid<class T>(Grid<T> &GridRef);                     //Create grid copied from GridRef (reference)
Grid<class T>(Grid<T>* GridRef);                     //Create grid copied from GridRef (pointer)
Grid<class T>(Grid<T> &GridRef, GridRegion Region);  //Create grid copied from a subregion Region of GridRef (reference)
Grid<class T>(Grid<T> &GridRef, Coarsegrain CG_Op);  //Create grid GridRef (reference) coarsegrained

//Create a grid from the array DataArr, of dimensions SizeX x SizeY x SizeZ with VectorSize components, either deep copy (new array allocated) or shallow copy (pointer copied), which will be deleted by the destructor if FreeShallowCopy = true
Grid<class T>(T* DataArr, bool ShallowCopy, bool FreeShallowCopy, long long SizeX,  long long SizeY, long long SizeZ, long long VectorSize);

```

A reference to the data at a particular coordinate can be obtained using the bracket operator (), as in example 1. This has multiple overloads.

```C++

T& Grid<class T>::operator()(long long X, long long Y, long long Z, long long n);  //Returns reference to value at X, Y, Z component n
T& Grid<class T>::operator()(int X, int Y, int Z, int n);
T& Grid<class T>::operator()(long long X, long long Y, long long Z);
T& Grid<class T>::operator()(int X, int Y, int Z);
T& Grid<class T>::operator()(long long X, long long Y);
T& Grid<class T>::operator()(int X, int Y);
T& Grid<class T>::operator()(long long X);
T& Grid<class T>::operator()(int X);

```

Likewise, two methods <i>Get</i> and <i>GetPtr</i> can be used to obtain the value and pointer to the value at a particular coordinate respectively. These have the same overloads as the () operator

```C++
T Grid<class T>::Get(long long X, long long Y, long long Z, long long n);      //Get value at X, Y, Z, component n

T* Grid<class T>::GetPtr(long long X, long long Y, long long Z, long long n);  //Get pointer to element at X, Y, Z, component n
```

A number of functions are also available to set the value at a particular coordinate. These have the same overloads as the Get functions, except that the extra parameter </i>Value</i> is always included as the last parameter

```C++
void Grid<class T>::Set(long long X, long long Y, long long Z, long long n, T Value); //Set value at X, Y, Z, component n
```

Values can also be set throughout the grid, or in specific regions using the <i>SetAll</i> function. Similarly, value replacement can be performed using <i>ReplaceValue</i>, and the function <i>ScaleValues</i> can be used to multiply all values by a given factor

```C++
void Grid<class T>::SetAll(T Value);                          //Set all to Value
void Grid<class T>::SetAll(T Value[]);                        //Set all to Value[], with n components
void Grid<class T>::SetAll(T Value, GridRegion Region);       //Set all to Value in Region
void Grid<class T>::SetAll(T Value[], GridRegion Region);     //Set all to Value[], with n components in Region      

void Grid<class T>::ReplaceValue(T ValueFind, T ValueReplace);                        //Replace all ValueFind with ValueReplace
void Grid<class T>::ReplaceValue(T ValueFind[], T ValueReplace[]);                    //Replace all ValueFind[] with ValueReplace[] (vectors with n components)
void Grid<class T>::ReplaceValue(T ValueFind, T ValueReplace, GridRegion Region);     //Replace all ValueFind with ValueReplace in Region
void Grid<class T>::ReplaceValue(T ValueFind[], T ValueReplace[], GridRegion Region); //Replace all ValueFind[] with ValueReplace[] (vectors with n components) in Region

void Grid<class T>::ScaleValues(T Scale);                         //Multiply all values by Scale
void Grid<class T>::ScaleValues(T Scale, GridRegion Region);      //Multiply all values by Scale in Region
```
#### 3. GridRegion

The <i>GridRegion</i> struct can be used to specify a 3D region and has the following constructors

```C++
GridRegion::GridRegion(long long XStart, long long XEnd, long long YStart, long long YEnd, long long ZStart, long long ZEnd);
GridRegion::GridRegion(int XStart, int XEnd, int YStart, int YEnd, int ZStart, int ZEnd);
GridRegion::GridRegion(long long XStart, long long XEnd, long long YStart, long long YEnd);
GridRegion::GridRegion(int XStart, int XEnd, int YStart, int YEnd);
GridRegion::GridRegion(long long XStart, long long XEnd);
GridRegion::GridRegion(int XStart, int XEnd);
GridRegion::GridRegion();
```
The default range for directions not specified in the constructor parameters is 0, 1.

Individual ranges can be set using

```C++
void GridRegion::SetRangeX(long long XStart, long long XEnd);  //Set range in X
void GridRegion::SetRangeY(long long YStart, long long YEnd);  //Set range in Y
void GridRegion::SetRangeZ(long long ZStart, long long ZEnd);  //Set range in Z
```

This struct can also be used to compute a number of useful metrics

```C++
long long GridRegion::SizeX();         //Size in X
long long GridRegion::SizeY();         //Size in Y
long long GridRegion::SizeZ();         //Size in Z

long long GridRegion::CountCells();    //Number of cells contained in region
bool GridRegion::IsZero();             //Size is 0 in all directions

bool GridRegion::Contains(long long x, long long y, long long z);  //Whether coordinate is within this region
bool GridRegion::TestOverlap(GridRegion& Region);  //Whether Region (ref) overlaps with this region
bool GridRegion::TestOverlap(GridRegion* Region);  //Whether Region (ptr) overlaps with this region
bool GridRegion::FitsWithin(GridRegion& Region);   //Whether Region (ref) fits entirely within this region
bool GridRegion::FitsWithin(GridRegion* Region);   //Whether Region (ptr) fits entirely within this region
```

New regions can be computed from the following operations

```C++
GridRegion GridRegion::Translate(long long dx, long long dy, long long dz);  //Return new GridRegion translated by dx, dy, dz
GridRegion GridRegion::Overlap(GridRegion& Region);         //Return new GridRegion equal to the overlap of Region (ref) with this GridRegion
GridRegion GridRegion::Overlap(GridRegion* Region);         //Return new GridRegion equal to the overlap of Region (ptr) with this GridRegion
```

#### 4. Grid Transformations

A useful capability of the <i>Grid</i> class is the ability to perform coarsegraining operations, extract subgrids and perform inversions in Cartesian directions.

### GridsMPI.h

### Decomposition.h

### Threads.h

This defines a class called <i>Threads</i> which creates local threads and a <i>Mutex</i> class for synchronisation. It provides similar capabilities to the C++11 thread library, however the <i>Threads</i> class has its own inter-thread synchronisation routines built in which simplifies coding for highly parallel applications.

#### 1. Example with thread creation and synchronisation

In this example, we use the <i>Threads</i> class to create a number of parallel threads, and demonstrate the synchronisation routines

```C++

int ThreadWork(Threads::Thread* thread){                             //Runs with multiple threads
        
    int Arg = *((int*)thread->Data);

    printf("Thread Id = %i with argument %i\n", thread->Id, Arg);
        
    thread->SyncThreads();                                           //Synchronise threads

    printf("Finished (%i)\n", thread->Id);

    return 0;
}


int main(){

    Threads Worker;           //Initialise a Threads class

    int Arg = 3;              //Data passed to thread

    Worker.RunThreads(4, ThreadWork, &Arg);   //Runs 'ThreadWork' with 4 threads

    return 0;
}

```
The Threads class is initialised and then the member function RunThreads is called which creates n threads. The ThreadWork function recieves a single parameter of type Threads::Thread which contains the following members

```C++
void* Threads::Thread::Data;                       //The argument passed by RunThreads(...)
int   Threads::Thread::Id;                         //Thread index (from 0 to n-1)
int   Threads::Thread::nThreads;                   //The number of threads running

void  Threads::Thread::SyncThreads();              //Block until all running threads reach this point
void  Threads::Thread::SyncThreads(int* nActive);  //Same as SyncThreads() but returns the current number of running threads in *nActive
```

These members allow each thread to know its index, the number of threads running as well as synchronize with one another. If some threads return before reaching SyncThreads(), then only running threads will wait, and the overload SyncThreads(int* nActive) can be used to determine how many threads are still running.

The Threads::RunThreads(...) function has multiple variants. RunThreads(...) functions block until all threads have returned, whereas the RunThreadsAsync(...) functions return immediately after creating the threads, which now run in parallel to the host thread.

```C++
bool Threads::RunThreads(int n, int (*Addr)(Threads::Thread*) );                 //Run n threads at Addr 
bool Threads::RunThreads(int n, int (*Addr)(Threads::Thread*), void* Data);      //Run n threads at Addr and pass Data as argument (copied to void* Threads::Thread::Data)

bool Threads::RunThreadsAsync(int n, int (*Addr)(Threads::Thread*) );            //Same as RunThreads but returns immediately
bool Threads::RunThreadsAsync(int n, int (*Addr)(Threads::Thread*), void* Data); 
```

The Threads class exposes two further synchronisation routines for the host thread, once RunThreadsAsync has been called. These are

```C++
bool Threads::IsRunning();                  //Check if any threads are still running
bool Threads::ThreadIsRunning(int tId);     //Check if thread with index tId is running

void Threads::WaitFinish();                 //Wait for all threads to complete
```

Finally, the host thread can obtain the int return value of any of the threads using

```C++
int Threads::ThreadReturnValue(int tId);    //Returns the return value of thread with index tId
```

#### 2. Mutex class

The <i>Threads.h</i> library also contains a simple mutex class with a lock function which can be used to control access to resources in a system with multiple threads. Its members are

```C++
Mutex::Mutex();                //Default constructor

void Mutex::Lock(bool Lock);   //Acquire the lock (Lock = true) or release the lock (Lock = false)
```

### Timer.h

This file provides a simple high resolution timer class called <i>SimulationTimer</i>, which compatible with both Linux and Windows.

#### 1. Example

```C++

int main(){

    SimulationTimer Timer;      //Initialise timer class, which automatically begins counting
    
    while(true){
    
    	sleep(1);
    
        double dtStep = Timer.GetTimeSinceLastStep();  //Returns time since GetTimeSinceLastStep() was last called in seconds
        double dtTot  = Timer.GetSimulationTime();     //Returns time since SimulationTimer class was intialised in seconds
        
	printf("Time since last check: %.6f s, total time: %.6f s\n", dtStep, dtTot);
    }

    return 0;
}

```

The following methods can be used to obtain time intervals


```C++

Timer::GetTimeSinceLastStep();      //Returns time since GetTimeSinceLastStep or GetTimeSinceLastStepms was last called in seconds
Timer::GetTimeSinceLastStepms();    //Returns time since GetTimeSinceLastStep or GetTimeSinceLastStepms was last called in milliseconds
Timer::GetSimulationTime();         //Returns time since SimulationTimer class was intialised in seconds

```

The timer class can also be reset using

```C++
Timer::Reset();    //Reset timer
```

A variety of different time formats can be printed to the console using

```C++
Timer::PrintTime(double Seconds);                      //Print time with default format (TimeFormat_Day_Hour_Min_Sec)
Timer::PrintTime(double Seconds, TimeFormat Format);   //Print time with chosen format from TimeFormat enumeration
```

Time formats can be chosen from the <i>TimeFormat</i> enumeration


```C++
enum TimeFormat{
    TimeFormat_ms,                //Milliseconds only e.g. "10ms"
    TimeFormat_Sec,               //Seconds only e.g. "10s"
    TimeFormat_Min,               //Minutes only e.g. "10m"
    TimeFormat_Hour,              //Hours only e.g. "10h"
    TimeFormat_Day,               //Days only e.g. "10d"

    TimeFormat_Sec_ms,            //Seconds and milliseconds e.g. "10s 305ms"
    TimeFormat_Min_Sec,           //Minutes and seconds e.g. "10m 30s"
    TimeFormat_Hour_Min_Sec,      //Hours, minutes and seconds e.g. "5h 10m 3s"
    TimeFormat_Day_Hour_Min_Sec,  //Days, hours, minutes and seconds e.g. "4d 7h 37m 9s"
};

```

### FilePaths.h

### Topology.cuh
