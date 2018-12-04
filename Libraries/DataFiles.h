#pragma once

#include <fstream>
#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#define _snprintf snprintf
#define _fseeki64 fseeko64
#define _ftelli64 ftello64
#endif

#include "Threads.h"			//Threading class
#include "DataFilesBase.h"		//Type conversions and array index counters and transformations
#include "Base.h"

#pragma warning( disable : 4996 )	//Potentially unsafe _snprintf and fopen

////////////////////////////////////////////////////////////////////////////////////////////

struct ArrayReadOp;			//Read operation for DataFileReader

struct ArrayWriteOp;		//Read operation for DataFileWriter

class DataFileReader;		//Reads arrays from files

class DataFileWriter;		//Writers arrays to files






////////////////////////////////////////////////////////////////////////////////////////////

typedef DataFileBase::DataType DataType;
typedef DataFileBase::Rotation Rotation;


	//Read operation for DataFileReader

struct ArrayReadOp{			//Defines an array read operation

	long long arrSize[3];			//File array dimensions (up to 3D)
	int nVector;					//Number of values per array element

	void* DataArray;				//Data array
	long long DataArrSize[3];		//Dimensions of output data array
	long long DataArrOffset[3];		//Offset writing into data array

	bool Bin;						//File data in bytes, else ASCII string
	DataType FileDataType;			//Data type in file
	DataType ArrayDataType;			//Data type of array elements

	Rotation Transform;				//Rotate grid: 0 = no rotation; 1 = XY -> YX or XYZ -> YZX; 2 = XYZ -> ZXY

private:
	void Set(long long _arrSize[3], int _nVector, void* _DataArray, long long _DataArrSize[3], long long _DataArrOffset[3], bool _Bin, DataType _FileType, DataType _ArrType, Rotation _Transform){

		arrSize[0] = _arrSize[0];
		arrSize[1] = _arrSize[1];
		arrSize[2] = _arrSize[2];

		nVector = _nVector;

		DataArray = _DataArray;

		DataArrSize[0] = _DataArrSize[0];
		DataArrSize[1] = _DataArrSize[1];
		DataArrSize[2] = _DataArrSize[2];

		DataArrOffset[0] = _DataArrOffset[0];
		DataArrOffset[1] = _DataArrOffset[1];
		DataArrOffset[2] = _DataArrOffset[2];

		Bin = _Bin;

		FileDataType = _FileType;
		ArrayDataType = _ArrType;

		Transform = _Transform;
	}

public:
	
	ArrayReadOp(){

		long long arrSz[3] = {0, 1, 1};
		long long Offset[3] = {0, 0, 0};

		Set(arrSz, 1, 0, arrSz, Offset, false, DataFileBase::Type_int32, DataFileBase::Type_int32, DataFileBase::RotateGrid_XYZ);
	}

	ArrayReadOp(long long Count, void* DataPtr, DataType arrType){

		long long arrSz[3] = {0, 1, 1};
		arrSz[0] = Count;

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, 1, DataPtr, arrSz, Offset, false, arrType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayReadOp(bool Binary, long long Count, void* DataPtr, DataType arrType){

		long long arrSz[3] = {0, 1, 1};
		arrSz[0] = Count;

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, 1, DataPtr, arrSz, Offset, Binary, arrType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayReadOp(bool Binary, long long Count, DataType FileType, void* DataPtr, DataType arrType){

		long long arrSz[3] = {0, 1, 1};
		arrSz[0] = Count;

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, 1, DataPtr, arrSz, Offset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayReadOp(bool Binary, long long Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType){

		long long Offset[3] = {0, 0, 0};

		Set(Count, VectorSz, DataPtr, Count, Offset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayReadOp(bool Binary, long long Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, Rotation _Transform){

		long long Offset[3] = {0, 0, 0};

		Set(Count, VectorSz, DataPtr, Count, Offset, Binary, FileType, arrType, _Transform);
	}

	ArrayReadOp(bool Binary, long long Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, long long DataArrSize[3], long long DataArrOffset[3]){

		Set(Count, VectorSz, DataPtr, DataArrSize, DataArrOffset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayReadOp(bool Binary, long long Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, long long DataArrSize[3], long long DataArrOffset[3], Rotation _Transform){

		Set(Count, VectorSz, DataPtr, DataArrSize, DataArrOffset, Binary, FileType, arrType, _Transform);
	}

	ArrayReadOp(bool Binary, int Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType){

		long long arrSz[3];

		arrSz[0] = Count[0];
		arrSz[1] = Count[1];
		arrSz[2] = Count[2];

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, VectorSz, DataPtr, arrSz, Offset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayReadOp(bool Binary, int Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, Rotation _Transform){

		long long arrSz[3];

		arrSz[0] = Count[0];
		arrSz[1] = Count[1];
		arrSz[2] = Count[2];

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, VectorSz, DataPtr, arrSz, Offset, Binary, FileType, arrType, _Transform);
	}

	ArrayReadOp(bool Binary, int Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, int DataArrSize[3], int DataArrOffset[3]){

		long long _FileArrSz[3];
		long long _DataArrSz[3];
		long long _DataArrOffset[3];

		for(int i=0; i<3; i++){
			_FileArrSz[i] = Count[i];
			_DataArrSz[i] = DataArrSize[i];
			_DataArrOffset[i] = DataArrOffset[i];
		}

		Set(_FileArrSz, VectorSz, DataPtr, _DataArrSz, _DataArrOffset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayReadOp(bool Binary, int Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, int DataArrSize[3], int DataArrOffset[3], Rotation _Transform){

		long long _FileArrSz[3];
		long long _DataArrSz[3];
		long long _DataArrOffset[3];

		for(int i=0; i<3; i++){
			_FileArrSz[i] = Count[i];
			_DataArrSz[i] = DataArrSize[i];
			_DataArrOffset[i] = DataArrOffset[i];
		}

		Set(_FileArrSz, VectorSz, DataPtr, _DataArrSz, _DataArrOffset, Binary, FileType, arrType, _Transform);
	}
};

struct ArrayWriteOp{		//Defines an array write operation

	//Data

	long long arrSize[3];			//File array dimensions (up to 3D)
	int nVector;					//Number of values per array element

	void* DataArray;				//Data array
	long long DataArrSize[3];		//Dimensions of output data array
	long long DataArrOffset[3];		//Offset writing into data array

	bool Bin;						//File data in bytes, else ASCII string
	DataType FileDataType;			//Data type in file
	DataType ArrayDataType;			//Data type of array elements

	Rotation Transform;				//Rotate grid: 0 = no rotation; 1 = XY -> YX or XYZ -> YZX; 2 = XYZ -> ZXY

	//Format

	char WriteHeaderName[256];		//VTK Header
	char WriteDatasetName[256];
	bool HasHeader;

	char WriteSpacer;				//Spacer between output values (ASCII)
	int nDecimalPlacesFlt;			//Number of decimal places to write float/double values

private:
	void Set(long long _arrSize[3], int _nVector, void* _DataArray, long long _DataArrSize[3], long long _DataArrOffset[3], bool _Bin, DataType _FileType, DataType _ArrType, Rotation _Transform){

		arrSize[0] = _arrSize[0];
		arrSize[1] = _arrSize[1];
		arrSize[2] = _arrSize[2];

		nVector = _nVector;

		DataArray = _DataArray;

		DataArrSize[0] = _DataArrSize[0];
		DataArrSize[1] = _DataArrSize[1];
		DataArrSize[2] = _DataArrSize[2];

		DataArrOffset[0] = _DataArrOffset[0];
		DataArrOffset[1] = _DataArrOffset[1];
		DataArrOffset[2] = _DataArrOffset[2];

		Bin = _Bin;

		FileDataType = _FileType;
		ArrayDataType = _ArrType;

		Transform = _Transform;

		SetFormatDefault();
	}

	void SetFormatDefault(){

		HasHeader = false;
		WriteHeaderName[0] = '\0';
		WriteDatasetName[0] = '\0';
		WriteSpacer = ' ';
		nDecimalPlacesFlt = 3;
	}

public:

	//Data

	ArrayWriteOp(){

		long long arrSz[3] = {0, 1, 1};
		long long Offset[3] = {0, 0, 0};

		Set(arrSz, 1, 0, arrSz, Offset, false, DataFileBase::Type_int32, DataFileBase::Type_int32, DataFileBase::RotateGrid_XYZ);
	}

	ArrayWriteOp(long long Count, void* DataPtr, DataType arrType){

		long long arrSz[3] = {0, 1, 1};
		arrSz[0] = Count;

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, 1, DataPtr, arrSz, Offset, false, arrType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayWriteOp(bool Binary, long long Count, void* DataPtr, DataType arrType){

		long long arrSz[3] = {0, 1, 1};
		arrSz[0] = Count;

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, 1, DataPtr, arrSz, Offset, Binary, arrType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayWriteOp(bool Binary, long long Count, DataType FileType, void* DataPtr, DataType arrType){

		long long arrSz[3] = {0, 1, 1};
		arrSz[0] = Count;

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, 1, DataPtr, arrSz, Offset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayWriteOp(bool Binary, long long Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType){

		long long Offset[3] = {0, 0, 0};

		Set(Count, VectorSz, DataPtr, Count, Offset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayWriteOp(bool Binary, long long Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, Rotation _Transform){

		long long Offset[3] = {0, 0, 0};

		Set(Count, VectorSz, DataPtr, Count, Offset, Binary, FileType, arrType, _Transform);
	}

	ArrayWriteOp(bool Binary, long long Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, long long DataArrSize[3], long long DataArrOffset[3]){

		Set(Count, VectorSz, DataPtr, DataArrSize, DataArrOffset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayWriteOp(bool Binary, long long Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, long long DataArrSize[3], long long DataArrOffset[3], Rotation _Transform){

		Set(Count, VectorSz, DataPtr, DataArrSize, DataArrOffset, Binary, FileType, arrType, _Transform);
	}

	ArrayWriteOp(bool Binary, int Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType){

		long long arrSz[3];

		arrSz[0] = Count[0];
		arrSz[1] = Count[1];
		arrSz[2] = Count[2];

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, VectorSz, DataPtr, arrSz, Offset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayWriteOp(bool Binary, int Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, Rotation _Transform){

		long long arrSz[3];

		arrSz[0] = Count[0];
		arrSz[1] = Count[1];
		arrSz[2] = Count[2];

		long long Offset[3] = {0, 0, 0};

		Set(arrSz, VectorSz, DataPtr, arrSz, Offset, Binary, FileType, arrType, _Transform);
	}

	ArrayWriteOp(bool Binary, int Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, int DataArrSize[3], int DataArrOffset[3]){

		long long _FileArrSz[3];
		long long _DataArrSz[3];
		long long _DataArrOffset[3];

		for(int i=0; i<3; i++){
			_FileArrSz[i] = Count[i];
			_DataArrSz[i] = DataArrSize[i];
			_DataArrOffset[i] = DataArrOffset[i];
		}

		Set(_FileArrSz, VectorSz, DataPtr, _DataArrSz, _DataArrOffset, Binary, FileType, arrType, DataFileBase::RotateGrid_XYZ);
	}

	ArrayWriteOp(bool Binary, int Count[3], int VectorSz, DataType FileType, void* DataPtr, DataType arrType, int DataArrSize[3], int DataArrOffset[3], Rotation _Transform){

		long long _FileArrSz[3];
		long long _DataArrSz[3];
		long long _DataArrOffset[3];

		for(int i=0; i<3; i++){
			_FileArrSz[i] = Count[i];
			_DataArrSz[i] = DataArrSize[i];
			_DataArrOffset[i] = DataArrOffset[i];
		}

		Set(_FileArrSz, VectorSz, DataPtr, _DataArrSz, _DataArrOffset, Binary, FileType, arrType, _Transform);
	}

	//Format

	void SetSpacer(char Spacer){
		WriteSpacer = Spacer;
	}

	void SetFltDecimalPlaces(int nDP){
		nDecimalPlacesFlt = nDP;
	}

	void AddVTKHeader(const char* HeaderName, const char* DatasetName){
		int l0 = min( (int)strlen(HeaderName), (int)sizeof(WriteHeaderName)-1 );
		int l1 = min( (int)strlen(DatasetName), (int)sizeof(WriteDatasetName)-1 );

		memcpy(WriteHeaderName, HeaderName, l0);
		memcpy(WriteDatasetName, DatasetName, l1);

		WriteHeaderName[l0] = '\0';
		WriteDatasetName[l1] = '\0';

		HasHeader = true;
	}

	void RemoveHeader(){
		HasHeader = false;
	}

};


	//Read arrays from files

class DataFileReader{

private:

	//File

	FILE* File;					//File, opened for read. = 0 if closed
	long long FileLength;		//File length

	char* buf;					//Read/write buffer
	long long bufSize;			//Read/write buffer length
	long long bufLength;		//Number of chars read from file into buffer
	long long bufPos;			//Read position in buffer

	//Status

	int ErrorCode;

	//Operation

	long long arrSize[3];		//File array dimensions (up to 3D)
	int nVector;				//Number of values per array element
	void* DataArray;			//Data array
	long long DataArrSize[3];	//Dimensions of output data array
	long long DataArrOffset[3];	//Offset writing into data array
	Rotation Transform;			//Rotate grid: 0 = no rotation; 1 = XY -> YX or XYZ -> YZX; 2 = XYZ -> ZXY
	bool Bin;					//File data in bytes, else ASCII string
	DataType FileDataType;		//Data type in file
	DataType ArrayDataType;		//Data type of array elements


		//Initialise and finalise class data

	void Initialise(long long bufSizeDefault){

		//File

		File = 0;
		buf = 0;
		bufLength = 0;
		bufPos = 0;

		bufSize = bufSizeDefault;
		buf = new char[bufSize];

		//Status

		ErrorCode = 0;

	}

	void Finalise(){

		delete[] buf;

	}


		//Set status

	void SetErrorCode(int eCode){

		ErrorCode = eCode;

	}


		//File operations

	long long GetFileLength(bool SetEnd){		//Get file length

		long long fPos = _ftelli64(File);		//Current position

		_fseeki64(File, 0, SEEK_END);			//Seek end of file

		long long flength = _ftelli64(File);	//Stream position of end of file

		if(!SetEnd)
			_fseeki64(File, fPos, SEEK_SET);	//Return to original stream position

		return flength;
	}

	void SeekRead(long long Pos, long long nBytes, char* Dest){		//Seek file position and read

		_fseeki64(File, Pos, SEEK_SET);

		fread(Dest, 1, nBytes, File);
	}


		//Sets an operation to perform

	bool SetReadOp(ArrayReadOp* Op){

		//Set class variables

		for(int i=0; i<3; i++){
			arrSize[i] = Op->arrSize[i];
			DataArrSize[i] = Op->DataArrSize[i];
			DataArrOffset[i] = Op->DataArrOffset[i];
		}

		nVector = Op->nVector;
		DataArray = Op->DataArray;
		Transform = Op->Transform;
		Bin = Op->Bin;
		FileDataType = Op->FileDataType;
		ArrayDataType = Op->ArrayDataType;

		//Validate

		if(File == 0){
			SetErrorCode(1);		//File open failed
			return false;
		}

		if(buf==0){
			SetErrorCode(4);		//Buffer allocation failed
			return false;
		}
	
		return true;
	}


		//Array read functions

	bool _ReadArray(){				//Reads in array from current file position

		long long nValues = arrSize[0]*arrSize[1]*arrSize[2];

		if(nValues == 0){
			SetErrorCode(0);		//Array size not properly set
			return true;
		}

		//File stream position

		long long flength = GetFileLength(false);						//Total file length
		long long fReadPos = _ftelli64(File) - bufLength + bufPos;		//Current buffered read position

		long long nRead = _ftelli64(File);		//Total data read from file into buffer

		const int vBufSz = 256;
		char vBuf[vBufSz];				//Current value buffer
		int vBufL = 0;
		vBuf[vBufSz-1] = '\0';

		int binTypeSz = DataTypeSize(FileDataType);		//For bin files, number of bytes per value

		long long arrOffset = DataTypeSize(ArrayDataType);

		long long nValuesTot = nValues * nVector;		//Number of values x number of elements in each

		//Array sizes and counters

		long long arrCount[4] = {0, 0, 0, 0};			//4D data counter

		long long arrSize4[4];
		arrSize4[0] = nVector;
		arrSize4[1] = arrSize[0];
		arrSize4[2] = arrSize[1];
		arrSize4[3] = arrSize[2];

		long long DataArrSize4[4];
		DataArrSize4[0] = nVector;
		DataArrSize4[1] = DataArrSize[0];
		DataArrSize4[2] = DataArrSize[1];
		DataArrSize4[3] = DataArrSize[2];

		//Read

		for(long long nCount=0; nCount<nValuesTot;   ){

			if(bufPos == bufLength){		//Read data into buffer

				bufLength = min( flength-nRead , bufSize);

				if(bufLength==0){
					SetErrorCode(2);
					return false;
				}

				fread(buf, 1, bufLength, File);

				nRead += bufLength;
				bufPos = 0;
			}

			//Obtain value

			char c = buf[bufPos++];

				//Binary

			if(Bin){	

				vBuf[vBufL++] = c;

				if(vBufL == binTypeSz){

					long long n = DataFileBase::arrIndexFwd(arrCount, arrSize4, DataArrSize4, DataArrOffset, Transform);

					void* arrPtr = (void*)&(((char*)DataArray)[n*arrOffset]);	//Array element to write to

					DataFileBase::StrToVal(vBuf, true, arrPtr, FileDataType, ArrayDataType);

					vBufL = 0;
					nCount++;
				}

				continue;
			}

				//ASCII

			bool vEndChar = (c==' '||c=='\r'||c=='\n'||c=='\t'||c==','||c==';'||c=='\0');	//End of value character
			bool LastChar = (bufPos==bufLength && nRead==flength);							//Last character in file

			if(LastChar && !vEndChar){		//Last character in file is not a value-terminating character
		
				if(vBufL == vBufSz-1){		//Buffer overflow
					SetErrorCode(3);
					return false;
				}

				vBuf[vBufL++] = c;			//Read in last character
			}

			if(vEndChar||LastChar){			//End of value

				if(vBufL==0)
					continue;

				vBuf[vBufL] = '\0';

				long long n = DataFileBase::arrIndexFwd(arrCount, arrSize4, DataArrSize4, DataArrOffset, Transform);

				void* arrPtr = (void*)&(((char*)DataArray)[n*arrOffset]);	//Array element to write to

				StrToVal(vBuf, false, arrPtr, FileDataType, ArrayDataType);

				vBufL = 0;
				nCount++;

			}else{

				if(vBufL == vBufSz-1){	//Buffer overflow
					SetErrorCode(3);
					return false;
				}

				vBuf[vBufL++] = c;		//Read next character
			}

		}

		return true;
	}

	bool _ReadArrayFromEnd(){		//Reads in array from current file position, in reverse

		long long nValues = arrSize[0]*arrSize[1]*arrSize[2];

		if(nValues == 0){
			SetErrorCode(0);		//Array size not properly set
			return true;
		}

		//File stream position

		long long flength = GetFileLength(false);			//Current file stream position
		long long fReadPos = flength - bufLength + bufPos;	//Current read position

		if(fReadPos==0){						//At beginning of file
			flength = GetFileLength(true);		//Set to end
			bufLength = 0;						//Reset buffer
			bufPos    = 0;
		}

		long long nRead = _ftelli64(File) + bufLength;		//Total data read from file into buffer

		const int vBufSz = 256;
		char vBuf[vBufSz];				//Current value buffer
		int vBufL = 0;
		vBuf[vBufSz-1] = '\0';

		int binTypeSz = DataTypeSize(FileDataType);		//For bin files, number of bytes per value

		long long arrOffset = DataTypeSize(ArrayDataType);

		long long nValuesTot = nValues * nVector;	//Number of values x number of elements in each

		long long arrCount[4];			//4D data counter
		arrCount[0] = nVector-1;
		arrCount[1] = arrSize[0]-1;
		arrCount[2] = arrSize[1]-1;
		arrCount[3] = arrSize[2]-1;

		long long arrSize4[4];			//4D array size
		arrSize4[0] = nVector;
		arrSize4[1] = arrSize[0];
		arrSize4[2] = arrSize[1];
		arrSize4[3] = arrSize[2];

		long long DataArrSize4[4];
		DataArrSize4[0] = nVector;
		DataArrSize4[1] = DataArrSize[0];
		DataArrSize4[2] = DataArrSize[1];
		DataArrSize4[3] = DataArrSize[2];

		//Read

		for(long long nCount=0; nCount<nValuesTot;   ){

			if(bufPos == 0){		//Read data into buffer

				bufLength = min( flength-nRead , bufSize);

				if(bufLength==0){
					SetErrorCode(2);
					return false;
				}

				SeekRead(flength-nRead-bufLength, bufLength, buf);

				nRead += bufLength;
				bufPos = bufLength;
			}

			//Obtain value

			char c = buf[(--bufPos)];

				//Binary

			if(Bin){	

				vBuf[binTypeSz - (vBufL++) - 1] = c;

				if(vBufL == binTypeSz){

					long long n = DataFileBase::arrIndexBack(arrCount, arrSize4, DataArrSize4, DataArrOffset, Transform);

					void* arrPtr = (void*)&(((char*)DataArray)[n*arrOffset]);	//Array element to write to

					StrToVal(vBuf, true, arrPtr, FileDataType, ArrayDataType);

					vBufL = 0;
					nCount++;
				}

				continue;
			}

				//ASCII

			bool vEndChar  = (c==' '||c=='\r'||c=='\n'||c=='\t'||c==','||c==';'||c=='\0');	//End of value character
			bool FirstChar = (bufPos==0 && nRead==flength);									//First character in file

			if(FirstChar && !vEndChar){		//First character in file is not a value-terminating character
		
				if(vBufL == vBufSz){	//Buffer overflow
					SetErrorCode(3);
					return false;
				}

				vBuf[vBufSz - (vBufL++) - 2] = c;	//Read in first character
			}

			if(vEndChar||FirstChar){		//End of value

				if(vBufL==0)
					continue;

				char* vStr = &vBuf[vBufSz-vBufL-1];							//Data element from file

				long long n = DataFileBase::arrIndexBack(arrCount, arrSize4, DataArrSize4, DataArrOffset, Transform);

				void* arrPtr = (void*)&(((char*)DataArray)[n*arrOffset]);	//Array element to write to

				StrToVal(vStr, false, arrPtr, FileDataType, ArrayDataType);

				vBufL = 0;
				nCount++;

			}else{

				if(vBufL == vBufSz){	//Buffer overflow
					SetErrorCode(3);
					return false;
				}

				vBuf[vBufSz - (vBufL++) - 2] = c;	//Read next character
			}

		}

		return true;
	}


public:

		//Initialise class

	DataFileReader(){

		long long BufSzDefault = 1048576;	//1MB read buffer

		Initialise(BufSzDefault);

	}

	DataFileReader(long long BufferSize){

		Initialise(BufferSize);

	}


		//Open and close file

	bool Open(char* FileName){

		return Open(FileName, false);
	}

	bool Open(char* FileName, bool SeekEnd){

		Close();						//Close current file

		File = fopen(FileName, "rb");

		if(File == 0){
			SetErrorCode(1);			//File open failed
			return false;
		}

		FileLength = GetFileLength(SeekEnd);			//Current file stream position

		SetErrorCode(0);

		return true;
	}

	void Close(){						//Close file

		if(File==0)
			return;

		fclose(File);
		File = 0;
	}

	
		//Read methods

	bool ReadArray(ArrayReadOp* Operation){		//Read array from file

		if(!SetReadOp(Operation))
			return false;

		return _ReadArray();

	}

	bool ReadArrayFromEnd(ArrayReadOp* Operation){		//Read array from end of file

		if(!SetReadOp(Operation))
			return false;

		return _ReadArrayFromEnd();

	}

	bool ReadDataBack(void* DataArray, long long Count, bool Binary, DataType FileType, DataType ArrayType){

		long long fPos = _ftelli64(File);
	
		int binTypeSz = DataTypeSize(FileType);		//For bin files, number of bytes per value

		long long arrOffset = DataTypeSize(ArrayType);

		const int vBufSz = 256;
		char vBuf[vBufSz];				//Current value buffer
		int vBufL = 0;
		vBuf[vBufSz-1] = '\0';

		for(long long n=0; n<Count;   ){

			if(bufPos == 0){		//Read data into buffer

				long long bufOld = bufLength;

				bufLength = min( fPos-bufOld , bufSize);		//Next buffer size

				if(bufLength==0){
					SetErrorCode(2);
					return false;
				}

				SeekRead( fPos-bufOld-bufLength, bufLength, buf);

				fPos -= bufOld;
				bufPos = bufLength;
			}

			//Obtain value

			char c = buf[(--bufPos)];

				//Binary

			if(Binary){	

				vBuf[binTypeSz - (vBufL++) - 1] = c;

				if(vBufL == binTypeSz){

					void* arrPtr = (void*)&(((char*)DataArray)[(Count-n-1)*arrOffset]);		//Array element to write to

					if(FileType != ArrayType){
						StrToVal(vBuf, true, arrPtr, FileType, ArrayType);		//Type conversion to array
					}else{
						for(int i=0; i<binTypeSz; i++)
							((char*)arrPtr)[i] = vBuf[i];		//Copy to array
					}

					vBufL = 0;
					n++;
				}

				continue;
			}

				//ASCII

			bool vEndChar  = (c==' '||c=='\r'||c=='\n'||c=='\t'||c==','||c==';'||c=='\0');		//End of value character
			bool FirstChar = (bufPos==0 && fPos==bufLength);									//First character in file

			if(FirstChar && !vEndChar){		//First character in file is not a value-terminating character
		
				if(vBufL == vBufSz){	//Buffer overflow
					SetErrorCode(3);
					return false;
				}

				vBuf[vBufSz - (vBufL++) - 2] = c;	//Read in first character
			}

			if(vEndChar||FirstChar){		//End of value

				if(vBufL==0)
					continue;

				char* vStr = &vBuf[vBufSz-vBufL-1];							//Data element from file

				void* arrPtr = (void*)&(((char*)DataArray)[(Count-n-1)*arrOffset]);	//Array element to write to

				StrToVal(vStr, false, arrPtr, FileType, ArrayType);

				vBufL = 0;
				n++;

			}else{

				if(vBufL == vBufSz){	//Buffer overflow
					SetErrorCode(3);
					return false;
				}

				vBuf[vBufSz - (vBufL++) - 2] = c;	//Read next character
			}

		}

		return true;
	}

	
		//Obtain status

	const char* ErrorString(){

		int eCode = LastError();

		const char* Strings[] = {
			"The file was successfully read",
			"The file could not be opened",
			"The file contained insufficient data",
			"The file contained invalid data",
			"Memory allocation error",
		};

		return Strings[eCode];
	}

	int LastError(){

		return ErrorCode;
	}


		//Finalise

	~DataFileReader(){

		Close();

		Finalise();

	}

};


//Write arrays to files

class DataFileWriter{

private:

	//File

	FILE* File;					//File, opened for write. = 0 if closed
	char* buf;					//Write buffer
	long long bufSize;			//Write buffer length
	long long bufPos;			//Read position in buffer

	//Status

	int ErrorCode;

	//Operation

	long long arrSize[3];			//File array dimensions (up to 3D)
	int nVector;					//Number of values per array element

	void* DataArray;				//Data array
	long long DataArrSize[3];		//Dimensions of output data array
	long long DataArrOffset[3];		//Offset writing into data array

	bool Bin;						//File data in bytes, else ASCII string
	DataType FileDataType;			//Data type in file
	DataType ArrayDataType;			//Data type of array elements

	Rotation Transform;				//Rotate grid: 0 = no rotation; 1 = XY -> YX or XYZ -> YZX; 2 = XYZ -> ZXY

	int nDecimalPlaces;				//Number of decimal places to write floating types
	char WriteSpacer;				//Spacer between output values (ASCII)

	char HeaderName[256];			//VTK Header
	char DatasetName[256];
	bool HasHeader;


		//Initialise and finalise class data

	void Initialise(long long bufSizeDefault){

		//File

		File = 0;
		buf = 0;
		bufPos = 0;

		bufSize = bufSizeDefault;
		buf = new char[bufSize];

		//Status

		ErrorCode = 0;

	}

	void Finalise(){

		delete[] buf;

	}


		//Set status

	void SetErrorCode(int eCode){

		ErrorCode = eCode;

	}


		//File operations

	void FlushBuffer(){					//Write the send buffer to file

		if(bufPos==0)
			return;

		fwrite(buf, 1, bufPos, File);

		bufPos = 0;

	}

	bool BufferedWrite(const char* Data, long long Length){

		bool Flushed = false;

		if(Length <= (bufSize-bufPos)){					//Data fits in buffer capacity

			memcpy(&buf[bufPos], Data, Length);

			bufPos+=Length;

			if(bufPos == bufSize){
				FlushBuffer();
				Flushed = true;
			}

		}else if(Length >= bufSize){					//Data larger than buffer

			FlushBuffer();
			Flushed = true;

			fwrite(Data, 1, Length, File);

		}else{											//Buffered write

			for(long long c=0; c<Length; ){
				
				long long n = min( bufSize-bufPos, Length-c );

				memcpy(&buf[bufPos], &Data[c], n);

				bufPos+=n;
				c+=n;

				if(bufPos == bufSize){
					FlushBuffer();
					Flushed = true;
				}
			}

		}

		return Flushed;
	}


		//Sets an operation to perform

	bool SetWriteOp(ArrayWriteOp* Op){

		//Set class variables

		for(int i=0; i<3; i++){
			arrSize[i] = Op->arrSize[i];
			DataArrSize[i] = Op->DataArrSize[i];
			DataArrOffset[i] = Op->DataArrOffset[i];
		}

		nVector = Op->nVector;
		DataArray = Op->DataArray;
		Transform = Op->Transform;
		Bin = Op->Bin;
		FileDataType = Op->FileDataType;
		ArrayDataType = Op->ArrayDataType;

		WriteSpacer = Op->WriteSpacer;
		nDecimalPlaces = Op->nDecimalPlacesFlt;

		HasHeader = Op->HasHeader;

		if(HasHeader){

		int l0 = min( (int)strlen(Op->WriteHeaderName)+1, (int)sizeof(HeaderName) );
		int l1 = min( (int)strlen(Op->WriteDatasetName)+1, (int)sizeof(DatasetName) );

		memcpy(HeaderName, Op->WriteHeaderName, l0);
		memcpy(DatasetName, Op->WriteDatasetName, l1);

		}

		//Validate

		if(File == 0){
			SetErrorCode(1);		//File open failed
			return false;
		}

		if(buf==0){
			SetErrorCode(2);		//Buffer allocation failed
			return false;
		}
	
		return true;
	}


		//Array write functions

	void _WriteVTKHeader(char* _HeaderName, char* _DataSetName, long long _GridSz[3], bool _Vectors){
	
		_WriteVTKHeader( _HeaderName, _DataSetName, _GridSz, false, _Vectors);
	}

	void _WriteVTKHeader(char* _HeaderName, char* _DataSetName, long long _GridSz[3], bool _Binary, bool _Vectors){

		long long nValues = _GridSz[0]*_GridSz[1]*_GridSz[2];

		char HeaderBuf[1024];

		int l = _snprintf(HeaderBuf, 1024,

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

				_HeaderName,
				_Binary ? "binary" : "ASCII",
				_GridSz[0], _GridSz[1], _GridSz[2],
				nValues,
				_Vectors ? "Vectors" : "Scalars",
				_DataSetName,
				!_Vectors ? "LOOKUP_TABLE default\n" : "" );

		BufferedWrite(HeaderBuf, l);

	}

	void _WriteVTKHeader(){

		long long nValues = arrSize[0]*arrSize[1]*arrSize[2];

		int arrTr[3];
		arrTr[0] = (Transform == DataFileBase::RotateGrid_XYZ) ? 0 : ( (Transform == DataFileBase::RotateGrid_YZX) ? 1 : 2 );
		arrTr[1] = (Transform == DataFileBase::RotateGrid_XYZ) ? 1 : ( (Transform == DataFileBase::RotateGrid_YZX) ? 2 : 0 );
		arrTr[2] = (Transform == DataFileBase::RotateGrid_XYZ) ? 2 : ( (Transform == DataFileBase::RotateGrid_YZX) ? 0 : 1 );

		char HeaderBuf[1024];

		int l = _snprintf(HeaderBuf, 1024,

			"# vtk DataFile Version 2.0\n"
			"%s\n"							//Header name
			"%s\n"
			"DATASET STRUCTURED_POINTS\n"
			"DIMENSIONS %lli %lli %lli\n"
			"ORIGIN 0 0 0\n"
			"SPACING 1 1 1\n"
			"POINT_DATA %lli\n"
			"%s %s float\n"					//Vectors/Scalars; Dataset name
			"%s"							//Lookup table
			,							

				HeaderName,
				Bin ? "binary" : "ASCII",
				arrSize[arrTr[0]], arrSize[arrTr[1]], arrSize[arrTr[2]],
				nValues,
				(nVector > 1) ? "Vectors" : "Scalars",
				DatasetName,
				(nVector == 1) ? "LOOKUP_TABLE default\n" : "" );

		BufferedWrite(HeaderBuf, l);

	}

	bool _WriteArray(){

		long long nValues = arrSize[0]*arrSize[1]*arrSize[2];

		if(nValues == 0){
			SetErrorCode(0);
			return true;
		}

		if(HasHeader){

			_WriteVTKHeader();

		}

		const int vBufSz = 256;
		char vBuf[vBufSz];				//Current value buffer
		int vBufL = 0;

		int binTypeSz       = DataTypeSize(FileDataType);		//For bin files, number of bytes per value
		long long arrOffset = DataTypeSize(ArrayDataType);

		long long nValuesTot = nValues * nVector;	//Number of values x number of elements in each

		long long arrCount[4] = {0, 0, 0, 0};		//4D data counter

		long long arrSize4[4];			//4D array size
		arrSize4[0] = nVector;
		arrSize4[1] = arrSize[0];
		arrSize4[2] = arrSize[1];
		arrSize4[3] = arrSize[2];

		long long DataArrSize4[4];		//4D grid array size
		DataArrSize4[0] = nVector;
		DataArrSize4[1] = DataArrSize[0];
		DataArrSize4[2] = DataArrSize[1];
		DataArrSize4[3] = DataArrSize[2];

		void* arr = DataArray;

		for(long long nCount=0; nCount<nValuesTot; nCount++){

			//Obtain array element and convert to file data type

			long long n = DataFileBase::arrIndexFwdWrite(arrCount, arrSize4, DataArrSize4, DataArrOffset, Transform);

			void* arrPtr = (void*)&(((char*)arr)[n*arrOffset]);	//Array element to read from

			DataFileBase::TypeVar vArr;				//Value as array type
			DataFileBase::TypeVar vFile;			//Value converted to file type

			DataFileBase::DataToTypeVar(arrPtr, &vArr, ArrayDataType);
			DataFileBase::ConvertType(&vArr, &vFile, ArrayDataType, FileDataType);

			if(Bin){	//Binary output

				BufferedWrite(vFile.byte, binTypeSz);		//Write to file through buffer

			}else{		//ASCII output

				DataFileBase::ValToStr(vBuf, vBufSz, &vBufL, &vFile, FileDataType, nDecimalPlaces);

				vBuf[vBufL++] = WriteSpacer;

				BufferedWrite(vBuf, vBufL);

			}

		}

		SetErrorCode(0);

		return true;
	}

public:

		//Initialise class

	DataFileWriter(){

		long long BufSzDefault = 1048576;	//1MB read buffer

		Initialise(BufSzDefault);

	}

	DataFileWriter(long long BufferSize){

		Initialise(BufferSize);

	}


		//Open and close file

	bool Open(const char* FileName){

		return Open(FileName, false);

	}

	bool Open(const char* FileName, bool Append){

		Close();						//Close current file

		File = fopen(FileName, ( Append ? "ab" : "wb" ) );

		if(File == 0){
			SetErrorCode(1);			//File open failed
			return false;
		}

		SetErrorCode(0);

		return true;
	}
	
	void Close(){						//Close file

		if(File==0)
			return;

		FlushBuffer();

		fclose(File);
		File = 0;
	}

	
		//Write methods

	bool WriteArray(ArrayWriteOp* Operation){		//Read array from file

		if(!SetWriteOp(Operation))
			return false;

		return _WriteArray();

	}

	bool Write(void* Bytes, long long Count){

		BufferedWrite((char*)Bytes, Count);

		return true;
	}

	bool Write(const char* String){

		BufferedWrite(String, strlen(String));

		return true;
	}

	bool Write(char c){

		BufferedWrite(&c, 1);

		return true;
	}

	bool Write(float value, int nSigFigs){

		DataFileBase::TypeVar vType;		//Value as generic type
		DataFileBase::DataToTypeVar(&value, &vType, DataFileBase::Type_float);

		char Str[32];
		int Strl;

		DataFileBase::ValToStr(Str, sizeof(Str), &Strl, &vType, DataFileBase::Type_float, nSigFigs);

		BufferedWrite(Str, Strl);

		return true;
	}

	bool Write(float value){

		return Write(value, 4);

	}

	bool Write(double value, int nSigFigs){

		DataFileBase::TypeVar vType;		//Value as generic type
		DataFileBase::DataToTypeVar(&value, &vType, DataFileBase::Type_double);

		char Str[32];
		int Strl;

		DataFileBase::ValToStr(Str, sizeof(Str), &Strl, &vType, DataFileBase::Type_double, nSigFigs);

		BufferedWrite(Str, Strl);

		return true;
	}

	bool Write(double value){

		return Write(value, 4);

	}

	bool Write(int value){

		DataFileBase::TypeVar vType;		//Value as generic type
		DataFileBase::DataToTypeVar(&value, &vType, DataFileBase::Type_int32);

		char Str[32];
		int Strl;

		DataFileBase::ValToStr(Str, sizeof(Str), &Strl, &vType, DataFileBase::Type_int32, 0);

		BufferedWrite(Str, Strl);

		return true;
	}

	bool Write(unsigned int value){

		DataFileBase::TypeVar vType;		//Value as generic type
		DataFileBase::DataToTypeVar(&value, &vType, DataFileBase::Type_uint32);

		char Str[32];
		int Strl;

		DataFileBase::ValToStr(Str, sizeof(Str), &Strl, &vType, DataFileBase::Type_uint32, 0);

		BufferedWrite(Str, Strl);

		return true;
	}

	bool Write(long long value){

		DataFileBase::TypeVar vType;		//Value as generic type
		DataFileBase::DataToTypeVar(&value, &vType, DataFileBase::Type_int64);

		char Str[32];
		int Strl;

		DataFileBase::ValToStr(Str, sizeof(Str), &Strl, &vType, DataFileBase::Type_int64, 0);

		BufferedWrite(Str, Strl);

		return true;
	}

	bool Write(unsigned long long value){

		DataFileBase::TypeVar vType;		//Value as generic type
		DataFileBase::DataToTypeVar(&value, &vType, DataFileBase::Type_uint64);

		char Str[32];
		int Strl;

		DataFileBase::ValToStr(Str, sizeof(Str), &Strl, &vType, DataFileBase::Type_uint64, 0);

		BufferedWrite(Str, Strl);

		return true;
	}

	bool WriteVTKHeader(char* HeaderName, char* DataSetName, long long GridSz[3], bool Vectors){

		_WriteVTKHeader(HeaderName, DataSetName, GridSz, Vectors);

		return true;
	}

	
		//Obtain status

	const char* ErrorString(){

		int eCode = LastError();

		const char* Strings[] = {
			"The file was successfully written",
			"The file could not be opened",
			"Memory allocation error",
		};

		return Strings[eCode];
	}

	int LastError(){

		return ErrorCode;
	}


		//Finalise

	~DataFileWriter(){

		Close();

		Finalise();

	}

};