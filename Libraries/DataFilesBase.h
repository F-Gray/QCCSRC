#pragma once

#pragma warning( disable : 4996 )	//Potentially unsafe _snprintf and fopen

//Defines methods used by the data file classes

namespace DataFiles{

	struct ArrayRange{
		long long X0;
		long long X1;
		long long Y0;
		long long Y1;
		long long Z0;
		long long Z1;

		ArrayRange(long long x0, long long x1, long long y0, long long y1, long long z0, long long z1){
			X0 = x0;
			X1 = x1;
			Y0 = y0;
			Y1 = y1;
			Z0 = z0;
			Z1 = z1;
		}
	};


};

namespace DataFileBase{

	enum DataType{
		Type_bool,
		Type_char,
		Type_uchar,
		Type_int16,
		Type_uint16,
		Type_int32,
		Type_uint32,
		Type_int64,
		Type_uint64,
		Type_float,
		Type_double
	};

	int DataTypeSize(DataType Type){

		const int TypeLength[] = { 
			sizeof(bool),
			sizeof(char),
			sizeof(unsigned char),
			sizeof(short),
			sizeof(unsigned short),
			sizeof(int),
			sizeof(unsigned int),
			sizeof(long long),
			sizeof(unsigned long long),
			sizeof(float),
			sizeof(double)
		};

		return TypeLength[Type];
	}

	enum Rotation{
		RotateGrid_XYZ,
		RotateGrid_YZX,
		RotateGrid_ZXY,
		RotateGridAndVector_XYZ,
		RotateGridAndVector_YZX,
		RotateGridAndVector_ZXY
	};
	
	union TypeVar{
		char byte[8];
		bool bool_1;
		char char_1;
		unsigned char uchar_1;
		unsigned short uint_2;
		signed   short int_2;
		unsigned int uint_4;
		signed   int int_4;
		unsigned long long uint_8;
		signed   long long int_8;
		float float_4;
		double double_8;
	};

	////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Data type conversions

	int DataTypeSize(DataType Type);			//Returns the size, in bytes, of the DataType Type
	int DataTypeSizeBits(DataType Type);		//Returns the size, in bits, of the DataType Type

	void DataToTypeVar(void* DataIn, TypeVar* DataOut, DataType TypeIn);						//Create generic TypeVar from DataIn of type TypeIn
	void TypeVarToData(void* DataOut, TypeVar* DataIn, DataType TypeOut);						//Write TypeVar to DataOut with type TypeOut
	void ConvertType(TypeVar* DataIn, TypeVar* DataOut, DataType TypeIn, DataType TypeOut);		//Convert TypeVar to different DataType
	
	bool StrToVal(char* Str, bool _bin, void* arrPtr, DataType _fileType, DataType _arrType);	//Convert ASCII string from _fileType to _arrType
	bool ValToStr(char* Str, int StrSz, int* nBytes, TypeVar* DataIn, DataType TypeIn, int nDPFlt);	//Convert DataIn of type TypeIn to string

	template<class FltType>
	int FltToStr(FltType f, char* Str, int StrSz, int nDP);		//Converts float type to string

	//Array index

	long long arrIndexFwd(long long* nCount, long long* arrSize, Rotation Transform);	//Obtains index of array and increments the 4D counter for array of size arrSize under rotation Transform
	long long arrIndexFwd(long long* nCount, long long* arrSize, long long* arrGridSize, long long* GridOffset, Rotation Transform);

	long long arrIndexBack(long long* nCount, long long* arrSize, Rotation Transform);	//Obtains index of array and decrements the 4D counter for array of size arrSize under rotation Transform
	long long arrIndexBack(long long* nCount, long long* arrSize, long long* arrGridSize, long long* GridOffset, Rotation Transform);

	long long arrIndexFwdWrite(long long* nCount, long long* arrSize, Rotation Transform);
	long long arrIndexFwdWrite(long long* nCount, long long* arrSize, long long* arrGridSize, long long* GridOffset, Rotation Transform);

	////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Array index

		//Transforms file to array

	long long arrIndexBack(long long* nCount, long long* arrSize, Rotation Transform){

		long long Offset[3] = {0, 0, 0};
		
		return arrIndexBack(nCount, arrSize, arrSize, Offset, Transform);
	}

	long long arrIndexBack(long long* nCount, long long* arrSize, long long* arrGridSize, long long* GridOffset, Rotation Transform){

		bool RotateXYZ = (Transform == RotateGrid_XYZ) || (Transform == RotateGridAndVector_XYZ);
		bool RotateYZX = (Transform == RotateGrid_YZX) || (Transform == RotateGridAndVector_YZX);

		// nCount[0] = n
		// nCount[1] = x
		// nCount[2] = y
		// nCount[3] = z

		//Rotate vector

		long long vi = nCount[0];
		
		if(Transform == RotateGridAndVector_YZX) vi = nCount[0]-1;		
		if(Transform == RotateGridAndVector_ZXY) vi = nCount[0]-2;

		while(vi < 0) vi += arrSize[0];

		//Transform according to offset

		long long nCountGrid[4];

		nCountGrid[1] = nCount[1] + GridOffset[0];	//X
		nCountGrid[2] = nCount[2] + GridOffset[1];	//Y
		nCountGrid[3] = nCount[3] + GridOffset[2];	//Z

		//Index

		long long n;

		if(RotateXYZ){																									// XYZ -> XYZ ; Z{ Y{ X{n} } } 
			n = arrGridSize[0]*(arrGridSize[1]*(arrGridSize[2]*nCountGrid[3] + nCountGrid[2]) + nCountGrid[1]) + vi;
		}else if(RotateYZX){																							// XYZ -> YZX ; X{ Z{ Y{n} } } 
			n = arrGridSize[0]*(arrGridSize[2]*(arrGridSize[3]*nCountGrid[1] + nCountGrid[3]) + nCountGrid[2]) + vi;
		}else{																											// XYZ -> ZXY ; Y{ X{ Z{n} } }
			n = arrGridSize[0]*(arrGridSize[3]*(arrGridSize[1]*nCountGrid[2] + nCountGrid[1]) + nCountGrid[3]) + vi;
		}

		nCount[0]--;
		if(nCount[0]<0){
			nCount[0] = arrSize[0]-1;
			nCount[1]--;
			if(nCount[1]<0){
				nCount[1] = arrSize[1]-1;
				nCount[2]--;
				if(nCount[2]<0){
					nCount[2] = arrSize[2]-1;
					nCount[3]--;
				}
			}
		}

		return n;
	}

	long long arrIndexFwd(long long* nCount, long long* arrSize, Rotation Transform){

		long long Offset[3] = {0, 0, 0};

		return arrIndexFwd(nCount, arrSize, arrSize, Offset, Transform);
	}

	long long arrIndexFwd(long long* nCount, long long* arrSize, long long* arrGridSize, long long* GridOffset, Rotation Transform){

		bool RotateXYZ = (Transform == RotateGrid_XYZ) || (Transform == RotateGridAndVector_XYZ);
		bool RotateYZX = (Transform == RotateGrid_YZX) || (Transform == RotateGridAndVector_YZX);

		// nCount[0] = n
		// nCount[1] = x
		// nCount[2] = y
		// nCount[3] = z

		//Rotate vector

		long long vi = nCount[0];
		
		if(Transform == RotateGridAndVector_YZX) vi = nCount[0]-1;		
		if(Transform == RotateGridAndVector_ZXY) vi = nCount[0]-2;

		while(vi < 0) vi += arrSize[0];

		//Transform according to offset

		long long nCountGrid[4];

		nCountGrid[1] = nCount[1] + GridOffset[0];	//X
		nCountGrid[2] = nCount[2] + GridOffset[1];	//Y
		nCountGrid[3] = nCount[3] + GridOffset[2];	//Z

		//Index

		long long n;

		if(RotateXYZ){																									// XYZ -> XYZ ; Z{ Y{ X{n} } } 
			n = arrGridSize[0]*(arrGridSize[1]*(arrGridSize[2]*nCountGrid[3] + nCountGrid[2]) + nCountGrid[1]) + vi;
		}else if(RotateYZX){																							// XYZ -> YZX ; X{ Z{ Y{n} } } 
			n = arrGridSize[0]*(arrGridSize[2]*(arrGridSize[3]*nCountGrid[1] + nCountGrid[3]) + nCountGrid[2]) + vi;
		}else{																											// XYZ -> ZXY ; Y{ X{ Z{n} } }
			n = arrGridSize[0]*(arrGridSize[3]*(arrGridSize[1]*nCountGrid[2] + nCountGrid[1]) + nCountGrid[3]) + vi;
		}

		nCount[0]++;
		if(nCount[0]>=arrSize[0]){		//Vector size
			nCount[0] = 0;
			nCount[1]++;
			if(nCount[1]>=arrSize[1]){		//X component
				nCount[1] = 0;
				nCount[2]++;
				if(nCount[2]>=arrSize[2]){		//Y component
					nCount[2] = 0;
					nCount[3]++;
				}
			}
		}

		return n;
	}

		//Transforms array to file

	long long arrIndexFwdWrite(long long* nCount, long long* arrSize, Rotation Transform){

		long long Offset[3] = {0, 0, 0};
		
		return arrIndexFwdWrite(nCount, arrSize, arrSize, Offset, Transform);
	}

	long long arrIndexFwdWrite(long long* nCount, long long* arrSize, long long* arrGridSize, long long* GridOffset, Rotation Transform){

		bool RotateXYZ = (Transform == RotateGrid_XYZ) || (Transform == RotateGridAndVector_XYZ);
		bool RotateYZX = (Transform == RotateGrid_YZX) || (Transform == RotateGridAndVector_YZX);

		int arrTr[3];	//Transform array dimensions to file dimensions
		arrTr[0] = RotateXYZ ? 1 : ( RotateYZX ? 2 : 3 );
		arrTr[1] = RotateXYZ ? 2 : ( RotateYZX ? 3 : 1 );
		arrTr[2] = RotateXYZ ? 3 : ( RotateYZX ? 1 : 2 );

		// nCount[0] = n
		// nCount[1] = x
		// nCount[2] = y
		// nCount[3] = z

		//Rotate vector

		long long vi = nCount[0];
		
		if(Transform == RotateGridAndVector_YZX) vi = nCount[0]+1;		
		if(Transform == RotateGridAndVector_ZXY) vi = nCount[0]+2;

		while(vi >= arrSize[0]) vi -= arrSize[0];

		//Index
		
		long long n;

		if(RotateXYZ){																																//Z{ Y{ X{n} } }
			n = arrGridSize[0]*(arrGridSize[1]*(arrGridSize[2]*(nCount[3] + GridOffset[2]) + (nCount[2] + GridOffset[1])) + (nCount[1] + GridOffset[0])) + vi;
		}else if(RotateYZX){																														//X{ Z{ Y{n} } }
			n = arrGridSize[0]*(arrGridSize[1]*(arrGridSize[2]*(nCount[2] + GridOffset[2]) + (nCount[1] + GridOffset[1])) + (nCount[3] + GridOffset[0])) + vi;
		}else{																																		//Y{ X{ Z{n} } }
			n = arrGridSize[0]*(arrGridSize[1]*(arrGridSize[2]*(nCount[1] + GridOffset[2]) + (nCount[3] + GridOffset[1])) + (nCount[2] + GridOffset[0])) + vi;
		}

		nCount[0]++;
		if(nCount[0]>=arrSize[0]){
			nCount[0] = 0;
			nCount[1]++;
			if(nCount[1]>=arrSize[arrTr[0]]){
				nCount[1] = 0;
				nCount[2]++;
				if(nCount[2]>=arrSize[arrTr[1]]){
					nCount[2] = 0;
					nCount[3]++;
				}
			}
		}

		return n;
	}

	//Data type conversions

	void ConvertType(TypeVar* DataIn, TypeVar* DataOut, DataType TypeIn, DataType TypeOut){

		if(TypeIn == Type_bool){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = DataIn->bool_1;
			else if(TypeOut == Type_char)   DataOut->char_1 = DataIn->char_1;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = DataIn->uchar_1;
			else if(TypeOut == Type_int16)  DataOut->int_4	= DataIn->bool_1 ? 1 : 0;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= DataIn->bool_1 ? 1 : 0;
			else if(TypeOut == Type_int32)  DataOut->int_4	= DataIn->bool_1 ? 1 : 0;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= DataIn->bool_1 ? 1 : 0;
			else if(TypeOut == Type_int64)  DataOut->int_8	= DataIn->bool_1 ? 1 : 0;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = DataIn->bool_1 ? 1 : 0;
			else if(TypeOut == Type_float)  DataOut->float_4  = DataIn->bool_1 ? 1.0f : 0;
			else if(TypeOut == Type_double) DataOut->double_8 = DataIn->bool_1 ? 1.0 : 0;
		}else
		if(TypeIn == Type_char){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->char_1!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = DataIn->char_1;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = DataIn->uchar_1;
			else if(TypeOut == Type_int16)  DataOut->int_4	= (short)DataIn->char_1;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= (unsigned short)DataIn->char_1;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->char_1;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->char_1;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->char_1;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->char_1;
			else if(TypeOut == Type_float)  DataOut->float_4 = (float)DataIn->char_1;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->char_1;
		}else
		if(TypeIn == Type_uchar){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->uchar_1!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)DataIn->uchar_1;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = DataIn->uchar_1;
			else if(TypeOut == Type_int16)  DataOut->int_4	= (short)DataIn->uchar_1;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= (unsigned short)DataIn->uchar_1;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->uchar_1;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->uchar_1;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->uchar_1;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->uchar_1;
			else if(TypeOut == Type_float)  DataOut->float_4 = (float)DataIn->uchar_1;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->uchar_1;
		}else
		if(TypeIn == Type_int16){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->int_2!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)DataIn->int_2;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = (unsigned char)DataIn->int_2;
			else if(TypeOut == Type_int16)  DataOut->int_4	= DataIn->int_2;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= (unsigned short)DataIn->int_2;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->int_2;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->int_2;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->int_2;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->int_2;
			else if(TypeOut == Type_float)  DataOut->float_4  = (float)DataIn->int_2;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->int_2;
		}else
		if(TypeIn == Type_uint16){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->uint_2!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)DataIn->uint_2;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = (unsigned char)DataIn->uint_2;
			else if(TypeOut == Type_int16)  DataOut->int_4	= (short)DataIn->uint_2;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= DataIn->uint_2;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->uint_2;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->uint_2;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->uint_2;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->uint_2;
			else if(TypeOut == Type_float)  DataOut->float_4  = (float)DataIn->uint_2;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->uint_2;
		}else
		if(TypeIn == Type_int32){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->int_4!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)DataIn->int_4;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = (unsigned char)DataIn->int_4;
			else if(TypeOut == Type_int16)  DataOut->int_4	= (short)DataIn->int_4;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= (unsigned short)DataIn->int_4;
			else if(TypeOut == Type_int32)  DataOut->int_4	= DataIn->int_4;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->int_4;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->int_4;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->int_4;
			else if(TypeOut == Type_float)  DataOut->float_4  = (float)DataIn->int_4;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->int_4;
		}else
		if(TypeIn == Type_uint32){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->uint_4!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)DataIn->uint_4;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = (unsigned char)DataIn->uint_4;
			else if(TypeOut == Type_int16)  DataOut->int_4	= (short)DataIn->uint_4;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= (unsigned short)DataIn->uint_4;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->uint_4;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= DataIn->uint_4;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->uint_4;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->uint_4;
			else if(TypeOut == Type_float)  DataOut->float_4  = (float)DataIn->uint_4;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->uint_4;
		}else
		if(TypeIn == Type_int64){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->int_8!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)DataIn->int_8;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = (unsigned char)DataIn->int_8;
			else if(TypeOut == Type_int16)  DataOut->int_4	= (short)DataIn->int_8;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= (unsigned short)DataIn->int_8;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->int_8;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->int_8;
			else if(TypeOut == Type_int64)  DataOut->int_8	= DataIn->int_8;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->int_8;
			else if(TypeOut == Type_float)  DataOut->float_4 = (float)DataIn->int_8;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->int_8;
		}else
		if(TypeIn == Type_uint64){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->uint_8!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)DataIn->uint_8;
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = (unsigned char)DataIn->uint_8;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (short)DataIn->uint_8;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned short)DataIn->uint_8;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->uint_8;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->uint_8;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->uint_8;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = DataIn->uint_8;
			else if(TypeOut == Type_float)  DataOut->float_4  = (float)DataIn->uint_8;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->uint_8;
		}else
		if(TypeIn == Type_float){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->float_4!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)((int)DataIn->float_4);
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = (unsigned char)((int)DataIn->float_4);
			else if(TypeOut == Type_int16)  DataOut->int_4	= (short)DataIn->float_4;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= (unsigned short)DataIn->float_4;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->float_4;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->float_4;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->float_4;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->float_4;
			else if(TypeOut == Type_float)  DataOut->float_4  = DataIn->float_4;
			else if(TypeOut == Type_double) DataOut->double_8 = (double)DataIn->float_4;
		}else
		if(TypeIn == Type_double){
				 if(TypeOut == Type_bool)   DataOut->bool_1 = (DataIn->double_8!=0);
			else if(TypeOut == Type_char)   DataOut->char_1 = (char)((int)DataIn->double_8);
			else if(TypeOut == Type_uchar)  DataOut->uchar_1 = (unsigned char)((int)DataIn->double_8);
			else if(TypeOut == Type_int16)  DataOut->int_4	= (short)DataIn->double_8;
			else if(TypeOut == Type_uint16) DataOut->uint_4	= (unsigned short)DataIn->double_8;
			else if(TypeOut == Type_int32)  DataOut->int_4	= (int)DataIn->double_8;
			else if(TypeOut == Type_uint32) DataOut->uint_4	= (unsigned int)DataIn->double_8;
			else if(TypeOut == Type_int64)  DataOut->int_8	= (long long)DataIn->double_8;
			else if(TypeOut == Type_uint64) DataOut->uint_8 = (unsigned long long)DataIn->double_8;
			else if(TypeOut == Type_float)  DataOut->float_4  = (float)DataIn->double_8;
			else if(TypeOut == Type_double) DataOut->double_8 = DataIn->double_8;
		}

	}

	void DataToTypeVar(void* DataIn, TypeVar* DataOut, DataType TypeIn){

			 if(TypeIn == Type_bool)   DataOut->bool_1 = *((bool*)DataIn);
		else if(TypeIn == Type_char)   DataOut->char_1 = *((char*)DataIn);
		else if(TypeIn == Type_uchar)  DataOut->uchar_1 = *((unsigned char*)DataIn);
		else if(TypeIn == Type_int16)  DataOut->int_2  = *((short*)DataIn);
		else if(TypeIn == Type_uint16) DataOut->uint_2 = *((unsigned short*)DataIn);
		else if(TypeIn == Type_int32)  DataOut->int_4  = *((int*)DataIn);
		else if(TypeIn == Type_uint32) DataOut->uint_4 = *((unsigned int*)DataIn);
		else if(TypeIn == Type_int64)  DataOut->int_8  = *((long long*)DataIn);
		else if(TypeIn == Type_uint64) DataOut->uint_8 = *((unsigned long long*)DataIn);
		else if(TypeIn == Type_float)  DataOut->float_4  = *((float*)DataIn);
		else if(TypeIn == Type_double) DataOut->double_8 = *((double*)DataIn);

	}

	void TypeVarToData(void* DataOut, TypeVar* DataIn, DataType TypeOut){

			 if(TypeOut == Type_bool)   *((bool*)DataOut) = DataIn->bool_1;
		else if(TypeOut == Type_char)   *((char*)DataOut) = DataIn->char_1;
		else if(TypeOut == Type_uchar)  *((unsigned char*)DataOut)		= DataIn->uchar_1;
		else if(TypeOut == Type_int16)  *((short*)DataOut)				= DataIn->int_2;
		else if(TypeOut == Type_uint16) *((unsigned short*)DataOut)		= DataIn->uint_2;
		else if(TypeOut == Type_int32)  *((int*)DataOut)				= DataIn->int_4;
		else if(TypeOut == Type_uint32) *((unsigned int*)DataOut)		= DataIn->uint_4;
		else if(TypeOut == Type_int64)  *((long long*)DataOut)			= DataIn->int_8;
		else if(TypeOut == Type_uint64) *((unsigned long long*)DataOut) = DataIn->uint_8;
		else if(TypeOut == Type_float)  *((float*)DataOut)	= DataIn->float_4;
		else if(TypeOut == Type_double) *((double*)DataOut) = DataIn->double_8;

	}

	bool StrToVal(char* Str, bool _bin, void* arrPtr, DataType _fileType, DataType _arrType){

		int l = (int)strlen(Str);

		int arrSz  = DataTypeSize(_arrType);
		int fdatSz = DataTypeSize(_fileType);

		TypeVar vBin;	//Data in _fileType format

		//Interpret as _fileType
		
		if(_bin){	//Binary

			for(int i=0; i<fdatSz; i++)
				vBin.byte[i] = Str[i];

		}else{		//ASCII

			if(_fileType == Type_bool)		  vBin.bool_1 = (atoi(Str)!=0);
			else if(_fileType == Type_char)   vBin.char_1 = Str[0];
			else if(_fileType == Type_uchar)  vBin.uchar_1 = (unsigned char)atoi(Str);
			else if(_fileType == Type_int16)  vBin.int_2  = (short)atoi(Str);
			else if(_fileType == Type_uint16) vBin.uint_2 = (unsigned short)atoi(Str);
			else if(_fileType == Type_int32)  vBin.int_4  = atoi(Str);
			else if(_fileType == Type_uint32) vBin.uint_4 = (unsigned int)atoi(Str);
			else if(_fileType == Type_int64)  vBin.int_8  = (long long)atoll(Str);
			else if(_fileType == Type_uint64) vBin.uint_8 = (unsigned long long)atoll(Str);
			else if(_fileType == Type_float)  vBin.float_4  = (float)atof(Str);
			else if(_fileType == Type_double) vBin.double_8 = atof(Str);

		}

		//Convert from _fileType to _arrType

		TypeVar vArr;

		ConvertType(&vBin, &vArr, _fileType, _arrType);		//Convert from file type to array type

		TypeVarToData(arrPtr, &vArr, _arrType);	//Put data into array element
		
		return true;

	}

	bool ValToStr(char* Str, int StrSz, int* nBytes, TypeVar* DataIn, DataType TypeIn, int nDPFlt){

		int n = 1;

		switch(TypeIn){

			case Type_bool:
				Str[0] = DataIn->bool_1 ? '1' : '0';
				break;
			case Type_char:
				Str[0] = DataIn->char_1;
				break;
			case Type_uchar:
				n = _snprintf(Str, StrSz, "%i", (int)DataIn->uchar_1);
				break;
			case Type_int16:
				n = _snprintf(Str, StrSz, "%i", DataIn->int_2);
				break;
			case Type_uint16:
				n = _snprintf(Str, StrSz, "%u", DataIn->uint_2);
				break;
			case Type_int32:
				n = _snprintf(Str, StrSz, "%i", DataIn->int_4);
				break;
			case Type_uint32:
				n = _snprintf(Str, StrSz, "%u", DataIn->uint_4);
				break;
			case Type_int64:
				n = _snprintf(Str, StrSz, "%lli", DataIn->int_8);
				break;
			case Type_uint64:
				n = _snprintf(Str, StrSz, "%llu", DataIn->uint_8);
				break;
			case Type_float:
				n = FltToStr(DataIn->float_4, Str, StrSz, nDPFlt);
				break;
			case Type_double:
				n = FltToStr(DataIn->double_8, Str, StrSz, nDPFlt);
				break;

		}

		*nBytes = n;

		return true;
	}

	template<class FltType>
	int FltToStr(FltType f, char* Str, int StrSz, int nDP){

		int n;

		if(f==0){
			Str[0] = '0';
			n = 1;
		}else{
			n = _snprintf(Str, StrSz, "%.*g", nDP, f);
		}

		return n;
	}	//*/

/*	template<class FltType>
	int FltToStr(FltType f, char* Str, int StrSz, int nDP){

		int n;


		n = _snprintf(Str, StrSz, "%.*e", nDP, f);
		

		return n;
	}	//*/

};