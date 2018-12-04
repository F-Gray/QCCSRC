#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

//Input file reading

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

struct InputValueStruct{
	const char* DataName;	//Name in input file
	void* Ptr;				//Ptr to variable
	int   PtrLength;		//Max data length
	int   DataType;			//One of datatype values eg DataType_Bool
	bool  Required;			//Throws error if value not found in input file
};

struct InputValueStructEx{
	const char* DataName;				//Name in input file
	void* Ptr;							//Ptr to variable
	int   DataSz;						//Max data length of one data entry
	int   Count;						//Max number of values
	int   DataType;						//One of datatype values eg DataType_Bool
	InputValueStructEx* DataStruct;		//If type is DataType_InputValueStructEx, Ptr to another InputValueStructEx defining the data within the class at Ptr
	int   DataStructCount;				//If type is DataType_InputValueStructEx, the number of entries in DataStruct
	int*  CountRead;					//Returns the number of values read in
	bool  Required;						//Throws error if value not found in input file
	const char* Description;			//Description of the field
};

int ReadInputFile(const char* InputFileName, InputValueStruct* Values, int nValues);

/* Usage

1) Create an array of InputValueStruct with values to be read in eg:

-----------------------------------------------------------------------------------------------------------------

InputValueStruct InputFileData[] = {		//Defines entries in the input file to be read in

	//	{ [Property Name], [Variable Ptr], [Data Length], [Data Type (eg DataType_Int)], [must be specified?] },

		{ "NumberOfThreads", &LocalThreadNum, sizeof(LocalThreadNum), DataType_Int, false },

		{ "NLatticeX", &NLatticeX, sizeof(NLatticeX), DataType_Int, true },
		{ "NLatticeY", &NLatticeY, sizeof(NLatticeY), DataType_Int, true },
		{ "NLatticeZ", &NLatticeZ, sizeof(NLatticeZ), DataType_Int, true },

		{ "ResolutionMicron", &Resolution, sizeof(Resolution), DataType_Double, true },
		{ "Viscosity", &Viscosity, sizeof(Viscosity), DataType_Double, true },

};

-----------------------------------------------------------------------------------------------------------------

2) Read in using

	int ReadInputFile(char* InputFileName, InputValueStruct* Values, int nValues);

3) Return values:

	0 - read in succesfully
	1 - file could not be opened
	2 - some required values were not read in

4) [If using MPI] Distribute input file data amongst threads

void DistributeInputFileMPI(InputValueStruct* Values, int nValues, int SrcThread, MPI_Comm Comm)

*/


//Class

class InputFileReader{

public:
	enum InputFileStatus{

		Success,
		FileOpenFailed,
		RequiredValuesNotRead
	};

private:
	int Mode;

#ifdef MPI_INCLUDED

	int _SrcThread;
	MPI_Comm MPIComm;

#endif

	InputFileStatus Status;
	vector<int> InDataSet;

public:
	InputFileReader(){

		memset(this, 0, sizeof(InputFileReader));

		Mode = 0;
		Status = Success;
	}

#ifdef MPI_INCLUDED
	InputFileReader(MPI_Comm Comm, int SrcThread){

		Mode = 1;		//MPI mode

		_SrcThread = SrcThread;
		MPIComm = Comm;
	}
#endif

	~InputFileReader(){

	}

	bool ReadInputFile(const char* InputFileName, InputValueStructEx* Values, int nValues){

		fflush(stdout);

		Status = Success;

		bool ret;

#ifdef MPI_INCLUDED

		if(Mode == 1){

			int tId;
			MPI_Comm_rank(MPIComm, &tId);	//Get thread ID (index from 0 -> nThreads-1)

			if(tId == _SrcThread){

				ret = _ReadInputFile(InputFileName, Values, nValues);

			}

			DistributeInputFileMPI(Values, nValues, 0);

			MPI_Bcast(&ret, 1, MPI_BYTE, _SrcThread, MPIComm);
			MPI_Bcast(&Status, sizeof(Status), MPI_BYTE, _SrcThread, MPIComm);

		}

#endif

		if(Mode == 0){

			ret = _ReadInputFile(InputFileName, Values, nValues);

		}

		return ret;
	}

	InputFileStatus GetErrorCode(){

		return Status;
	}

	const char* GetErrorStr(InputFileStatus Status){

		const char* Str[] = {
			"The operation completed successfully",
			"The file could not be opened",
			"Required values were not found",
		};

		return Str[(int)Status];
	}

private:

	bool _ReadInputFile(const char* InputFileName, InputValueStructEx* Values, int nValues){

		if(nValues <= 0)
			return true;

		ifstream InFile;
		InFile.open(InputFileName);
		InFile.unsetf(ios::skipws);

		if(InFile.fail()){
			Status = FileOpenFailed;
			return false;				//File couldn't be opened
		}

		InDataSet = vector<int>(nValues, 0);		//Whether each value has been set

		string Buf;

		bool ret = ReadData(InFile, Buf, true, Values, nValues, InDataSet, 0);

		bool AllRead = true;

		for(int i=0; i<nValues; i++){

		//	printf("Entry '%s' %s\n", Values[i].DataName, InDataSet[i] ? "set" : "not set");

			if(!InDataSet[i] && Values[i].Required){
				AllRead = false;
				Status = RequiredValuesNotRead;
			}

		}

		return (ret && AllRead);
	}

	bool ReadData(ifstream& InFile, string& Data, bool FMode, InputValueStructEx* Values, int nValues, vector<int>& InDataRead, size_t BaseAddr){

		string Buf;
		string Property;
		string Value;

		Buf.reserve(1024);
		Property.reserve(256);
		Value.reserve(1024);

		int pos = 0;

		//Obtain entry

		while(true){

			Buf.clear();
			Property.clear();
			Value.clear();

			int nBrackets = 0;
			bool InQuotes = 0;
			bool Comment = false;

			while(true){

				//Next character

				char c;

				if(FMode){
					if(!InFile.get(c)) break;
				}else{
					if(pos<Data.length()) c = Data[pos++];
					else break;
				}

				if(!Comment && !InQuotes && c == '/'){	//Comment "//"

					char c2 = ' ';

					if(FMode){
						InFile.get(c2);
						InFile.unget();
					}else{
						if(pos+1 < Data.length()) c2 = Data[pos+1];
					}

					if(c2 == '/')
						Comment = true;
				}

				if(c == '"' && !Comment)
					InQuotes = !InQuotes;

				if(c == '\n'){

					Comment = false;
					InQuotes = false;

					if(nBrackets == 0)		//End of entry
						break;
				}

				if(c == '{' && !InQuotes && !Comment)
					nBrackets++;

				if(c == '}' && !InQuotes && !Comment)
					nBrackets = max(nBrackets - 1, 0);

				Buf.append(1, c);
			}

		/*	printf("Obtain from buf: %s\n", Buf.c_str());
			fflush(stdout);		//*/

			int r = ObtainValues(Buf, Property, Value);

			if(r > 0){
				if((FMode && InFile.eof()) || (!FMode && pos==Data.length())) break;		//End of file/buffer
				continue;
			}

			int Ind = FindEntry(Property, Values, nValues);

			if(Ind == -1){
				if((FMode && InFile.eof()) || (!FMode && pos==Data.length())) break;		//End of file/buffer
				continue;
			}

		/*	printf("---------------\n");
			printf("%s", Buf.c_str());
			printf("\n");
			printf("Name: '%s' [%i]\n", Property.c_str(), Ind);
			printf("Data: %s\n", Value.c_str());	
			fflush(stdout);	//*/

			int Count = InterpretValue(Value, &Values[Ind], BaseAddr);

			if(Count > 0){
				InDataRead[Ind] = 1;
			}
		}

		return true;
	}

	int ObtainValues(string& Buf, string& Property, string& Value){

		size_t SepInd = Buf.find_first_of("=\n");

		if(SepInd == string::npos) return 2;	//Invalid line
		if(Buf[SepInd] == '\n')    return 2;

		//Property name

		string Prop = Buf.substr(0, SepInd);
			
		size_t NameStart = Prop.find_first_not_of(" \t");
		size_t NameEnd   = Prop.find_first_of(" \t=", NameStart);

		if(NameEnd-NameStart <= 0) return 2;

		Property = Prop.substr(NameStart, NameEnd-NameStart);

		//Data

		size_t VStart = Buf.find_first_not_of(" \t", SepInd+1);

		if(VStart == string::npos) return 2;
		if(Buf.length()-VStart <= 0) return 2;

		Value = Buf.substr(VStart, Buf.length() - VStart);

		return 0;
	}

	int FindEntry(string& Name, InputValueStructEx* Values, int nValues){

		for(int i=0; i<nValues; i++){

			if(Name.compare(Values[i].DataName) == 0)
				return i;
		}

		return -1;
	}

	int InterpretValue(string& Value, InputValueStructEx* ValueStruct, size_t BaseAddr){

		int Count = 0;

		if(ValueStruct->DataType != DataType_InputValueStructEx){	//Simple entry

			string Element;

			for(int Pos=0, EndPos; ReadElement(Value, Element, Pos, &EndPos) == 0; Pos = EndPos){		//Read next element

				//Value

				char* Ptr = BaseAddr + (char*)ValueStruct->Ptr + (size_t)ValueStruct->DataSz*Count;		//Entry address
				
				switch(ValueStruct->DataType){

					case DataType_Short:			//Interpret as short
						*((short*)Ptr) = (short)atoi(Element.c_str());
						break;
					case DataType_Bool:				//Interpret as bool
						*((bool*)Ptr) = ReadBool(Element.c_str());
						break;
					case DataType_Int:				//Interpret as int
						*((int*)Ptr) = atoi(Element.c_str());
						break;
					case DataType_LongLongInt:		//Interpret as long long
						*((long long*)Ptr) = atoll(Element.c_str());
						break;
					case DataType_Float:			//Interpret as float
						*((float*)Ptr) = (float)atof(Element.c_str());
						break;
					case DataType_Double:			//Interpret as double
						*((double*)Ptr) = atof(Element.c_str());
						break;
					case DataType_String:			//Copy over string
						memcpy((char*)Ptr, Element.c_str(), min(ValueStruct->DataSz-1, (int)Element.length()));
						((char*)Ptr)[ min(ValueStruct->DataSz-1, (int)Element.length()) ] = '\0';
						break;
				}

				Count++;

				if(Count == ValueStruct->Count)
					break;
			}

			if(ValueStruct->CountRead != 0){
				char* CountPtr = BaseAddr + (char*)ValueStruct->CountRead;
				*((int*)CountPtr) = Count;
			}

		}else{			//Sub struct of values

			string SubValue;
			
			int nBrackets = 0;

			int DataLevel = -1;
			int DataPos = -1;
			int DataLen;

			int Len = (int)Value.length();

			int Pos = 0;

			while(Count < ValueStruct->Count){

				int ExitCode = 1;

				for(int i=Pos; i<Len; i++){

					if(Value[i] == '/' && i<Len-1 && Value[i+1] == '/'){	//Comment "//"

						size_t EoL = Value.find_first_of("\r\n", i+2);

						if(EoL == string::npos)
							break;

						i = (int)EoL + 1;	//Skip comment

					}else if(Value[i] == '{'){

						nBrackets++;

						if(nBrackets == DataLevel){
							DataPos = i+1;
						}

					}else if(Value[i] == '}'){

						nBrackets = max(0, nBrackets-1); 

						if(nBrackets == DataLevel-1){
							DataLen = i - DataPos;
							Pos = i+1;				//Restart character
							ExitCode = 0;
							break;
						}

					}else if(DataLevel == -1 && Value[i] != ' ' && Value[i] != '\t' && Value[i] != '\r' && Value[i] != '\n'){		//First significant character

						DataPos = i;
						DataLevel = nBrackets;
					}
				}

				if(ExitCode == 1)
					break;

				SubValue = Value.substr(DataPos, DataLen);

				vector<int> InDataSet(ValueStruct->DataStructCount, 0);		//Whether each value has been set

				ifstream InFile;

			//	printf("SubValue: %s\n", SubValue.c_str());

				size_t Addr = (size_t)ValueStruct->Ptr + (size_t)ValueStruct->DataSz*Count;

				ReadData(InFile, SubValue, false, ValueStruct->DataStruct, ValueStruct->DataStructCount, InDataSet, Addr);

				Count++;
			}

			if(ValueStruct->CountRead != 0){
				char* CountPtr = BaseAddr + (char*)ValueStruct->CountRead;
				*((int*)CountPtr) = Count;
			}
		}

		return Count;
	}

	int ReadElement(string& Value, string& Element, int StartPos, int* EndPos){

		bool InQuotes = false;
		int Stage = 0;			// 0 = searching for start of value; 1 = reading value

		int ReadPos = 0;

		size_t Len = Value.length();

		for(int i=StartPos; i<Len; i++){

			char c = Value[i];

		//	printf("Value: %c\n", c);

			if(!InQuotes && Value[i] == '/' && i<Len-1 && Value[i+1] == '/'){	//Comment "//"			

				if(Stage == 1){			//End of value
					*EndPos = i;
					Element = Value.substr(ReadPos, i-ReadPos);
					break;
				}

				size_t EoL = Value.find_first_of("\r\n", i+2);

				if(EoL == string::npos)
					break;

				i = (int)EoL + 1;	//Skip comment

			}else if(c == '\"'){		//Quotation mark

				if(Stage == 0){

					InQuotes = true;
					Stage = 1;				//Reading value
					ReadPos = i+1;			//First character of value

				}else if(Stage == 1){		//End of value

					*EndPos = i+1;
					Element = Value.substr(ReadPos, i-ReadPos);
					break;
				}

			}else if((c == '\r' || c == '\n') || (!InQuotes && (c==',' || c==';' || c=='\t' || c==' ' || c=='{' || c=='}'))){		//End of line or list separator

				if(Stage == 1){		//End of element
					*EndPos = i;
					Element = Value.substr(ReadPos, i-ReadPos);
					break;
				}

			}else if(i == Len-1){		//Last character

				if(Stage == 0){
					ReadPos = i;	//First and last character of value
					Stage = 1;
				}

				if(Stage == 1){
					*EndPos = i+1;
					Element = Value.substr(ReadPos, i+1-ReadPos);
					break;
				}

			}else{

				if(Stage == 0){
					Stage = 1;
					ReadPos = i;			//First character of value
				}
			}
		}

		if(Stage == 1)		//Value read in
			return 0;

		return 1;	//Value not read
	}

	bool CompareStr(const char* DataName, const char* Property){

		int PropL = (int)strlen(Property);
		int DataL = (int)strlen(DataName);

		if(PropL!=DataL) return false;

		for(int i=0; i<PropL; i++){

			if(DataName[i]=='\0' && Property[i]=='\0') return true;
			if(tolower(DataName[i])!=Property[i]) return false;

		}

		return true;
	}

	bool ReadBool(const char* Value){

		if(Value[1]=='\0'){

			if(Value[0]=='0') return false; 
			if(Value[0]=='1') return true; 

			if(Value[0]=='n') return false; 
			if(Value[0]=='y') return true; 

			return false;
		}
	
		if(CompareStr("true",Value))  return true;
		if(CompareStr("false",Value)) return false;

		if(CompareStr("yes",Value)) return true;
		if(CompareStr("no",Value))  return false;

		if(CompareStr("on",Value))  return true;
		if(CompareStr("off",Value)) return false;

		return false;
	}

#ifdef MPI_INCLUDED		//Defined if mpi.h is included

	void DistributeInputFileMPI(InputValueStructEx* Values, int nValues, size_t BaseAddr){

		for(int i=0; i<nValues; i++){

			char* Ptr = BaseAddr + (char*)Values[i].Ptr;
			int Size = Values[i].DataSz * Values[i].Count;

			MPI_Bcast(Ptr, Size, MPI_BYTE, _SrcThread, MPIComm);

			if(Values[i].CountRead != 0){

				char* CountPtr = BaseAddr + (char*)Values[i].CountRead;
				
				MPI_Bcast(CountPtr, 1, MPI_INT, _SrcThread, MPIComm);

			}
		}

	}

#endif

};



int ReadProperty(ifstream& File, char* Property, char* Value);
bool CompareStr(const char* DataName, const char* Property);
bool ReadBool(const char* Value);

#ifdef MPI_INCLUDED		//Defined if mpi.h is included
void DistributeInputFileMPI(InputValueStruct* Values, int nValues, int SrcThread, MPI_Comm Comm);
#endif

int ReadInputFile(const char* InputFileName, InputValueStruct* Values, int nValues){

	if(nValues <= 0)
		return 0;

	ifstream InFile;
	InFile.open(InputFileName);
	InFile.unsetf(ios::skipws);

	if(InFile.fail())
		return 1;			//File couldn't be opened

	bool* InDataSet = new bool[nValues];	//Whether each value has been set
	memset(InDataSet, 0, nValues);

	char Property[256];
	char Value[1024];

	while(true){

		int r = ReadProperty(InFile, Property, Value);

		if(r==0 || r==1){	//Valid input line

			for(int i=0; i<nValues; i++){

				if(!CompareStr(Values[i].DataName, Property)) continue;

				switch(Values[i].DataType){

					case DataType_Short:			//Interpret as short
						*((short*)(Values[i].Ptr)) = (short)atoi(Value);
						break;
					case DataType_Bool:				//Interpret as bool
						*((bool*)(Values[i].Ptr)) = ReadBool(Value);
						break;
					case DataType_Int:				//Interpret as int
						*((int*)(Values[i].Ptr)) = atoi(Value);
						break;
					case DataType_LongLongInt:		//Interpret as long long
						*((long long*)(Values[i].Ptr)) = atoll(Value);
						break;
					case DataType_Float:
						*((float*)(Values[i].Ptr)) = (float)atof(Value);
						break;
					case DataType_Double:			//Interpret as double
						*((double*)(Values[i].Ptr)) = atof(Value);
						break;
					case DataType_String:			//Copy over string
						memcpy(Values[i].Ptr, Value, min(Values[i].PtrLength,(int)sizeof(Value)));
						((char*)(Values[i].Ptr))[Values[i].PtrLength-1] = '\0';
						break;
				}
				
				InDataSet[i] = true;
				
			}
			

		}

		if(r==1 || r==3)	//End of file
			break;
		
	}

	int nRequired = 0;

	for(int i=0; i<nValues; i++){

		if(Values[i].Required && !InDataSet[i])		//Required value not set
			nRequired++;
		
	}

	if(nRequired!=0){

		cout << "Input file error: the following value";

		if(nRequired>1)
			cout << "s were ";
		else
			cout << " was ";
		
		cout << "left unspecified." << endl;

		for(int i=0; i<nValues; i++){

			if(Values[i].Required && !InDataSet[i])		//Required value not set
				cout << '\t' << Values[i].DataName << endl;

		}

		delete[] InDataSet;
		return 2;
	}

	delete[] InDataSet;
	return 0;		//Complete read-in
}

int ReadProperty(ifstream& File, char* Property, char* Value){

	int PropertyL = 0;
	int ValueL = 0;

	char c;
	int Stage = 0;					//Reading through property (0), value (1), or comment (2); Invalid line (-1)
	bool InQuotes = false;			//Reading a value in quotes (keep white spaces)
	while(true){

		File >> c;

		if(File.eof()){				//Reached end of file
			if(Stage==0 || Stage==-1)
				return 3;			//Reached end of file with no value
			
			Property[PropertyL] = '\0';
			Value[ValueL] = '\0';
			return 1;				//Reached end of file 
		}

		if(c=='\r' || c=='\t') continue; 
		if(c=='\n'){
			if(Stage==0 || Stage==-1)
				return 2;			//Reached end of invalid line with no value
			
			Property[PropertyL] = '\0';
			Value[ValueL] = '\0';
			return 0;				//Finished 
		}
		if(c=='/' && (Stage==0||Stage==1)){
			c = ' ';
			File >> c;
			if(c=='/'){				//Reached comment
				if(Stage == 0)		//Invalid line
					Stage = -1;
				else				//End of value
					Stage = 2;
			}else{
				c = '/';
				File.unget();
			}
		}

		if(Stage == 0){				//Reading in property name

			if(c==' ') continue;	//Ignore white-space

			if(c=='='){
				Stage = 1;			//Move onto value
			}else{
				Property[PropertyL] = tolower(c);
				PropertyL++;
			}

		}else if(Stage == 1){		//Reading in value

			if(!InQuotes && c==' ') continue;	//Ignore white-space outside quotes

			if(c=='"'){
				if(!InQuotes){
					InQuotes = true;
				}else{
					InQuotes = false;
					Stage = 2;
				}
			}else{
				Value[ValueL] = c;
				ValueL++;
			}

		}

	}

}

bool CompareStr(const char* DataName, const char* Property){

	for(int i=0; true; i++){

		if(DataName[i]=='\0' && Property[i]=='\0') return true;
		if(tolower(DataName[i])!=Property[i]) return false;

	}

}

bool ReadBool(const char* Value){

	if(Value[1]=='\0'){

		if(Value[0]=='0') return false; 
		if(Value[0]=='1') return true; 

		if(Value[0]=='n') return false; 
		if(Value[0]=='y') return true; 

		return false;
	}
	
	if(CompareStr("true",Value))  return true;
	if(CompareStr("false",Value)) return false;

	if(CompareStr("yes",Value)) return true;
	if(CompareStr("no",Value))  return false;

	if(CompareStr("on",Value))  return true;
	if(CompareStr("off",Value)) return false;

	return false;
}

#ifdef MPI_INCLUDED		//Defined if mpi.h is included

void DistributeInputFileMPI(InputValueStruct* Values, int nValues, int SrcThread, MPI_Comm Comm){

	for(int i=0; i<nValues; i++){
		MPI_Bcast(Values[i].Ptr, Values[i].PtrLength, MPI_BYTE, SrcThread, Comm);
	}

}

#endif