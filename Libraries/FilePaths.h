#pragma once

#include <string>
#include <iostream>

using namespace std;

class FilePath{

public:
	enum DataType{
			Type_bool,
			Type_char,
			Type_int32,
			Type_uint32,
			Type_int64,
			Type_uint64,
			Type_float,
			Type_double,
			Type_string
		};

	struct ValueSub{
		const char* Text;
		const void* DataPtr;
		DataType Type;
		int nDecimalPlaces;
		bool floatExp;
	};

	enum SubOption{
		Sub_File,
		Sub_Folder,
		Sub_File_Folder
	};

private:
	const static int StrSz = 1024;
	char Folder[StrSz];
	char File[StrSz];
	char Buffer[StrSz];

	void _SetFolder(const char* folder){

		int l = (int)strlen(folder);
		int folderL = min(l, StrSz-2);

		if(l == 0){
			Folder[0] = '\0';
			return;
		}

		memcpy(Folder, folder, folderL);

		char fs = '/';					//Add / or \ at end if nesc.

		for(int i=0; i<folderL; i++){
			if(Folder[i] == '\\'){
				fs = '\\';
				break;
			}
		}

		if(Folder[folderL-1] != '\\' && Folder[folderL-1] != '/'){

			Folder[folderL++] = fs;

		}

		Folder[folderL] = '\0';
		
	}

	void _SetFile(const char* file){

		int l = (int)strlen(file);

		memcpy(File, file, min(l,StrSz-2));
		File[min(l,StrSz-2)] = '\0';

	}

	static bool ValToStr(char* Str, int StrSz, int* nBytes, ValueSub* Sub){

		int n = 1;

		const void* Data = Sub->DataPtr;
		DataType TypeIn = Sub->Type;

		switch(TypeIn){

			case Type_bool:
				Str[0] = *((bool*)Data) ? '1' : '0';
				break;
			case Type_char:
				Str[0] = *((char*)Data);
				break;
			case Type_int32:
				n = _snprintf(Str, StrSz, "%i", *((int*)Data));
				break;
			case Type_uint32:
				n = _snprintf(Str, StrSz, "%u", *((unsigned int*)Data));
				break;
			case Type_int64:
				n = _snprintf(Str, StrSz, "%lli", *((long long*)Data));
				break;
			case Type_uint64:
				n = _snprintf(Str, StrSz, "%llu", *((unsigned long long*)Data));
				break;
			case Type_float:
				n = _snprintf(Str, StrSz, (Sub->floatExp ? "%.*e" : "%.*f") , Sub->nDecimalPlaces, *((float*)Data));
				break;
			case Type_double:
				n = _snprintf(Str, StrSz, (Sub->floatExp ? "%.*e" : "%.*f") , Sub->nDecimalPlaces, *((double*)Data));
				break;
			case Type_string:
				n = _snprintf(Str, StrSz, "%s", (char*)Data);
				break;

		}

		*nBytes = n;

		return true;
	}

public:

	FilePath(){

		_SetFolder("");

		_SetFile("");

	}

	FilePath(const char* file){

		_SetFolder("");

		_SetFile(file);

	}

	FilePath(const char* file, ValueSub* Subs, int nSubs){

		_SetFolder("");
		
		ReplaceStrings(Buffer, StrSz, file, Subs, nSubs);

		_SetFile(Buffer);

	}

	FilePath(const char* folder, const char* file){

		_SetFolder(folder);

		_SetFile(file);

	}

	FilePath(const char* folder, const char* file, ValueSub* Subs, int nSubs, SubOption Option){

		if(Option == Sub_File_Folder || Option ==  Sub_Folder){

			ReplaceStrings(Buffer, StrSz, folder, Subs, nSubs);

			_SetFolder(Buffer);

		}else{

			_SetFolder(folder);

		}

		if(Option == Sub_File_Folder || Option ==  Sub_File){

			ReplaceStrings(Buffer, StrSz, file, Subs, nSubs);

			_SetFile(Buffer);

		}else{

			_SetFile(file);

		}

	}

	void SetFile(const char* file){

		_SetFile(file);

	}

	void SetFile(const char* file, ValueSub* Subs, int nSubs){
		
		ReplaceStrings(Buffer, StrSz, file, Subs, nSubs);

		_SetFile(Buffer);

	}

	void SetFolder(const char* folder){

		_SetFolder(folder);

	}

	void SetFolder(const char* folder, ValueSub* Subs, int nSubs){
		
		ReplaceStrings(Buffer, StrSz, folder, Subs, nSubs);

		_SetFolder(Buffer);

	}

	int ObtainPath(char* Dest, int DestSz){

		if(File[0] == '\\' || File[0] == '/' || File[0] == '.' || File[1] == ':'){		//File specifies path

			int fl = min( DestSz-1 , (int)strlen(File) );

			if(fl < 1)
				return 0;

			memcpy(Dest, File, fl);
			Dest[fl] = '\0';

			return fl;
		}

		int fileL = (int)strlen(File);
		int foldL = (int)strlen(Folder);

		int l = min( DestSz-1 , fileL+foldL );

		int nFolder = min( DestSz-1, foldL );
		int nFile   = min( DestSz-1-nFolder, fileL );

		if(nFolder > 0){
			memcpy(Dest, Folder, nFolder);
		}

		if(nFile > 0){
			memcpy(&Dest[nFolder], File, nFile);
		}

		Dest[l] = '\0';

		return l;
	}

	static void ReplaceStrings(char* OutStr, int OutStrSz, const char* filename, ValueSub* Subs, int nSubs){

		char** DataStrings = new char*[nSubs];		//Strings representing each data value

		for(int i=0; i<nSubs; i++){

			DataStrings[i] = new char[256];

			int nChar;

			ValToStr(DataStrings[i], 256, &nChar, &(Subs[i]));

			DataStrings[i][nChar] = '\0';
		}

		//Replace in file name

		int filenameL = (int)strlen(filename);

		int OutStrInd = 0;
		bool ExceedStr = false;

		for(int i=0; i<filenameL; i++){		//Each character in filename

			bool Replaced = false;

			for(int c=0; c<nSubs; c++){		//Each substitution value

				bool Replace = true;		

				int RepStrL = (int)strlen(Subs[c].Text);		//Length of substitution text

				if(i+RepStrL > filenameL)
					continue;

				for(int i2=0; i2<RepStrL; i2++){	//Find matching text
					if(tolower(filename[i+i2]) != tolower(Subs[c].Text[i2])){
						Replace = false;
						break;
					}
				}

				if(Replace){
					Replaced = true;

					int l = (int)strlen(DataStrings[c]);

					if(OutStrInd+l >= OutStrSz){
						ExceedStr = true;
						break;
					}

					memcpy(&OutStr[OutStrInd], DataStrings[c], l);

					OutStrInd += l;
					i += (RepStrL-1);

					break;
				}

			}

			if(!Replaced){
				OutStr[OutStrInd++] = filename[i];
			}

			if(ExceedStr || OutStrInd>=(OutStrSz-1))
				break;

		}

		OutStr[OutStrInd] = '\0';

		//Finalise

		for(int i=0; i<nSubs; i++)
			delete[] DataStrings[i];

		delete[] DataStrings;

	}

};