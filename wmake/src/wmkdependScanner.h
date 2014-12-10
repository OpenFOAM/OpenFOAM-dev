/*---------------------------------*- C++ -*---------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

@file wmkdependParser.atg

Description
    An attributed Coco/R grammar to parse C/C++, Fortran and Java files
    for include and import statements.

SourceFiles
    generated

\*---------------------------------------------------------------------------*/
// This file was generated with Coco/R C++ (10 Mar 2010)
// http://www.ssw.uni-linz.ac.at/coco/
// with these defines:
//     - FORCE_UTF8


#ifndef COCO_wmkdependSCANNER_H__
#define COCO_wmkdependSCANNER_H__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cwchar>
#include <string>
#include <fstream>
#include <iostream>

namespace wmake {

// * * * * * * * * * * * Miscellaneous String Routines * * * * * * * * * * * //

//! Simple lower-case string transformation
template<class StringT>
inline void coco_string_toLower(StringT& str)
{
	for
	(
		typename StringT::iterator iter = str.begin();
		iter != str.end();
		++iter
	)
	{
		if (*iter >= 'A' && *iter <= 'Z')
		{
			*iter += ('a' - 'A');   // lower-case
		}
	}
}


//! Simple string hashing function
template<class StringT>
inline int coco_string_hash(const StringT& str)
{
	int h = 0;
	for
	(
		typename StringT::const_iterator iter = str.begin();
		iter != str.end();
		++iter
	)
	{
		h = (h * 7) ^ *iter;
	}
	return h < 0 ? -h : h;
}


//
// String conversions
// ~~~~~~~~~~~~~~~~~~

//! Convert wide string to double
inline double coco_string_toDouble(const wchar_t* str)
{
	return str ? wcstod(str, NULL) : 0;
}


//! Convert wide string to long
inline long coco_string_toLong(const wchar_t* str)
{
	return str ? wcstol(str, NULL, 10) : 0;
}


//! A byte string (restricted to 8bit values) by copying str
std::string coco_stdString(const wchar_t* str);

//! A byte string (restricted to 8bit values) by copying str,
//! up to length characters long
std::string coco_stdString(const wchar_t* str, unsigned length);

//! A byte substring (restricted to 8bit values) of str,
//! starting at index and length characters long
std::string coco_stdString(const wchar_t* str, unsigned index, unsigned length);

//! A UTF8 byte string by copying str
std::string coco_stdStringUTF8(const wchar_t* str);

//! A UTF8 byte string by copying str, up to length characters long
std::string coco_stdStringUTF8(const wchar_t* str, unsigned length);

//! A UTF8 byte substring, starting at index and length characters long
std::string coco_stdStringUTF8(const wchar_t* str, unsigned index, unsigned length);

// * * * * * * * * * * * *  End of String Routines * * * * * * * * * * * * * //



/*---------------------------------------------------------------------------*\
                            Class Token Declaration
\*---------------------------------------------------------------------------*/
/*!
 * @brief Scanner Token
 *
 * @note since each Token is allocated by the internal heap mechanism,
 * the destructor does not clean up the val member.
 */
class Token
{
public:
	int kind;       //!< token kind
	int pos;        //!< token position in the source text (starting at 0)
	int col;        //!< token column (starting at 1)
	int line;       //!< token line (starting at 1)
	wchar_t* val;   //!< token value (normally allocated from the internal heap)
	Token *next;    //!< Peek tokens are kept in linked list

	int length() const;    //!< The length of val, or 0 if val is NULL

	//! Construct null Token, optionally with pointer to a string value
	Token(wchar_t* value = 0);
	~Token();       //!< Destructor - does not cleanup val member

	//! Token val as byte string (restricted to 8bit values)
	inline std::string toString() const
	{
		return coco_stdString(val);
	}

	//! Token val as byte string (restricted to 8bit values), up to length characters long
	inline std::string toString(unsigned length) const
	{
		return coco_stdString(val, length);
	}

	//! Token val as byte string (restricted to 8bit values), starting at index and length characters long
	inline std::string toString(unsigned index, unsigned length) const
	{
		return coco_stdString(val, index, length);
	}

	//! Token val as UTF8 byte string
	inline std::string toStringUTF8() const
	{
		return coco_stdStringUTF8(val);
	}

	//! Token val as UTF8 byte string, up to length characters long
	inline std::string toStringUTF8(unsigned length) const
	{
		return coco_stdStringUTF8(val, length);
	}

	//! Token val as UTF8 byte substring, starting at index and length characters long
	inline std::string toStringUTF8(unsigned index, unsigned length) const
	{
		return coco_stdStringUTF8(this->val, index, length);
	}

};


/*---------------------------------------------------------------------------*\
                           Class Buffer Declaration
\*---------------------------------------------------------------------------*/
/*!
 * @brief Scanner Buffer
 *
 * This Buffer supports the following cases:
 * -# seekable stream (file)
 *    -# whole stream in buffer
 *    -# part of stream in buffer
 * -# non seekable stream (network, console)
 */
class Buffer
{
	unsigned char *buf; //!< input buffer
	int bufCapacity;    //!< capacity of buf
	int bufLen;         //!< length of buffer
	int bufPos;         //!< current position in buffer
	int bufStart;       //!< position of first byte in buffer relative to input stream
	int fileLen;        //!< length of input stream (may change if the stream is no file)
	FILE* cStream;      //!< input stdio stream (normally seekable)
	std::istream* stdStream;  //!< STL std stream (seekable)
	bool isUserStream_;       //!< was the stream opened by the user?

	int ReadNextStreamChunk();
	bool CanSeek() const; //!< true if stream can be seeked otherwise false

protected:
	Buffer(Buffer*);    //!< for the UTF8Buffer

public:
	//! max unicode characters is 0xFFFF (16bit storage)
	static const int MaxChar = 65535;
	static const int EoF = MaxChar + 1;

	//! Copy buffer contents from constant character string
	Buffer(const char* chars, int len);

	//! Copy buffer contents from constant character string
	Buffer(const unsigned char* chars, int len);

	//! @brief Attach buffer to a stdio stream.
	//! User streams are not closed in the destructor
	Buffer(FILE*, bool isUserStream = true);

	//! @brief Attach buffer to an STL standard stream
	//! User streams are not closed in the destructor
	explicit Buffer(std::istream*, bool isUserStream = true);

	//! Close stream (but not user streams) and free buf (if any)
	virtual ~Buffer();

	virtual void Close();   //!< Close stream (but not user streams)
	virtual int Read();     //!< Get character from stream or buffer
	virtual int Peek();     //!< Peek character from stream or buffer

	virtual int GetPos() const;
	virtual void SetPos(int value);
	virtual bool isUTF8() const;  //!< Return false - buffer is not UTF8
};


/*---------------------------------------------------------------------------*\
                         Class UTF8Buffer Declaration
\*---------------------------------------------------------------------------*/
//! A Scanner Buffer variant that decodes UTF-8 characters into 16bit unicode
class UTF8Buffer : public Buffer
{
public:
	UTF8Buffer(Buffer* b) : Buffer(b) {}
	virtual int Read();
	virtual bool isUTF8() const;  //!< Return true - buffer is UTF8
};


/*---------------------------------------------------------------------------*\
                         Class StartStates Declaration
\*---------------------------------------------------------------------------*/
//! maps characters (integers) to start states of tokens as a HashTable
class StartStates
{
	//! HashTable entry
	struct Entry
	{
		int key;        //<! The lookup key
		int val;        //<! The data
		Entry *next;    //<! Pointer next Entry in sub-list

		Entry(int k, int v, Entry *n=0)
		:
			key(k), val(v), next(n)
		{}
	};

	static const int size_ = 128;   //<! fixed HashTable size
	Entry **table_;

public:
	StartStates()
	:
		table_(new Entry*[size_])
	{
		memset(table_, 0, size_*sizeof(Entry*));
	}

	virtual ~StartStates()
	{
		for (int i = 0; i < size_; ++i)
		{
			Entry *e = table_[i];
			while (e)
			{
				Entry *next = e->next;
				delete e;
				e = next;
			}
		}
		delete[] table_;
	}

	void set(int key, int val)
	{
		const int hashIndex = unsigned(key) % size_;
		table_[hashIndex] = new Entry(key, val, table_[hashIndex]);
	}

	int state(int key)
	{
		Entry *e = table_[unsigned(key) % size_];
		while (e && e->key != key) e = e->next;
		return e ? e->val : 0;
	}
};


/*---------------------------------------------------------------------------*\
                         Class KeywordMap Declaration
\*---------------------------------------------------------------------------*/
//! maps strings to integers (identifiers to keyword kinds) as a HashTable
class KeywordMap
{
	//! HashTable entry
	struct Entry
	{
		const std::wstring key;  //<! The lookup key
		int val;                 //<! The data
		Entry *next;             //<! Pointer next Entry in sub-list

		Entry(const std::wstring& k, int v, Entry *n=0)
		:
			key(k), val(v), next(n)
		{}
	};

	static const int size_ = 128;   //<! fixed HashTable size
	Entry **table_;

public:
	KeywordMap()
	:
		table_(new Entry*[size_])
	{
		memset(table_, 0, size_*sizeof(Entry*));
	}

	virtual ~KeywordMap()
	{
		for (int i = 0; i < size_; ++i)
		{
			Entry *e = table_[i];
			while (e)
			{
				Entry *next = e->next;
				delete e;
				e = next;
			}
		}
		delete[] table_;
	}

	void set(const std::wstring& key, int val)
	{
		const int hashIndex = coco_string_hash(key) % size_;
		table_[hashIndex] = new Entry(key, val, table_[hashIndex]);
	}

	int get(const std::wstring& key, int defaultVal)
	{
		Entry *e = table_[coco_string_hash(key) % size_];
		while (e && e->key != key) e = e->next;
		return e ? e->val : defaultVal;
	}
};


/*---------------------------------------------------------------------------*\
                           Class Scanner Declaration
\*---------------------------------------------------------------------------*/
//! A Coco/R Scanner
class Scanner
{
	static const int maxT = 10;
	static const int noSym = 10;
	static const int eofSym = 0;    //!< end-of-file token id
	static const char EOL = '\n';   //!< end-of-line character

	void *firstHeap;  //!< the start of the heap management
	void *heap;       //!< the currently active block
	void *heapTop;    //!< the top of the heap
	void **heapEnd;   //!< the end of the last heap block

	StartStates start;   //!< A map of start states for particular characters
	KeywordMap keywords; //!< A hash of keyword literals to token kind

	Token *t;         //!< current token
	wchar_t *tval;    //!< text of current token
	int tvalLength;   //!< maximum capacity (length) for tval
	int tlen;         //!< length of tval

	Token *tokens;    //!< list of tokens already peeked (first token is a dummy)
	Token *pt;        //!< current peek token

	int ch;           //!< current input character
	int pos;          //!< byte position of current character
	int line;         //!< line number of current character
	int col;          //!< column number of current character
	int oldEols;      //!< the number of EOLs that appeared in a comment

	void CreateHeapBlock();       //!< add a heap block, freeing unused ones
	Token* CreateToken();         //!< fit token on the heap
	void AppendVal(Token* tok);   //!< adjust tok->val to point to the heap and copy tval into it
	void SetScannerBehindT();

	void Init();      //!< complete the initialization for the constructors
	void NextCh();    //!< get the next input character into ch
	void AddCh();     //!< append the character ch to tval
	bool Comment0();
	bool Comment1();
	Token* NextToken();  //!< get the next token

public:
	//! The scanner buffer
	Buffer *buffer;

	//! Attach scanner to an existing character buffer
	Scanner(const char* chars, int len);

	//! Attach scanner to an existing character buffer
	Scanner(const unsigned char* chars, int len);

	//! Attach scanner to an existing open file handle
	Scanner(FILE*);

#ifdef _WIN32
	//! Open a file for reading and attach scanner - Windows only
	explicit Scanner(const std::wstring& fileName);
#endif

	//! Open a file for reading and attach scanner
	explicit Scanner(const std::string& fileName);

	//! Attach scanner to an existing open STL standard stream
	explicit Scanner(std::istream&);

	~Scanner();        //!< free heap and allocated memory
	Token* Scan();     //!< get the next token (possibly a token already seen during peeking)
	Token* Peek();     //!< peek for the next token, ignore pragmas
	void ResetPeek();  //!< ensure that peeking starts at the current scan position

	int  Line() const;     //!< Return the current line
	void Line(int lineNo); //!< Define the starting line for reporting errors

}; // end Scanner

} // End namespace

#endif // COCO_wmkdependSCANNER_H__

// ************************************************************************* //
