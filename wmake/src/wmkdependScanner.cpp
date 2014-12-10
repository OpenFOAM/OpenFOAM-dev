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


#include <sstream>

// io.h and fcntl are used to ensure binary read from streams on windows
#ifdef _WIN32
# include <io.h>
# include <fcntl.h>
#endif

#include <climits>

#include "wmkdependScanner.h"

// values for the file stream buffering
#define MIN_BUFFER_LENGTH 1024        // 1KB
#define MAX_BUFFER_LENGTH (64*MIN_BUFFER_LENGTH)   // 64KB
// value for the heap management
#define HEAP_BLOCK_SIZE   (64*1024)   // 64KB


namespace wmake {

// * * * * * * * * * * * Miscellaneous String Routines * * * * * * * * * * * //

//
// string handling, byte character
//

std::string coco_stdString(const wchar_t* str)
{
    return str ? coco_stdString(str, 0, wcslen(str)) : std::string();
}


std::string coco_stdString(const wchar_t* str, unsigned length)
{
    return coco_stdString(str, 0, length);
}


std::string coco_stdString(const wchar_t* str, unsigned index, unsigned length)
{
    const unsigned len = (str && *str) ? length : 0;
    std::string dst;
    dst.reserve(len);

    for (unsigned i = 0; i < len; ++i)
    {
        dst += char(str[index+i] & 0xFF);
    }

    return dst;
}


std::string coco_stdStringUTF8(const wchar_t* str)
{
    return str ? coco_stdStringUTF8(str, 0, wcslen(str)) : std::string();
}


std::string coco_stdStringUTF8(const wchar_t* str, unsigned length)
{
    return coco_stdStringUTF8(str, 0, length);
}


std::string coco_stdStringUTF8(const wchar_t* str, unsigned index, unsigned length)
{
    const unsigned len = (str && *str) ? length : 0;
    std::string dst;
    dst.reserve(len);

    for (unsigned i = 0; i < len; ++i)
    {
        wchar_t wc = str[index+i];

        if (!(wc & ~0x0000007F))
        {
            // 0x00000000 - 0x0000007F [min. 8bit storage, 1-byte encoding)
            // 0aaaaaaa
            dst += char(wc);
        }
        else if (!(wc & ~0x000007FF))
        {
            // 0x00000080 - 0x000007FF [min. 16bit storage, 2-byte encoding]
            // 110bbbaa 10aaaaaa
            dst += char(0xC0 | ((wc >> 6) & 0x1F));
            dst += char(0x80 | ((wc) & 0x3F));
        }
        else if (!(wc & ~0x0000FFFF))
        {
            // 0x00000800 - 0x0000FFFF [min. 16bit storage, 3-byte encoding]
            // 1110bbbb 10bbbbaa 10aaaaaa
            dst += char(0xE0 | ((wc >> 12) & 0x0F));
            dst += char(0x80 | ((wc >> 6) & 0x3F));
            dst += char(0x80 | ((wc) & 0x3F));
        }
        else if (!(wc & ~0x001FFFFF))
        {
            // 0x00010000 - 0x001FFFFF [min. 24bit storage, 4-byte encoding]
            // 11110ccc 10ccbbbb 10bbbbaa 10aaaaaa
            dst += char(0xF0 | ((wc >> 18) & 0x07));
            dst += char(0x80 | ((wc >> 12) & 0x3F));
            dst += char(0x80 | ((wc >> 6) & 0x3F));
            dst += char(0x80 | ((wc) & 0x3F));
        }
//
// Not (yet) used - wchar_t storage is limited to 16bit on windows
// This also corresponds to the unicode BMP (Basic Multilingual Plane)
//
//        else if (!(wc & ~0x03FFFFFF))
//        {
//            // 0x00200000 - 0x03FFFFFF [min. 32bit storage, 5-byte encoding]
//            // 111110dd 10cccccc 10ccbbbb 10bbbbaa 10aaaaaa
//            dst += char(0xF8 | ((wc >> 24) & 0x03));
//            dst += char(0x80 | ((wc >> 18) & 0x3F));
//            dst += char(0x80 | ((wc >> 12) & 0x3F));
//            dst += char(0x80 | ((wc >> 6) & 0x3F));
//            dst += char(0x80 | ((wc) & 0x3F));
//        }
//        else if (!(wc & ~0x7FFFFFFF))
//        {
//            // 0x04000000 - 0x7FFFFFFF [min. 32bit storage, 6-byte encoding]
//            // 1111110d 10dddddd 10cccccc 10ccbbbb 10bbbbaa 10aaaaaa
//            dst += char(0xFC | ((wc >> 30) & 0x01));
//            dst += char(0x80 | ((wc >> 24) & 0x3F));
//            dst += char(0x80 | ((wc >> 18) & 0x3F));
//            dst += char(0x80 | ((wc >> 12) & 0x3F));
//            dst += char(0x80 | ((wc >> 6) & 0x3F));
//            dst += char(0x80 | ((wc) & 0x3F));
//        }
//
        else
        {
            // report anything unknown/invalid as replacement character U+FFFD
            dst += char(0xEF);
            dst += char(0xBF);
            dst += char(0xBD);
        }
    }

    return dst;
}


// * * * * * * * * * * * *  End of String Routines * * * * * * * * * * * * * //


Token::Token(wchar_t* value)
:
    kind(0),
    pos(0),
    col(0),
    line(0),
    val(value),
    next(NULL)
{}


Token::~Token()
{}


int Token::length() const
{
    return val ? wcslen(val) : 0;
}


// ----------------------------------------------------------------------------
// Buffer Implementation
// ----------------------------------------------------------------------------

Buffer::Buffer(Buffer* b)
:
	buf(b->buf),
	bufCapacity(b->bufCapacity),
	bufLen(b->bufLen),
	bufPos(b->bufPos),
	bufStart(b->bufStart),
	fileLen(b->fileLen),
	cStream(b->cStream),
	stdStream(b->stdStream),
	isUserStream_(b->isUserStream_)
{
	// avoid accidental deletion on any of these members
	b->buf = NULL;
	b->cStream = NULL;
	b->stdStream = NULL;
}


Buffer::Buffer(const char* chars, int len)
:
	buf(new unsigned char[len]),
	bufCapacity(len),
	bufLen(len),
	bufPos(0),
	bufStart(0),
	fileLen(len),
	cStream(NULL),
	stdStream(NULL),
	isUserStream_(false)
{
	memcpy(this->buf, chars, len*sizeof(char));
}


Buffer::Buffer(const unsigned char* chars, int len)
:
	buf(new unsigned char[len]),
	bufCapacity(len),
	bufLen(len),
	bufPos(0),
	bufStart(0),
	fileLen(len),
	cStream(NULL),
	stdStream(NULL),
	isUserStream_(false)
{
	memcpy(this->buf, chars, len*sizeof(char));
}


Buffer::Buffer(FILE* ifh, bool isUserStream)
:
	buf(NULL),
	bufCapacity(0),
	bufLen(0),
	bufPos(0),
	bufStart(0),
	fileLen(0),
	cStream(ifh),
	stdStream(NULL),
	isUserStream_(isUserStream)
{
// ensure binary read on windows
#ifdef _WIN32
	_setmode(_fileno(cStream), _O_BINARY);
#endif

	if (CanSeek())
	{
		fseek(cStream, 0, SEEK_END);
		fileLen = ftell(cStream);
		fseek(cStream, 0, SEEK_SET);
		bufLen = (fileLen < MAX_BUFFER_LENGTH) ? fileLen : MAX_BUFFER_LENGTH;
		bufStart = INT_MAX; // nothing in the buffer so far
	}

	bufCapacity = (bufLen > 0) ? bufLen : MIN_BUFFER_LENGTH;
	buf = new unsigned char[bufCapacity];
	if (fileLen > 0) SetPos(0);          // setup buffer to position 0 (start)
	else bufPos = 0; // index 0 is already after the file, thus Pos = 0 is invalid
	if (bufLen == fileLen && CanSeek()) Close();
}


Buffer::Buffer(std::istream* istr, bool isUserStream)
:
	buf(NULL),
	bufCapacity(0),
	bufLen(0),
	bufPos(0),
	bufStart(0),
	fileLen(0),
	cStream(NULL),
	stdStream(istr),
	isUserStream_(isUserStream)
{
#if _WIN32
	// TODO: ensure binary read on windows?
#endif
}


Buffer::~Buffer()
{
	Close();
	if (buf)
	{
		delete[] buf;
		buf = NULL;
	}
}


void Buffer::Close()
{
	if (!isUserStream_)
	{
		if (cStream)
		{
			fclose(cStream);
			cStream = NULL;
		}
		if (stdStream)
		{
			delete stdStream;
			stdStream = 0;
		}
	}
}


int Buffer::Read()
{
	if (stdStream)
	{
		int ch = stdStream->get();
		if (stdStream->eof())
		{
			return EoF;
		}
		return ch;
	}

	if (bufPos < bufLen) {
		return buf[bufPos++];
	}
	else if (GetPos() < fileLen) {
		SetPos(GetPos()); // shift buffer start to Pos
		return buf[bufPos++];
	}
	else if (cStream && !CanSeek() && (ReadNextStreamChunk() > 0)) {
		return buf[bufPos++];
	}

	return EoF;
}

bool Buffer::isUTF8() const
{
	return false;
}

int UTF8Buffer::Read()
{
	int ch;
	do {
		ch = Buffer::Read();
		// until we find a utf8 start (0xxxxxxx or 11xxxxxx)
	} while (ch != EoF && ch >= 128 && ((ch & 0xC0) != 0xC0));
	if (ch < 128 || ch == EoF) {
		// nothing to do, first 127 chars are identical in ASCII and UTF8
		// 0xxxxxxx or end of file character
	}
	else if ((ch & 0xF0) == 0xF0) {
		// 0x00010000 - 0x001FFFFF [min. 24bit storage, 4-byte encoding]
		// 11110ccc 10ccbbbb 10bbbbaa 10aaaaaa
		// CAUTION: this should probably be disallowed since it overflows
		// wchar_t on windows and overflows the max (0xFFFF) used here
		int c1 = ch & 0x07; ch = Buffer::Read();
		int c2 = ch & 0x3F; ch = Buffer::Read();
		int c3 = ch & 0x3F; ch = Buffer::Read();
		int c4 = ch & 0x3F;
		ch = (((((c1 << 6) | c2) << 6) | c3) << 6) | c4;
	}
	else if ((ch & 0xE0) == 0xE0) {
		// 0x00000800 - 0x0000FFFF [min. 16bit storage, 3-byte encoding]
		// 1110bbbb 10bbbbaa 10aaaaaa
		int c1 = ch & 0x0F; ch = Buffer::Read();
		int c2 = ch & 0x3F; ch = Buffer::Read();
		int c3 = ch & 0x3F;
		ch = (((c1 << 6) | c2) << 6) | c3;
	}
	else if ((ch & 0xC0) == 0xC0) {
		// 0x00000080 - 0x000007FF [min. 16bit storage, 2-byte encoding]
		// 110bbbaa 10aaaaaa
		int c1 = ch & 0x1F; ch = Buffer::Read();
		int c2 = ch & 0x3F;
		ch = (c1 << 6) | c2;
	}
	return ch;
}


bool UTF8Buffer::isUTF8() const
{
	return true;
}


int Buffer::Peek()
{
	int curPos = GetPos();
	int ch = Read();
	SetPos(curPos);
	return ch;
}


int Buffer::GetPos() const
{
	if (stdStream)
	{
		return stdStream->tellg();
	}

	return bufPos + bufStart;
}


void Buffer::SetPos(int value)
{
	if (stdStream)
	{
		stdStream->seekg(value, std::ios::beg);
		return;
	}

	if ((value >= fileLen) && cStream && !CanSeek())
	{
		// Wanted position is after buffer and the stream
		// is not seek-able e.g. network or console,
		// thus we have to read the stream manually till
		// the wanted position is in sight.
		while ((value >= fileLen) && (ReadNextStreamChunk() > 0))
		{}
	}

	if ((value < 0) || (value > fileLen))
	{
		fwprintf(stderr, L"--- buffer out of bounds access, position: %d\n", value);
		::exit(1);
	}

	if ((value >= bufStart) && (value < (bufStart + bufLen))) // already in buffer
	{
		bufPos = value - bufStart;
	}
	else if (cStream) // must be swapped in
	{
		fseek(cStream, value, SEEK_SET);
		bufLen = fread(buf, sizeof(char), bufCapacity, cStream);
		bufStart = value; bufPos = 0;
	}
	else
	{
		bufPos = fileLen - bufStart; // make Pos return fileLen
	}
}


//
// Read the next chunk of bytes from the stream, increases the buffer
// if needed and updates the fields fileLen and bufLen.
// Returns the number of bytes read.
//
int Buffer::ReadNextStreamChunk()
{
	int freeLen = bufCapacity - bufLen;
	if (freeLen == 0)
	{
		// in the case of a growing input stream
		// we can neither seek in the stream, nor can we
		// foresee the maximum length, thus we must adapt
		// the buffer size on demand.
		bufCapacity = bufLen * 2;
		unsigned char *newBuf = new unsigned char[bufCapacity];
		memcpy(newBuf, buf, bufLen*sizeof(char));
		delete[] buf;
		buf = newBuf;
		freeLen = bufLen;
	}
	int read = fread(buf + bufLen, sizeof(char), freeLen, cStream);
	if (read > 0)
	{
		fileLen = bufLen = (bufLen + read);
		return read;
	}
	// end of stream reached
	return 0;
}


bool Buffer::CanSeek() const
{
	return cStream && (ftell(cStream) != -1);
}

// ----------------------------------------------------------------------------
// Scanner Implementation
// ----------------------------------------------------------------------------

Scanner::Scanner(const char* buf, int len)
:
	buffer(new Buffer(buf, len))
{
	Init();
}


Scanner::Scanner(const unsigned char* buf, int len)
:
	buffer(new Buffer(buf, len))
{
	Init();
}


Scanner::Scanner(FILE* ifh)
:
	buffer(new Buffer(ifh, true))
{
	Init();
}


#ifdef _WIN32
Scanner::Scanner(const std::wstring& fileName)
{
	FILE* ifh;

	if ((ifh = _wfopen(fileName.c_str(), L"rb")) == NULL)
	{
		fwprintf(stderr, L"--- Cannot open file %ls\n", fileName.c_str());
		::exit(1);
	}
	buffer = new Buffer(ifh, false);
	Init();
}
#endif


Scanner::Scanner(const std::string& fileName)
{
	FILE* ifh;
	if ((ifh = fopen(fileName.c_str(), "rb")) == NULL)
	{
		fwprintf(stderr, L"--- Cannot open file %s\n", fileName.c_str());
		::exit(1);
	}
	buffer = new Buffer(ifh, false);
	Init();
}


Scanner::Scanner(std::istream& istr)
:
	buffer(new Buffer(&istr, true))
{
	Init();
}


Scanner::~Scanner()
{
	char* cur = reinterpret_cast<char*>(firstHeap);

#ifdef COCO_DEBUG_HEAP
	fwprintf(stderr, L"~Scanner:\n");
#endif

	while (cur)
	{
		cur = *(reinterpret_cast<char**>(cur + HEAP_BLOCK_SIZE));
		free(firstHeap);
#ifdef COCO_DEBUG_HEAP
		fwprintf
		(
			stderr, L"    free %p -> %p\n",
			firstHeap,
			reinterpret_cast<char*>(firstHeap) + HEAP_BLOCK_SIZE
		);
#endif
		firstHeap = cur;
	}
	delete[] tval;
	delete buffer;
}


void Scanner::Init()
{
	for (int i = 36; i <= 36; ++i) start.set(i, 7);
	for (int i = 65; i <= 90; ++i) start.set(i, 7);
	for (int i = 95; i <= 95; ++i) start.set(i, 7);
	for (int i = 97; i <= 122; ++i) start.set(i, 7);
	start.set(34, 1);
	start.set(39, 4);
	start.set(35, 11);
	start.set(10, 12);
	start.set(59, 13);
	start.set(Buffer::EoF, -1);

	keywords.set(L"include", 6);
	keywords.set(L"import", 8);

	tvalLength = 128;
	tval = new wchar_t[tvalLength]; // text of current token
	tlen = 0;
	tval[tlen] = 0;

	// HEAP_BLOCK_SIZE byte heap + pointer to next heap block
	heap = malloc(HEAP_BLOCK_SIZE + sizeof(void*));
	firstHeap = heap;
	heapEnd =
		reinterpret_cast<void**>
		(reinterpret_cast<char*>(heap) + HEAP_BLOCK_SIZE);
	*heapEnd = 0;
	heapTop = heap;
	if (sizeof(Token) > HEAP_BLOCK_SIZE)
	{
		fwprintf(stderr, L"--- Too small HEAP_BLOCK_SIZE\n");
		::exit(1);
	}
#ifdef COCO_DEBUG_HEAP
	fwprintf
	(
		stderr, L"Scanner::init: firstHeap %p -> %p\n",
		firstHeap,
		reinterpret_cast<char*>(firstHeap) + HEAP_BLOCK_SIZE
	);
#endif

	pos = -1; line = 1; col = 0;
	oldEols = 0;
	NextCh();
	if (ch == 0xEF)   // check optional byte order mark for UTF-8
	{                 // Windows-specific magic
		NextCh(); int ch1 = ch;
		NextCh(); int ch2 = ch;
		if (ch1 != 0xBB || ch2 != 0xBF)
		{
			fwprintf(stderr, L"Illegal byte order mark at start of file");
			::exit(1);
		}
		Buffer *oldBuf = buffer;
		buffer = new UTF8Buffer(oldBuf); col = 0;
		delete oldBuf; oldBuf = NULL;
		NextCh();
	}
	else
	{
		// FORCE_UTF8 was defined
		// use UTF8Buffer without relying on a byte order mark.
		Buffer *oldBuf = buffer;
		buffer = new UTF8Buffer(oldBuf); col = 0;
		delete oldBuf; oldBuf = NULL;
	}

	pt = tokens = CreateToken(); // first token is a dummy
}


void Scanner::NextCh()
{
	if (oldEols > 0)
	{
		ch = EOL;
		oldEols--;
	}
	else
	{
		pos = buffer->GetPos();
		ch = buffer->Read(); col++;
		// replace isolated '\r' by '\n' in order to make
		// eol handling uniform across Windows, Unix and Mac
		if (ch == '\r' && buffer->Peek() != '\n') ch = EOL;
		if (ch == EOL) { line++; col = 0; }
	}
}


void Scanner::AddCh()
{
	if (tlen >= tvalLength)
	{
		tvalLength *= 2;
		wchar_t *newBuf = new wchar_t[tvalLength];
		memcpy(newBuf, tval, tlen*sizeof(wchar_t));
		delete[] tval;
		tval = newBuf;
	}
	if (ch != Buffer::EoF)
	{
		tval[tlen++] = ch;
		NextCh();
	}
}



bool Scanner::Comment0() {
	int level = 1, pos0 = pos, line0 = line, col0 = col;
	NextCh();
	if (ch == '/') {
		NextCh();
		while (true) {
			if (ch == 10) {
				level--;
				if (level == 0) { oldEols = line - line0; NextCh(); return true; }
				NextCh();
			} else if (ch == buffer->EoF) return false;
			else NextCh();
		}
	} else {
		buffer->SetPos(pos0); NextCh(); line = line0; col = col0;
	}
	return false;
}

bool Scanner::Comment1() {
	int level = 1, pos0 = pos, line0 = line, col0 = col;
	NextCh();
	if (ch == '*') {
		NextCh();
		while (true) {
			if (ch == '*') {
				NextCh();
				if (ch == '/') {
					level--;
					if (level == 0) { oldEols = line - line0; NextCh(); return true; }
					NextCh();
				}
			} else if (ch == '/') {
				NextCh();
				if (ch == '*') {
					level++; NextCh();
				}
			} else if (ch == buffer->EoF) return false;
			else NextCh();
		}
	} else {
		buffer->SetPos(pos0); NextCh(); line = line0; col = col0;
	}
	return false;
}

void Scanner::CreateHeapBlock()
{
	char* cur = reinterpret_cast<char*>(firstHeap);

#ifdef COCO_DEBUG_HEAP
	fwprintf(stderr, L"CreateHeapBlock: tokens %p\n", tokens);
#endif

	// release unused blocks
	while
	(
	    (reinterpret_cast<char*>(tokens) < cur)
	 || (reinterpret_cast<char*>(tokens) > (cur + HEAP_BLOCK_SIZE))
	)
	{
		cur = *(reinterpret_cast<char**>(cur + HEAP_BLOCK_SIZE));
#ifdef COCO_DEBUG_HEAP
		fwprintf
		(
			stderr, L"    free %p -> %p\n",
			firstHeap,
			reinterpret_cast<char*>(firstHeap) + HEAP_BLOCK_SIZE
		);
#endif
		free(firstHeap);
		firstHeap = cur;
	}

	// HEAP_BLOCK_SIZE byte heap + pointer to next heap block
	void* newHeap = malloc(HEAP_BLOCK_SIZE + sizeof(void*));
	*heapEnd = newHeap;
	heapEnd =
		reinterpret_cast<void**>
		(reinterpret_cast<char*>(newHeap) + HEAP_BLOCK_SIZE);
	*heapEnd = 0;
	heap = newHeap;
	heapTop = heap;
#ifdef COCO_DEBUG_HEAP
	fwprintf
	(
		stderr, L"    malloc %p -> %p\n",
		newHeap,
		reinterpret_cast<char*>(newHeap) + HEAP_BLOCK_SIZE
	);
#endif
}


Token* Scanner::CreateToken()
{
	const int reqMem = sizeof(Token);
	if
	(
	    (reinterpret_cast<char*>(heapTop) + reqMem)
	 >= reinterpret_cast<char*>(heapEnd)
	)
	{
		CreateHeapBlock();
	}
	// token 'occupies' heap starting at heapTop
	Token* tok = reinterpret_cast<Token*>(heapTop);
	// increment past this part of the heap, which is now used
	heapTop =
		reinterpret_cast<void*>
		(reinterpret_cast<char*>(heapTop) + reqMem);
	tok->val  = NULL;
	tok->next = NULL;
	return tok;
}


void Scanner::AppendVal(Token* tok)
{
	const int reqMem = (tlen + 1) * sizeof(wchar_t);
	if
	(
	    (reinterpret_cast<char*>(heapTop) + reqMem)
	 >= reinterpret_cast<char*>(heapEnd)
	)
	{
		if (reqMem > HEAP_BLOCK_SIZE)
		{
			fwprintf(stderr, L"--- Too long token value\n");
			::exit(1);
		}
		CreateHeapBlock();
	}

	// add text value from heap
	tok->val = reinterpret_cast<wchar_t*>(heapTop);

	// increment past this part of the heap, which is now used
	heapTop =
		reinterpret_cast<void*>
		(reinterpret_cast<char*>(heapTop) + reqMem);

	// copy the currently parsed tval into the token
	wcsncpy(tok->val, tval, tlen);
	tok->val[tlen] = '\0';
}


Token* Scanner::NextToken()
{
	while
	(
	    ch == ' '
	 || ch == 9
	) NextCh();
	if ((ch == '/' && Comment0()) || (ch == '/' && Comment1())) return NextToken();
	int recKind = noSym;
	int recEnd = pos;
	t = CreateToken();
	t->pos = pos; t->col = col; t->line = line;
	int state = start.state(ch);
	tlen = 0; AddCh();

	switch (state)
	{
		case -1: { t->kind = eofSym; break; } // NextCh already done
		case 0: {
			case_0:
			if (recKind != noSym) {
				tlen = recEnd - t->pos;
				SetScannerBehindT();
			}
			t->kind = recKind; break;
		} // NextCh already done
		case 1:
			case_1:
			if (ch <= 9 || (ch >= 11 && ch <= 12) || (ch >= 14 && ch <= '!') || (ch >= '#' && ch <= '[') || (ch >= ']' && ch <= 65535)) {AddCh(); goto case_1;}
			else if (ch == '"') {AddCh(); goto case_3;}
			else if (ch == 92) {AddCh(); goto case_2;}
			else {goto case_0;}
		case 2:
			case_2:
			if ((ch >= ' ' && ch <= '~')) {AddCh(); goto case_1;}
			else {goto case_0;}
		case 3:
			case_3:
			{t->kind = 1; break;}
		case 4:
			case_4:
			if (ch <= 9 || (ch >= 11 && ch <= 12) || (ch >= 14 && ch <= '!') || (ch >= '#' && ch <= '&') || (ch >= '(' && ch <= '[') || (ch >= ']' && ch <= 65535)) {AddCh(); goto case_4;}
			else if (ch == 39) {AddCh(); goto case_8;}
			else if (ch == 92) {AddCh(); goto case_5;}
			else {goto case_0;}
		case 5:
			case_5:
			if ((ch >= ' ' && ch <= '~')) {AddCh(); goto case_4;}
			else {goto case_0;}
		case 6:
			case_6:
			{t->kind = 4; break;}
		case 7:
			case_7:
			recEnd = pos; recKind = 3;
			if (ch == '$' || (ch >= '0' && ch <= '9') || (ch >= 'A' && ch <= 'Z') || ch == '_' || (ch >= 'a' && ch <= 'z')) {AddCh(); goto case_7;}
			else if (ch == '.') {AddCh(); goto case_9;}
			else {t->kind = 3; std::wstring literal(tval, tlen); t->kind = keywords.get(literal, t->kind); break;}
		case 8:
			case_8:
			recEnd = pos; recKind = 2;
			if (ch <= 9 || (ch >= 11 && ch <= 12) || (ch >= 14 && ch <= '!') || (ch >= '#' && ch <= '&') || (ch >= '(' && ch <= '[') || (ch >= ']' && ch <= 65535)) {AddCh(); goto case_4;}
			else if (ch == 39) {AddCh(); goto case_8;}
			else if (ch == 92) {AddCh(); goto case_5;}
			else {t->kind = 2; break;}
		case 9:
			case_9:
			if (ch == '$' || (ch >= 'A' && ch <= 'Z') || ch == '_' || (ch >= 'a' && ch <= 'z')) {AddCh(); goto case_10;}
			else if (ch == '*') {AddCh(); goto case_6;}
			else {goto case_0;}
		case 10:
			case_10:
			recEnd = pos; recKind = 3;
			if (ch == '$' || (ch >= '0' && ch <= '9') || (ch >= 'A' && ch <= 'Z') || ch == '_' || (ch >= 'a' && ch <= 'z')) {AddCh(); goto case_10;}
			else if (ch == '.') {AddCh(); goto case_9;}
			else {t->kind = 3; std::wstring literal(tval, tlen); t->kind = keywords.get(literal, t->kind); break;}
		case 11:
			{t->kind = 5; break;}
		case 12:
			{t->kind = 7; break;}
		case 13:
			{t->kind = 9; break;}
	}
	AppendVal(t);
	return t;
}


void Scanner::SetScannerBehindT()
{
	buffer->SetPos(t->pos);
	NextCh();
	line = t->line; col = t->col;
	for (int i = 0; i < tlen; i++) NextCh();
}


// get the next token (possibly a token already seen during peeking)
Token* Scanner::Scan()
{
	if (tokens->next == NULL) {
		pt = tokens = NextToken();
	}
	else {
		pt = tokens = tokens->next;
	}
	return tokens;
}


// peek for the next token, ignore pragmas
Token* Scanner::Peek()
{
	do
	{
		if (pt->next == NULL)
		{
			pt->next = NextToken();
		}
		pt = pt->next;
	} while (pt->kind > maxT);   // skip pragmas

	return pt;
}


// make sure that peeking starts at the current scan position
void Scanner::ResetPeek()
{
	pt = tokens;
}


int Scanner::Line() const
{
	return line;
}


void Scanner::Line(int lineNo)
{
	line = lineNo;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
