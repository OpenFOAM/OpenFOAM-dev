/*---------------------------------------------------------------------------*\
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

Class
    Foam::IOstream

Description
    An IOstream is an abstract base class for all input/output systems; be
    they streams, files, token lists etc.

    The basic operations are construct, close, read token, read primitive
    and read binary block.  In addition version control and line number
    counting is incorporated.  Usually one would use the read primitive
    member functions, but if one were reading a stream on unknown data
    sequence one can read token by token, and then analyse.

SourceFiles
    IOstream.C

\*---------------------------------------------------------------------------*/

#ifndef IOstream_H
#define IOstream_H

#include "char.H"
#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "scalar.H"
#include "fileName.H"
#include "InfoProxy.H"

#include <iostream>

#if __GNUC__ < 3
#   define ios_base ios
#endif

using std::ios_base;
using std::istream;
using std::ostream;

using std::cin;
using std::cout;
using std::cerr;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOstream Declaration
\*---------------------------------------------------------------------------*/

class IOstream
{

public:

    // Public data types

        //- Enumeration for whether the stream open or closed
        enum streamAccess
        {
            OPENED,
            CLOSED
        };

        //- Enumeration for the format of data in the stream
        enum streamFormat
        {
            ASCII,
            BINARY
        };

        //- Ostream operator
        friend Ostream& operator<<(Ostream& os, const streamFormat& sf);

        //- Version number type
        class versionNumber
        {
            //- The version number
            scalar versionNumber_;

            //- The version number as an integer
            int index_;


        public:

            // Constructors

                //- Construct from number
                versionNumber(const scalar num)
                :
                    versionNumber_(num),
                    index_(numberToIndex(num))
                {}

                //- Construct from Istream
                versionNumber(Istream& is)
                :
                    versionNumber_(readScalar(is)),
                    index_(numberToIndex(versionNumber_))
                {}


            // Member functions

                //- Convert a version number into an index
                int numberToIndex(const scalar num) const
                {
                    return int(10*num + SMALL);
                }

                //- Return major version
                int majorVersion() const
                {
                    return int(versionNumber_);
                }

                //- Return minor version
                int minorVersion() const
                {
                    return int(10.0*(versionNumber_ - majorVersion()));
                }

                //- Return the versionNumber as a character string
                string str() const;


            // Member operators

                //- Are these versionNumbers the same?
                bool operator==(const versionNumber& vn)
                {
                    return index_ == vn.index_;
                }

                //- Are these versionNumbers different?
                bool operator!=(const versionNumber& vn)
                {
                    return index_ != vn.index_;
                }

                //- Is this version older than the one given
                bool operator<(const versionNumber& vn)
                {
                    return index_ < vn.index_;
                }

                //- Is this version the same as or older than the one given
                bool operator<=(const versionNumber& vn)
                {
                    return index_ <= vn.index_;
                }

                //- Is this version newer than the one given
                bool operator>(const versionNumber& vn)
                {
                    return index_ > vn.index_;
                }

                //- this version the same as or newer than the one given
                bool operator>=(const versionNumber& vn)
                {
                    return index_ >= vn.index_;
                }


            //- Ostream operator
            friend Ostream& operator<<(Ostream& os, const versionNumber& vn);
        };


        //- Enumeration for the format of data in the stream
        enum compressionType
        {
            UNCOMPRESSED,
            COMPRESSED
        };


    // Public static data

        //- Original version number
        static const versionNumber originalVersion;

        //- Current version number
        static const versionNumber currentVersion;

        //- Default precision
        static unsigned int precision_;


private:

    // Private data

        //- Name of the stream
        static fileName name_;

        streamFormat format_;
        versionNumber version_;
        compressionType compression_;

        streamAccess openClosed_;
        ios_base::iostate ioState_;


protected:

    // Protected data

        label lineNumber_;


    // Protected Member Functions

        // Access

            //- Set stream opened
            void setOpened()
            {
                openClosed_ = OPENED;
            }

            //- Set stream closed
            void setClosed()
            {
                openClosed_ = CLOSED;
            }

            //- Set stream state
            void setState(ios_base::iostate state)
            {
                ioState_ = state;
            }

            //- Set stream to be good
            void setGood()
            {
                ioState_ = ios_base::iostate(0);
            }


public:

    // Constructors

        //- Construct setting format and version
        IOstream
        (
            streamFormat format,
            versionNumber version,
            compressionType compression=UNCOMPRESSED
        )
        :
            format_(format),
            version_(version),
            compression_(compression),
            openClosed_(CLOSED),
            ioState_(ios_base::iostate(0)),
            lineNumber_(0)
        {
            setBad();
        }


    //- Destructor
    virtual ~IOstream()
    {}


    // Member functions

        // Access

            //- Return the name of the stream
            //  Useful for Fstream to return the filename
            virtual const fileName& name() const
            {
                return name_;
            }

            //- Return non-const access to the name of the stream
            //  Useful to alter the stream name
            virtual fileName& name()
            {
                return name_;
            }


        // Check

            //- Check IOstream status for given operation
            //  print IOstream state if error has occured
            virtual bool check(const char* operation) const;

            //- Check IOstream status for given operation
            //  print IOstream state if error has occured and exit
            void fatalCheck(const char* operation) const;

            //- Return true if stream has been opened
            bool opened() const
            {
                return openClosed_ == OPENED;
            }

            //- Return true if stream is closed
            bool closed() const
            {
                return openClosed_ == CLOSED;
            }

            //- Return true if next operation might succeed
            bool good() const
            {
                return ioState_ == 0;
            }

            //- Return true if end of input seen
            bool eof() const
            {
                return ioState_ & ios_base::eofbit;
            }

            //- Return true if next operation will fail
            bool fail() const
            {
                return ioState_ & (ios_base::badbit | ios_base::failbit);
            }

            //- Return true if stream is corrupted
            bool bad() const
            {
                return ioState_ & ios_base::badbit;
            }

            //- Return non-zero if the stream has not failed
            operator void*() const
            {
                return fail()
                    ? reinterpret_cast<void*>(0)
                    : reinterpret_cast<void*>(-1);
            }

            //- Return true if the stream has failed
            bool operator!() const
            {
                return fail();
            }


        // Stream state functions

            //- Return stream format of given format name
            static streamFormat formatEnum(const word&);

            //- Return current stream format
            streamFormat format() const
            {
                return format_;
            }

            //- Set the stream format
            streamFormat format(const streamFormat fmt)
            {
                streamFormat fmt0 = format_;
                format_ = fmt;
                return fmt0;
            }

            //- Set the stream format from word
            streamFormat format(const word& fmt)
            {
                streamFormat fmt0 = format_;
                format_ = formatEnum(fmt);
                return fmt0;
            }

            //- Return the stream version
            versionNumber version() const
            {
                return version_;
            }

            //- Set the stream version
            versionNumber version(const versionNumber ver)
            {
                versionNumber ver0 = version_;
                version_ = ver;
                return ver0;
            }

            //- Return compression of given compression name
            static compressionType compressionEnum(const word&);

            //- Return the stream compression
            compressionType compression() const
            {
                return compression_;
            }

            //- Set the stream compression
            compressionType compression(const compressionType cmp)
            {
                compressionType cmp0 = compression_;
                compression_ = cmp;
                return cmp0;
            }

            //- Set the stream compression from word
            compressionType compression(const word& cmp)
            {
                compressionType cmp0 = compression_;
                compression_ = compressionEnum(cmp);
                return cmp0;
            }

            //- Return current stream line number
            label lineNumber() const
            {
                return lineNumber_;
            }

            //- Return current stream line number
            label& lineNumber()
            {
                return lineNumber_;
            }

            //- Set the stream line number
            label lineNumber(const label ln)
            {
                label ln0 = lineNumber_;
                lineNumber_ = ln;
                return ln0;
            }

            //- Return flags of stream
            virtual ios_base::fmtflags flags() const = 0;

            //- Return the default precision
            static unsigned int defaultPrecision()
            {
                return precision_;
            }

            //- Reset the default precision (and return old precision)
            static unsigned int defaultPrecision(unsigned int p)
            {
                unsigned int precision0 = precision_;
                precision_ = p;
                return precision0;
            }

            //- Set stream to have reached eof
            void setEof()
            {
                ioState_ |= ios_base::eofbit;
            }

            //- Set stream to have failed
            void setFail()
            {
                ioState_ |= ios_base::failbit;
            }

            //- Set stream to be bad
            void setBad()
            {
                ioState_ |= ios_base::badbit;
            }

            //- Set flags of stream
            virtual ios_base::fmtflags flags(const ios_base::fmtflags f) = 0;

            //- Set flags of stream
            ios_base::fmtflags setf(const ios_base::fmtflags f)
            {
                return flags(flags() | f);
            }

            //- Set flags of given field of stream
            ios_base::fmtflags setf
            (
                const ios_base::fmtflags f,
                const ios_base::fmtflags mask
            )
            {
                return flags((flags() & ~mask) | (f & mask));
            }

            //- Unset flags of stream
            void unsetf(const ios_base::fmtflags uf)
            {
                flags(flags()&~uf);
            }


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream&) const;

            //- Check given stream state bits
            void print(Ostream&, const int streamState) const;


        // Info

            //- Return info proxy.
            //  Used to print IOstream information to a stream
            InfoProxy<IOstream> info() const
            {
                return *this;
            }
};


Ostream& operator<<(Ostream& os, const IOstream::streamFormat& sf);
Ostream& operator<<(Ostream& os, const IOstream::versionNumber& vn);


// --------------------------------------------------------------------
// ------ Manipulators (not taking arguments)
// --------------------------------------------------------------------

typedef IOstream& (*IOstreamManip)(IOstream&);

//- operator<< handling for manipulators without arguments
inline IOstream& operator<<(IOstream& io, IOstreamManip f)
{
    return f(io);
}


inline IOstream& dec(IOstream& io)
{
    io.setf(ios_base::dec, ios_base::dec|ios_base::hex|ios_base::oct);
    return io;
}

inline IOstream& hex(IOstream& io)
{
    io.setf(ios_base::hex, ios_base::dec|ios_base::hex|ios_base::oct);
    return io;
}

inline IOstream& oct(IOstream& io)
{
    io.setf(ios_base::oct, ios_base::dec|ios_base::hex|ios_base::oct);
    return io;
}

inline IOstream& fixed(IOstream& io)
{
    io.setf(ios_base::fixed, ios_base::floatfield);
    return io;
}

inline IOstream& scientific(IOstream& io)
{
    io.setf(ios_base::scientific, ios_base::floatfield);
    return io;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
