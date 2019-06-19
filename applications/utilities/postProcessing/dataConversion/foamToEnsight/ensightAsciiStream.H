/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::ensightAsciiStream

Description

SourceFiles
    ensightAsciiStream.C

\*---------------------------------------------------------------------------*/

#ifndef ensightAsciiStream_H
#define ensightAsciiStream_H

#include "ensightStream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightAsciiStream Declaration
\*---------------------------------------------------------------------------*/

class ensightAsciiStream
:
    public ensightStream
{
    // Private Data

        //- Description of data_
        OFstream str_;


public:

    // Constructors

        //- Construct from components
        ensightAsciiStream(const fileName& f, const Time& runTime)
        :
            ensightStream(f),
            str_
            (
                f,
                runTime.writeFormat(),
                runTime.writeVersion(),
                IOstream::UNCOMPRESSED
            )
        {

            str_.setf(ios_base::scientific, ios_base::floatfield);
            str_.precision(5);
        }

        //- Disallow default bitwise copy construction
        ensightAsciiStream(const ensightAsciiStream&) = delete;


    //- Destructor
    virtual ~ensightAsciiStream()
    {}


    // Member Functions

        virtual bool ascii() const
        {
            return true;
        }

        virtual void write(const char* c)
        {
            str_ << c << nl;
        }

        virtual void write(const int v)
        {
            str_ << setw(10) << v << nl;
        }

        virtual void write(const scalarField& sf)
        {
            forAll(sf, i)
            {
                if (mag(sf[i]) >= scalar(floatScalarVSmall))
                {
                    str_ << setw(12) << sf[i] << nl;
                }
                else
                {
                    str_ << setw(12) << scalar(0) << nl;
                }
            }
        }

        virtual void write(const List<int>& sf)
        {
            forAll(sf, i)
            {
                str_ << setw(10) << sf[i];
            }
            str_<< nl;
        }

        virtual void writePartHeader(const label partI)
        {
            str_<< "part" << nl
                << setw(10) << partI << nl;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ensightAsciiStream&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
