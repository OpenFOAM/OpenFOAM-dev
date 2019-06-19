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
    Foam::ensightStream

Description
    Abstract base class for writing Ensight data

SourceFiles
    ensightStream.C

\*---------------------------------------------------------------------------*/

#ifndef ensightStream_H
#define ensightStream_H

#include "fileName.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class ensightStream Declaration
\*---------------------------------------------------------------------------*/

class ensightStream
{
    // Private Data

        const fileName name_;


public:

    // Constructors

        //- Construct from components
        ensightStream(const fileName& f)
        :
            name_(f)
        {}

        //- Disallow default bitwise copy construction
        ensightStream(const ensightStream&) = delete;


    //- Destructor
    virtual ~ensightStream()
    {}


    // Member Functions

        const fileName& name() const
        {
            return name_;
        }

        virtual bool ascii() const = 0;

        virtual void write(const char*) = 0;

        virtual void write(const int) = 0;

        virtual void write(const scalarField&) = 0;

        virtual void write(const List<int>&) = 0;

        virtual void writePartHeader(const label) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ensightStream&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
