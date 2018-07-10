/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fileFormats::NASedgeFormat

Description
    Nastran edge reader.

SourceFiles
    NASedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef NASedgeFormat_H
#define NASedgeFormat_H

#include "edgeMesh.H"
#include "NASCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class NASedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class NASedgeFormat
:
    public edgeMesh,
    public NASCore
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        NASedgeFormat(const NASedgeFormat&);

        //- Disallow default bitwise assignment
        void operator=(const NASedgeFormat&);

public:

    // Constructors

        //- Construct from file name
        NASedgeFormat(const fileName&);


    // Selectors

        //- Read file and return edge mesh
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>
            (
                new NASedgeFormat(name)
            );
        }


    //- Destructor
    virtual ~NASedgeFormat()
    {}


    // Member Functions

        //- Read from a file
        virtual bool read(const fileName&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
