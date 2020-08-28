/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::fileFormats::extendedFeatureEdgeMeshFormat

Description
    Provide a means of reading extendedFeatureEdgeMesh as featureEdgeMesh

SourceFiles
    extendedFeatureEdgeMeshFormat.C

\*---------------------------------------------------------------------------*/

#ifndef extendedFeatureEdgeMeshFormat_H
#define extendedFeatureEdgeMeshFormat_H

#include "edgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class extendedFeatureEdgeMeshFormat Declaration
\*---------------------------------------------------------------------------*/

class extendedFeatureEdgeMeshFormat
:
    public edgeMesh
{
public:

    // Constructors

        //- Construct from file name
        extendedFeatureEdgeMeshFormat(const fileName&);

        //- Disallow default bitwise copy construction
        extendedFeatureEdgeMeshFormat
        (
            const extendedFeatureEdgeMeshFormat&
        ) = delete;


    // Selectors

        //- Read file and return surface
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>
            (
                new extendedFeatureEdgeMeshFormat(name)
            );
        }


    //- Destructor
    virtual ~extendedFeatureEdgeMeshFormat()
    {}


    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extendedFeatureEdgeMeshFormat&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
