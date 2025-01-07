/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionObjects::cloudFvMeshFunctionObject

Description
    Base class for function objects that refer to an fvMesh and a cloud. Used,
    for example, to construct and write accumulated cloud properties as finite
    volume fields.

SourceFiles
    cloudFvMeshFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudFvMeshFunctionObject_H
#define functionObjects_cloudFvMeshFunctionObject_H

#include "fvMeshFunctionObject.H"
#include "cloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class cloudFvMeshFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class cloudFvMeshFunctionObject
:
    public fvMeshFunctionObject,
    public cloudFunctionObject
{
public:

    //- Runtime type information
    TypeName("cloudFvMeshFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        cloudFvMeshFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudFvMeshFunctionObject(const cloudFvMeshFunctionObject&) = delete;


    //- Destructor
    virtual ~cloudFvMeshFunctionObject();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudFvMeshFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
