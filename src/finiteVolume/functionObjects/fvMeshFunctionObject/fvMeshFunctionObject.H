/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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
    Foam::functionObjects::fvMeshFunctionObject

Description
    Specialisation of Foam::functionObject for an Foam::fvMesh, providing a
    reference to the Foam::fvMesh.

    If the selected region is not an Foam::fvMesh a Foam::FatalError will be
    generated.

See also
    Foam::functionObject
    Foam::functionObjects::objectRegistryFunctionObject
    Foam::functionObjects::regionFunctionObject

SourceFiles
    fvMeshFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fvMeshFunctionObject_H
#define functionObjects_fvMeshFunctionObject_H

#include "regionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class fvMeshFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class fvMeshFunctionObject
:
    public regionFunctionObject
{

protected:

    // Protected member data

        //- Reference to the fvMesh
        const fvMesh& mesh_;


    // Protected member functions

        //- Return the patch set corresponding to the patch selection entries
        //  in the given dictionary
        labelHashSet patchSet
        (
            const dictionary&,
            const bool optional = false
        ) const;


public:

    //- Runtime type information
    TypeName("fvMeshFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        fvMeshFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from the region objectRegistry
        fvMeshFunctionObject
        (
            const word& name,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        fvMeshFunctionObject(const fvMeshFunctionObject&) = delete;


    //- Destructor
    virtual ~fvMeshFunctionObject();


    // Member Functions

        //- Return a reference to the mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
