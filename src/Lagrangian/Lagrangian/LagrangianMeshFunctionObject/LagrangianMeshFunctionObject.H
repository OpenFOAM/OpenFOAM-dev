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
    Foam::functionObjects::LagrangianMeshFunctionObject

Description
    Base class for function objects relating to a Lagrangian mesh

SourceFiles
    LagrangianMeshFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_LagrangianMeshFunctionObject_H
#define functionObjects_LagrangianMeshFunctionObject_H

#include "objectRegistryFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class LagrangianMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class LagrangianMeshFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class LagrangianMeshFunctionObject
:
    public objectRegistryFunctionObject
{
private:

    // Private Data

        //- Reference to the Lagrangian mesh
        const LagrangianMesh& mesh_;


public:

    //- Runtime type information
    TypeName("LagrangianMeshFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        LagrangianMeshFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const word& LagrangianKeyword = "Lagrangian"
        );

        //- Disallow default bitwise copy construction
        LagrangianMeshFunctionObject
        (
            const LagrangianMeshFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~LagrangianMeshFunctionObject();


    // Member Functions

        //- Access the mesh
        inline const LagrangianMesh& mesh() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianMeshFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LagrangianMeshFunctionObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
