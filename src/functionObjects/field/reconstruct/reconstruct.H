/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::functionObjects::reconstruct

Description
    Calculates the reconstruction of a field; e.g., to construct a cell-centred
    velocity, U, from the face-centred flux, phi. The operation is limited to
    scalar and vector surface fields, and the output is a volume vector or
    tensor field.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    reconstruct.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_reconstruct_H
#define functionObjects_reconstruct_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class reconstruct Declaration
\*---------------------------------------------------------------------------*/

class reconstruct
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the magnitude of the field and register the result
        template<class Type>
        inline bool calcReconstruction();

        //- Calculate the reconstruction field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("reconstruct");


    // Constructors

        //- Construct from Time and dictionary
        reconstruct
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reconstruct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
