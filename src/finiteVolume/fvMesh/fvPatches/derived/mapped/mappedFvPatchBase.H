/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::mappedFvPatchBase

Description
    Base class for fv patches that provide interpolative mapping between two
    globally conforming fv patches

SourceFiles
    mappedFvPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFvPatchBase_H
#define mappedFvPatchBase_H

#include "mappedFvPatchBaseBase.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class mappedFvPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedFvPatchBase
:
    public mappedFvPatchBaseBase
{
    // Private Data

        //- Poly patch mapper base
        const mappedPatchBase& mapper_;


    // Private Member Functions

        //- Map/interpolate a neighbour patch field to this patch
        template<class Type>
        tmp<Field<Type>> fromNeighbour(const Field<Type>& nbrFld) const;

        //- Map/interpolate a patch field to the neighbour patch
        template<class Type>
        tmp<Field<Type>> toNeighbour(const Field<Type>& fld) const;


public:

    //- Runtime type information
    TypeName("mappedFvPatchBase");


    // Constructors

        //- Construct from a patch
        mappedFvPatchBase(const fvPatch& patch);


    //- Destructor
    virtual ~mappedFvPatchBase();


    // Member Functions

        //- Map/interpolate
        FOR_ALL_FIELD_TYPES
        (
            DEFINE_MAPPED_FV_PATCH_BASE_FROM_AND_TO_NEIGHBOUR,
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
