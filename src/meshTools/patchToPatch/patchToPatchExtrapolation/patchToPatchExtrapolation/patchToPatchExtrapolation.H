/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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
    Foam::patchToPatchExtrapolation

Description
    Base class for patch-to-patch extrapolations. Extrapolations fill in values
    for target faces that were not successfully interpolated because they fall
    outside the bounds of the source patch.

SourceFiles
    patchToPatchExtrapolation.C

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchExtrapolation_H
#define patchToPatchExtrapolation_H

#include "PackedBoolList.H"
#include "polyPatch.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class patchToPatchExtrapolation Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchExtrapolation
{
protected:

    // Protected Data

        //- List of uncoupled faces
        labelList uncoupledFaces_;

        //- Is extrapolation occurring?
        bool extrapolation_;


    // Protected Member Functions

        //- Update the extrapolation flag
        virtual void update(const PackedBoolList& faceCoupleds);


public:

    //- Run-time type information
    TypeName("patchToPatchExtrapolation");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        patchToPatchExtrapolation,
        word,
        (),
        ()
    );

    // Constructors

        //- Default construct
        patchToPatchExtrapolation();


    //- Destructor
    virtual ~patchToPatchExtrapolation();


    // Selector

        //- Select from name
        static autoPtr<patchToPatchExtrapolation> New(const word& type);


    // Member Functions

        //- Compute the extrapolation addressing if necessary
        virtual void update
        (
            const polyPatch& patch,
            const PackedBoolList& faceCoupleds
        ) = 0;

        //- Extrapolate the given field if necessary
        #define declareExtrapolateType(Type, nullArg) \
            virtual void extrapolate(Field<Type>& fld) const = 0;
        FOR_ALL_FIELD_TYPES(declareExtrapolateType);
        #undef declareExtrapolateType
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
