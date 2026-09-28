/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    An unimplemented extrapolation procedure for patch-to-patch interpolations.
    This simply checks that the interpolation is complete and that
    extrapolation is not necessary and raises an error if that is not the case.

SourceFiles
    patchToPatchExtrapolation.C

\*---------------------------------------------------------------------------*/

#ifndef nonePatchToPatchExtrapolation_H
#define nonePatchToPatchExtrapolation_H

#include "patchToPatchExtrapolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatchExtrapolations
{

/*---------------------------------------------------------------------------*\
                            Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public patchToPatchExtrapolation
{
public:

    //- Run-time type information
    TypeName("none");


    // Constructors

        //- Construct null
        none();


    //- Destructor
    virtual ~none();


    // Member Functions

        //- Compute the extrapolation addressing if necessary
        void update
        (
            const polyPatch& patch,
            const PackedBoolList& faceCoupleds
        );

        //- Extrapolate the given field if necessary
        #define declareExtrapolateType(Type, nullArg)                          \
            virtual void extrapolate(Field<Type>& fld) const                   \
            {                                                                  \
                NotImplemented;                                                \
            }
        FOR_ALL_FIELD_TYPES(declareExtrapolateType);
        #undef declareExtrapolateType
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchExtrapolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
