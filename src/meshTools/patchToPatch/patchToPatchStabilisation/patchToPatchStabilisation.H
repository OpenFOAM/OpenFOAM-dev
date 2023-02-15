/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::patchToPatchStabilisation

Description
    Stabilisation data and routines for patch-to-patch interpolations

SourceFiles
    patchToPatchStabilisation.C
    patchToPatchStabilisationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchStabilisation_H
#define patchToPatchStabilisation_H

#include "patchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class patchToPatchStabilisation Declaration
\*---------------------------------------------------------------------------*/

class patchToPatchStabilisation
{
    // Private Data

        //- Is stabilisation occurring?
        bool stabilisation_;

        //- For each face, the local face used to stabilise the interpolation
        labelList localStabilisationCells_;

        //- Map from faces to local stabilisation faces
        autoPtr<distributionMap> stabilisationMapPtr_;


public:

    //- Run-time type information
    TypeName("patchToPatchStabilisation");


    // Constructors

        //- Construct null
        patchToPatchStabilisation();


    //- Destructor
    virtual ~patchToPatchStabilisation();


    // Member Functions

        //- Compute the stabilisation addressing if necessary
        void update
        (
            const polyPatch& patch,
            const PackedBoolList& faceCoupleds
        );

        //- Stabilise the given field if necessary
        template<class Type>
        void stabilise(Field<Type>& fld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchToPatchStabilisationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
