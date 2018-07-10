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
    Foam::nonuniformTransformCyclicPointPatch

Description
    Cyclic patch with slip constraint

SourceFiles
    nonuniformTransformCyclicPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonuniformTransformCyclicPointPatch_H
#define nonuniformTransformCyclicPointPatch_H

#include "cyclicPointPatch.H"
#include "nonuniformTransformCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class nonuniformTransformCyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class nonuniformTransformCyclicPointPatch
:
    public cyclicPointPatch
{

public:

    //- Runtime type information
    TypeName(nonuniformTransformCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        nonuniformTransformCyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        )
        :
            cyclicPointPatch(patch, bm)
        {}


    // Destructor

        virtual ~nonuniformTransformCyclicPointPatch()
        {}


    // Member Functions

        //- Return point unit normals.
        virtual const vectorField& pointNormals() const;

        //- Accumulate the effect of constraint direction of this patch
        virtual void applyConstraint
        (
            const label pointi,
            pointConstraint&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
