/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::symmetryPlaneFvPatch

Description
    Symmetry-plane patch.

SourceFiles
    symmetryPlaneFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryPlaneFvPatch_H
#define symmetryPlaneFvPatch_H

#include "fvPatch.H"
#include "symmetryPlanePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class symmetryPlaneFvPatch Declaration
\*---------------------------------------------------------------------------*/

class symmetryPlaneFvPatch
:
    public fvPatch
{
    // Private data

        //- Local reference cast into the symmetryPlane patch
        const symmetryPlanePolyPatch& symmetryPlanePolyPatch_;


public:

    //- Runtime type information
    TypeName(symmetryPlanePolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        symmetryPlaneFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            symmetryPlanePolyPatch_
            (
                refCast<const symmetryPlanePolyPatch>(patch)
            )
        {}


    // Member Functions

        //- Return symmetry plane normal
        const vector& n() const
        {
            return symmetryPlanePolyPatch_.n();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
