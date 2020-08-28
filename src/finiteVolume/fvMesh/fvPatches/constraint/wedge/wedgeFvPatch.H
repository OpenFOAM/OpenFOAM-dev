/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::wedgeFvPatch

Description
    Wedge front and back plane patch.

SourceFiles
    wedgeFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFvPatch_H
#define wedgeFvPatch_H

#include "fvPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class wedgeFvPatch Declaration
\*---------------------------------------------------------------------------*/

class wedgeFvPatch
:
    public fvPatch
{
    // Private Data

        const wedgePolyPatch& wedgePolyPatch_;


public:

    //- Runtime type information
    TypeName(wedgePolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        wedgeFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm);


    // Member Functions

        // Access

            //- Return face transformation tensor
            const tensor& faceT() const
            {
                return wedgePolyPatch_.faceT();
            }

            //- Return neighbour-cell transformation tensor
            const tensor& cellT() const
            {
                return wedgePolyPatch_.cellT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
