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
    Foam::filmFvPatch

Description
    Foam::filmFvPatch

SourceFiles
    filmFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef filmFvPatch_H
#define filmFvPatch_H

#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class filmFvPatch Declaration
\*---------------------------------------------------------------------------*/

class filmFvPatch
:
    public fvPatch
{
    // Private member data

        mutable scalarField deltaCoeffs_;


public:

    //- Runtime type information
    TypeName("film");


    // Constructors

        //- Construct from components
        filmFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm)
        {}


    // Member Functions

        // Access functions for demand driven data

            //- Return the face - cell distance coefficient
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre distance coefficient is returned
            virtual const scalarField& deltaCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
