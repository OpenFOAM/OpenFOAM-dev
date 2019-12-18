/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::cellQuality

Description
    Class calculates cell quality measures.

SourceFiles
    cellQuality.C

\*---------------------------------------------------------------------------*/

#ifndef cellQuality_H
#define cellQuality_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class cellQuality Declaration
\*---------------------------------------------------------------------------*/

class cellQuality
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;


public:

    // Constructors

        //- Construct from mesh
        cellQuality(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        cellQuality(const cellQuality&) = delete;


    // Member Functions

        //- Return cell non-orthogonality
        tmp<scalarField> nonOrthogonality() const;

        //- Return cell skewness
        tmp<scalarField> skewness() const;

        //- Return face non-orthogonality
        tmp<scalarField> faceNonOrthogonality() const;

        //- Return face skewness
        tmp<scalarField> faceSkewness() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellQuality&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
