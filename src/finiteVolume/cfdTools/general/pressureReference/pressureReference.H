/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::pressureReference

Description
    Provides controls for the pressure reference in closed-volume simulations

SourceFiles
    pressureReferenceI.H
    pressureReference.C

\*---------------------------------------------------------------------------*/

#ifndef pressureReference_H
#define pressureReference_H

#include "dimensionedScalar.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class pressureReference Declaration
\*---------------------------------------------------------------------------*/

class pressureReference
{
    // Private Data

        //- Optional cell in which the reference pressure is set
        label refCell_;

        //- Optional pressure reference level
        scalar refValue_;


public:

    // Constructors

        //- Construct from the simple/pimple sub-dictionary
        //  supporting the pressure and derived pressure fields (e.g. p_rgh)
        pressureReference
        (
            const volScalarField& p,
            const volScalarField& pRef,
            const dictionary& dict,
            const bool pRefRequired = true
        );

        //- Construct from the simple/pimple sub-dictionary
        pressureReference
        (
            const volScalarField& p,
            const dictionary& dict,
            const bool pRefRequired = true
        );


    // Member Functions

        //- Return the cell in which the reference pressure is set
        inline label refCell() const;

        //- Return the pressure reference level
        inline scalar refValue() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pressureReferenceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
