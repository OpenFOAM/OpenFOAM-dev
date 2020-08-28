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
    Foam::laminarFlameSpeedModels::constant

Description
    Constant laminar flame speed model.

SourceFiles
    selectLaminarFlameSpeed.C

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "laminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{

/*---------------------------------------------------------------------------*\
                          Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public laminarFlameSpeed
{
    // Private Data

        dimensionedScalar Su_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary and psiuReactionThermo
        constant
        (
            const dictionary&,
            const psiuReactionThermo&
        );

        //- Disallow default bitwise copy construction
        constant(const constant&) = delete;


    //- Destructor
    virtual ~constant();


    // Member Operators

        //- Return the laminar flame speed [m/s]
        tmp<volScalarField> operator()() const;

        //- Disallow default bitwise assignment
        void operator=(const constant&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarFlameSpeedModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
