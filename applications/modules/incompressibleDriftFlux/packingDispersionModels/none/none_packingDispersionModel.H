/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::packingDispersionModels::none

Description
    Dummy packing dispersion model

Usage
    Example usage:
    \verbatim
        packingDispersionModel none;
    \endverbatim

SourceFiles
    none.C

\*---------------------------------------------------------------------------*/

#ifndef none_packingDispersionModel_H
#define none_packingDispersionModel_H

#include "packingDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace packingDispersionModels
{

/*---------------------------------------------------------------------------*\
                           Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public packingDispersionModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from relativeVelocity
        none
        (
            const relativeVelocityModel& relativeVelocity
        );

        //- Construct from components
        none
        (
            const dictionary& dict,
            const relativeVelocityModel& relativeVelocity
        );


    //- Destructor
    ~none();


    // Member Functions

        //- Return the derivative of the packing stress w.r.t. phase-fraction
        virtual tmp<volScalarField> sigmaPrime() const;

        //- Return the diffusion coefficient of the dispersed phase
        virtual tmp<volScalarField> Dd() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace packingDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
