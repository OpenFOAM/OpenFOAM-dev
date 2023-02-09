/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
    Foam::diffusiveMassTransferModels::sphericalDiffusiveMassTransfer

Description
    Model which applies an analytical solution for mass transfer from the
    surface of a sphere to the fluid within the sphere.

SourceFiles
    sphericalDiffusiveMassTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalDiffusiveMassTransfer_H
#define sphericalDiffusiveMassTransfer_H

#include "diffusiveMassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusiveMassTransferModels
{

/*---------------------------------------------------------------------------*\
               Class sphericalDiffusiveMassTransfer Declaration
\*---------------------------------------------------------------------------*/

class sphericalDiffusiveMassTransfer
:
    public diffusiveMassTransferModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;


public:

    //- Runtime type information
    TypeName("spherical");


    // Constructors

        //- Construct from a dictionary and an interface
        sphericalDiffusiveMassTransfer
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~sphericalDiffusiveMassTransfer();


    // Member Functions

        //- The implicit mass transfer coefficient
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusiveMassTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
