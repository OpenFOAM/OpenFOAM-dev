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
    Foam::heatTransferModels::Prandtl

Description
    This simple model creates a heat transfer coefficient in proportion with
    the corresponding drag model's momentum transfer coefficient. A
    user-defined Prandtl number and a harmonic average of the phases'
    specific heats are used to specify the constant of proportionality.

    This model has no physical basis. It exists primarily for testing
    purposes. It has the advantage of being applicable to any interface,
    including those representing segregated configurations.

    Example usage:
    \verbatim
    heatTransfer
    {
        gas_segregatedWith_liquid
        {
            type            Prandtl;
            Pr              0.7;
        }
    }
    \endverbatim

SourceFiles
    Prandtl.C

\*---------------------------------------------------------------------------*/

#ifndef Prandtl_H
#define Prandtl_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                            Class Prandtl Declaration
\*---------------------------------------------------------------------------*/

class Prandtl
:
    public heatTransferModel
{
    // Private Data

        //- Interface pointer
        const autoPtr<phaseInterface> interfacePtr_;

        //- Interface
        const phaseInterface& interface_;

        //- Prandtl number
        const dimensionedScalar Pr_;


public:

    //- Runtime type information
    TypeName("Prandtl");


    // Constructors

        //- Construct from a dictionary and an interface
        Prandtl
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~Prandtl();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
