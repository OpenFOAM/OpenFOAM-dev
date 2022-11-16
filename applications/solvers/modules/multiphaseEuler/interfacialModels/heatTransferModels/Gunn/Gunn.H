/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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
    Foam::heatTransferModels::Gunn

Description
    Gunn (1978) correlation for fixed and fluidised beds with Re < 10^5
    and continuous phase fraction between 0.35 and 1.

    Reference:
    \verbatim
        Gunn, D. J. (1978).
        Transfer of heat or mass to particles in fixed and fluidised beds.
        International Journal of Heat and Mass Transfer, 21(4), 467-476.
    \endverbatim

SourceFiles
    Gunn.C

\*---------------------------------------------------------------------------*/

#ifndef Gunn_H
#define Gunn_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                            Class Gunn Declaration
\*---------------------------------------------------------------------------*/

class Gunn
:
    public heatTransferModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;


public:

    //- Runtime type information
    TypeName("Gunn");


    // Constructors

        //- Construct from a dictionary and an interface
        Gunn
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~Gunn();


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
