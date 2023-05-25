/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::heatTransferModels::timeScaleFiltered

Description
    A time scale filtering wrapper around an underlying heat transfer model
    intended for simulations involving dispersed phase with a very small
    diameter. The model limits the heat transfer coefficient based on the
    relaxation time given.

SourceFiles
    timeScaleFiltered.C

\*---------------------------------------------------------------------------*/

#ifndef timeScaleFilteredHeatTransfer_H
#define timeScaleFilteredHeatTransfer_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                Class timeScaleFilteredHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class timeScaleFilteredHeatTransfer
:
    public heatTransferModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;

        //- Pointer to the underlying heat transfer model
        autoPtr<heatTransferModel> heatTransferModel_;

        //- The minimum relaxation time
        const dimensionedScalar minRelaxTime_;


public:

    //- Runtime type information
    TypeName("timeScaleFiltered");


    // Constructors

        //- Construct from a dictionary and an interface
        timeScaleFilteredHeatTransfer
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~timeScaleFilteredHeatTransfer();


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
