/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

Typedef
    Foam::compressible::turbulenceModel

Typedef
    Foam::compressible::RASModel

Typedef
    Foam::compressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for compressible flow
    based on the standard laminar transport package.

SourceFiles
    turbulentFluidThermoModels.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentFluidThermoModel_H
#define turbulentFluidThermoModel_H

#include "SpecificCompressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "EddyDiffusivity.H"
#include "RASModel.H"
#include "LESModel.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressible
    {
        typedef ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermo> >
            turbulenceModel;

        typedef SpecificCompressibleTurbulenceModel
        <
            RASModel<EddyDiffusivity<turbulenceModel> >
        > RASModel;

        typedef SpecificCompressibleTurbulenceModel
        <
            LESModel<EddyDiffusivity<turbulenceModel> >
        > LESModel;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
