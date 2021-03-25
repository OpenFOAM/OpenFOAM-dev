/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    Foam::phaseCompressible::momentumTransportModel

Typedef
    Foam::phaseCompressible::laminarModel

Typedef
    Foam::phaseCompressible::RASModel

Typedef
    Foam::phaseCompressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for phase-compressible flow
    based on the phase-fluidThermo transport package.

SourceFiles
    phaseKinematicMomentumTransportModel.C
    phaseKinematicMomentumTransportModels.C

\*---------------------------------------------------------------------------*/

#ifndef phaseKinematicMomentumTransportModel_H
#define phaseKinematicMomentumTransportModel_H

#include "PhaseIncompressibleMomentumTransportModel.H"
#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "kinematicTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace phaseIncompressible
    {
        typedef PhaseIncompressibleMomentumTransportModel
            <kinematicTransportModel>
            momentumTransportModel;
        typedef laminarModel<momentumTransportModel> laminarModel;
        typedef RASModel<momentumTransportModel> RASModel;
        typedef LESModel<momentumTransportModel> LESModel;

        template<class BasicPhaseIncompressibleMomentumTransportModel>
        autoPtr<BasicPhaseIncompressibleMomentumTransportModel> New
        (
            const volScalarField& alpha,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const typename BasicPhaseIncompressibleMomentumTransportModel::
                transportModel& transport
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseKinematicMomentumTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
