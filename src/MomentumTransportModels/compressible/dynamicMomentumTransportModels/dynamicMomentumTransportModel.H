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
    Foam::compressible::momentumTransportModel

Typedef
    Foam::compressible::laminarModel

Typedef
    Foam::compressible::RASModel

Typedef
    Foam::compressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for compressible flow
    based on the standard laminar transport package.

SourceFiles
    dynamicMomentumTransportModel.C
    dynamicMomentumTransportModels.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMomentumTransportModel_H
#define dynamicMomentumTransportModel_H

#include "CompressibleMomentumTransportModel.H"
#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "dynamicTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressible
    {
        typedef CompressibleMomentumTransportModel<dynamicTransportModel>
            momentumTransportModel;
        typedef laminarModel<momentumTransportModel> laminarModel;
        typedef RASModel<momentumTransportModel> RASModel;
        typedef LESModel<momentumTransportModel> LESModel;

        template<class BasicCompressibleMomentumTransportModel>
        autoPtr<BasicCompressibleMomentumTransportModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const typename BasicCompressibleMomentumTransportModel::
                transportModel& transport
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicMomentumTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
