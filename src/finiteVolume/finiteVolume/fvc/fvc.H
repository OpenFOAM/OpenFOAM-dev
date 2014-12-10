/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Namespace
    Foam::fvc

Description
    Namespace of functions to calculate explicit derivatives.

\*---------------------------------------------------------------------------*/

#ifndef fvc_H
#define fvc_H

#include "fv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceInterpolate.H"
#include "fvcVolumeIntegrate.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcAverage.H"
#include "fvcReconstruct.H"
#include "fvcDdt.H"
#include "fvcDDt.H"
#include "fvcD2dt2.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcGrad.H"
#include "fvcMagSqrGradGrad.H"
#include "fvcSnGrad.H"
#include "fvcCurl.H"
#include "fvcLaplacian.H"
#include "fvcSup.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
