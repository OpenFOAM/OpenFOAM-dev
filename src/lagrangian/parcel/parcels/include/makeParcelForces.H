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

\*---------------------------------------------------------------------------*/

#ifndef makeParcelForces_H
#define makeParcelForces_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SphereDragForce.H"
#include "NonSphereDragForce.H"
#include "SchillerNaumannDragForce.H"
#include "WenYuDragForce.H"
#include "ErgunWenYuDragForce.H"
#include "PlessisMasliyahDragForce.H"

#include "SaffmanMeiLiftForce.H"
#include "TomiyamaLiftForce.H"

#include "GravityForce.H"
#include "NonInertialFrameForce.H"
#include "ParamagneticForce.H"
#include "PressureGradientForce.H"
#include "SRFForce.H"
#include "VirtualMassForce.H"

#include "ScaledForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelForces(CloudType)                                            \
                                                                               \
    makeParticleForceModel(CloudType);                                         \
    makeParticleForceModelType(SphereDragForce, CloudType);                    \
    makeParticleForceModelType(NonSphereDragForce, CloudType);                 \
    makeParticleForceModelType(SchillerNaumannDragForce, CloudType);           \
    makeParticleForceModelType(WenYuDragForce, CloudType);                     \
    makeParticleForceModelType(ErgunWenYuDragForce, CloudType);                \
    makeParticleForceModelType(PlessisMasliyahDragForce, CloudType);           \
    makeParticleForceModelType(SaffmanMeiLiftForce, CloudType);                \
    makeParticleForceModelType(TomiyamaLiftForce, CloudType);                  \
    makeParticleForceModelType(GravityForce, CloudType);                       \
    makeParticleForceModelType(NonInertialFrameForce, CloudType);              \
    makeParticleForceModelType(ParamagneticForce, CloudType);                  \
    makeParticleForceModelType(PressureGradientForce, CloudType);              \
    makeParticleForceModelType(SRFForce, CloudType);                           \
    makeParticleForceModelType(VirtualMassForce, CloudType);                   \
    makeParticleForceModelType(ScaledForce, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
