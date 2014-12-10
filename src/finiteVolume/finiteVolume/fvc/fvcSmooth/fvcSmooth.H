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

InNamespace
    Foam::fvc

Description
    Provides functions smooth spread and sweep which use the FaceCellWave
    algorithm to smooth and redistribute the first field argument.

    smooth: smooths the field by ensuring the values in neighbouring cells are
            at least coeff* the cell value.

    spread: redistributes the field by spreading the maximum value within the
            region defined by the value (being between alphaMax and alphaMin)
            and gradient of alpha (where the difference between the values in
            neighbouring cells is larger than alphaDiff).

    sweep: redistributes the field by sweeping the maximum value where the
           gradient of alpha is large (where the difference between the values
           in neighbouring cells is larger than alphaDiff) away from that
           starting point of the sweep.

SourceFiles
    fvcSmooth.C

\*---------------------------------------------------------------------------*/

#ifndef fvcSmooth_H
#define fvcSmooth_H

#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvc
{
    void smooth
    (
        volScalarField& field,
        const scalar coeff
    );

    void spread
    (
        volScalarField& field,
        const volScalarField& alpha,
        const label nLayers,
        const scalar alphaDiff = 0.2,
        const scalar alphaMax = 0.99,
        const scalar alphaMin = 0.01
    );

    void sweep
    (
        volScalarField& field,
        const volScalarField& alpha,
        const label nLayers,
        const scalar alphaDiff = 0.2
    );
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
