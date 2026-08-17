/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "fviCurl.H"
#include "fviGrad.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvi
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolInternalField<Type>>
curl(const VolField<Type>& vf)
{
    word nameCurlVf = "curl(" + vf.name() + ')';

    // Gausses theorem curl
    // tmp<VolInternalField<Type>> tcurlVf =
    // fvi::surfaceIntegrateExtrapolate(vf.mesh().Sf() ^ fvi::interpolate(vf));

    // Calculate curl as the Hodge dual of the skew-symmetric part of grad
    tmp<VolInternalField<Type>> tcurlVf =
        2.0*(*skew(fvi::grad(vf, nameCurlVf)));

    tcurlVf.ref().rename(nameCurlVf);

    return tcurlVf;
}


template<class Type>
tmp<VolInternalField<Type>>
curl(const tmp<VolField<Type>>& tvf)
{
    tmp<VolInternalField<Type>> Curl(fvi::curl(tvf()));
    tvf.clear();
    return Curl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
