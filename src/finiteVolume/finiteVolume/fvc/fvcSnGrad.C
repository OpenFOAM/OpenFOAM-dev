/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fvcSnGrad.H"
#include "fvMesh.H"
#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<SurfaceField<Type>>
snGrad
(
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::snGradScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().snGrad(name)
    )().snGrad(vf);
}


template<class Type>
tmp<SurfaceField<Type>>
snGrad
(
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<SurfaceField<Type>> SnGrad
    (
        fvc::snGrad(tvf(), name)
    );
    tvf.clear();
    return SnGrad;
}


template<class Type>
tmp<SurfaceField<Type>>
snGrad
(
    const VolField<Type>& vf
)
{
    return fvc::snGrad(vf, "snGrad(" + vf.name() + ')');
}


template<class Type>
tmp<SurfaceField<Type>>
snGrad
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<SurfaceField<Type>> SnGrad
    (
        fvc::snGrad(tvf())
    );
    tvf.clear();
    return SnGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
