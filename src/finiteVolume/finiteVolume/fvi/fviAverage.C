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

#include "fviAverage.H"
#include "fviSurfaceIntegrate.H"
#include "fvMesh.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvi
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolInternalField<Type>>
average(const SurfaceField<Type>& ssf)
{
    const fvMesh& mesh = ssf.mesh()();

    tmp<VolInternalField<Type>> taverage
    (
        VolInternalField<Type>::New
        (
            "average("+ssf.name()+')',
            mesh,
            dimensioned<Type>("0", ssf.dimensions(), Zero)
        )
    );

    if (!mesh.nGeometricD())
    {
        return taverage;
    }

    VolInternalField<Type>& av = taverage.ref();

    av.primitiveFieldRef() =
    (
        surfaceSum(mesh.magSf()*ssf)().primitiveField()
       /surfaceSum(mesh.magSf())().primitiveField()
    );

    return taverage;
}


template<class Type>
tmp<VolInternalField<Type>>
average(const tmp<SurfaceField<Type>>& tssf)
{
    tmp<VolInternalField<Type>> taverage
    (
        fvi::average(tssf())
    );
    tssf.clear();
    return taverage;
}


template<class Type>
tmp<VolInternalField<Type>>
average(const VolField<Type>& vtf)
{
    return fvi::average(linearInterpolate(vtf));
}


template<class Type>
tmp<VolInternalField<Type>>
average(const tmp<VolField<Type>>& tvtf)
{
    tmp<VolInternalField<Type>> taverage
    (
        fvi::average(tvtf())
    );
    tvtf.clear();
    return taverage;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
