/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fvcVolumeIntegrate.H"
#include "fvMesh.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>>
volumeIntegrate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return vf.mesh().V()*vf.primitiveField();
}


template<class Type>
tmp<Field<Type>>
volumeIntegrate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    tmp<Field<Type>> tvivf = tvf().mesh().V()*tvf().primitiveField();
    tvf.clear();
    return tvivf;
}


template<class Type>
tmp<Field<Type>> volumeIntegrate(const DimensionedField<Type, volMesh>& df)
{
    return df.mesh().V()*df.field();
}


template<class Type>
tmp<Field<Type>>
volumeIntegrate(const tmp<DimensionedField<Type, volMesh>>& tdf)
{
    tmp<Field<Type>> tdidf = tdf().mesh().V()*tdf().field();
    tdf.clear();
    return tdidf;
}


template<class Type>
dimensioned<Type>
domainIntegrate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return dimensioned<Type>
    (
        "domainIntegrate(" + vf.name() + ')',
        dimVol*vf.dimensions(),
        gSum(fvc::volumeIntegrate(vf))
    );
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    dimensioned<Type> integral = domainIntegrate(tvf());
    tvf.clear();
    return integral;
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const DimensionedField<Type, volMesh>& df
)
{
    return dimensioned<Type>
    (
        "domainIntegrate(" + df.name() + ')',
        dimVol*df.dimensions(),
        gSum(fvc::volumeIntegrate(df))
    );
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const tmp<DimensionedField<Type, volMesh>>& tdf
)
{
    dimensioned<Type> integral = domainIntegrate(tdf());
    tdf.clear();
    return integral;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
