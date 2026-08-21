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

#include "fviVolumeIntegrate.H"
#include "fvMesh.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvi
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>> volumeIntegrate(const DimensionedField<Type, fvMesh>& df)
{
    return df.mesh().V().primitiveField()*df.primitiveField();
}


template<class Type>
tmp<Field<Type>>
volumeIntegrate(const tmp<DimensionedField<Type, fvMesh>>& tdf)
{
    tmp<Field<Type>> tdidf
    (
        tdf().mesh().V().primitiveField()*tdf().primitiveField()
    );
    tdf.clear();
    return tdidf;
}


template<class Type>
dimensioned<Type> domainIntegrate(const DimensionedField<Type, fvMesh>& df)
{
    return dimensioned<Type>
    (
        "domainIntegrate(" + df.name() + ')',
        dimensions::volume*df.dimensions(),
        gSum(fvi::volumeIntegrate(df))
    );
}


template<class Expression, class>
ElementType<Expression> domainIntegrate(const Expression& e)
{
    const fvMesh& mesh = expression::getFirst<expression::Mesh<fvMesh>>(e);

    return ElementType<Expression>
    (
        "domainIntegrate(" + expression::name(e) + ')',
        dimensions::volume*expression::access(e, dimensions::invalid),
        gSum
        (
            mesh.V().primitiveField()
           *expression::access
            (
                e,
                expression::Value(),
                expression::Base()
            )
        )
    );
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const tmp<DimensionedField<Type, fvMesh>>& tdf
)
{
    dimensioned<Type> integral = domainIntegrate(tdf());
    tdf.clear();
    return integral;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
