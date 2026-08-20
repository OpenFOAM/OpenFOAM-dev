/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
dimensioned<Type> average
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df
)
{
    dimensioned<Type> Average
    (
        "average(" + df.name() + ')',
        df.dimensions(),
        gAverage(df.primitiveField())
    );

    return Average;
}


template
<
    class Type,
    class GeoMesh,
    template<class> class PrimitiveField,
    template<class> class PrimitiveField2
>
dimensioned<Type> weightedAverage
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& df,
    const DimensionedField<scalar, GeoMesh, PrimitiveField2>& weightField
)
{
    return
    (
        dimensioned<Type>
        (
            "weightedAverage(" + df.name() + ',' + weightField.name() + ')',
            df.dimensions(),
            gSum(weightField.primitiveField()*df.primitiveField())
           /gSum(weightField.primitiveField())
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
