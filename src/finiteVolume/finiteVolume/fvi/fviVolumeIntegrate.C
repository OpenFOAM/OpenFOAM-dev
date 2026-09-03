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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvi::volumeIntegrate
(
    const VolInternalField<Type>& vf
)
{
    return vf.mesh().V().primitiveField()*vf.primitiveField();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvi::volumeIntegrate
(
    const tmp<VolInternalField<Type>>& tvf
)
{
    tmp<Field<Type>> tvivf(volumeIntegrate(tvf()));
    tvf.clear();
    return tvivf;
}


template<class Type>
Foam::dimensioned<Type> Foam::fvi::domainIntegrate
(
    const VolInternalField<Type>& vf
)
{
    return dimensioned<Type>
    (
        "domainIntegrate(" + vf.name() + ')',
        dimensions::volume*vf.dimensions(),
        gSum(volumeIntegrate(vf))
    );
}


template<class Type>
Foam::dimensioned<Type> Foam::fvi::domainIntegrate
(
    const tmp<VolInternalField<Type>>& tvf
)
{
    dimensioned<Type> integral = domainIntegrate(tvf());
    tvf.clear();
    return integral;
}


template<class Expression, class>
Foam::ElementType<Expression> Foam::fvi::domainIntegrate(const Expression& e)
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


// ************************************************************************* //
