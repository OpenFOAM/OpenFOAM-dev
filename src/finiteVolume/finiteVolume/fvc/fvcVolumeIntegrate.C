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

#include "fvcVolumeIntegrate.H"
#include "fviVolumeIntegrate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvc::volumeIntegrate
(
    const VolField<Type>& vf
)
{
    return vf.mesh().V().primitiveField()*vf.primitiveField();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvc::volumeIntegrate
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<Field<Type>> tvivf(volumeIntegrate(tvf()));
    tvf.clear();
    return tvivf;
}


template<class Type>
Foam::dimensioned<Type> Foam::fvc::domainIntegrate(const VolField<Type>& vf)
{
    return fvi::domainIntegrate(vf.internalField());
}


template<class Type>
Foam::dimensioned<Type> Foam::fvc::domainIntegrate
(
    const tmp<VolField<Type>>& tvf
)
{
    dimensioned<Type> integral(domainIntegrate(tvf()));
    tvf.clear();
    return integral;
}


template<class Expression, class>
Foam::ElementType<Expression> Foam::fvc::domainIntegrate(const Expression& e)
{
    return fvi::domainIntegrate
    (
        expression::access(e, GeometricField_InternalField())
    );
}


// ************************************************************************* //
