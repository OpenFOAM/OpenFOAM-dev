/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "CloudLagrangianFieldSource.H"
#include "LagrangianModel.H"
#include "cloud.H"
#include "CloudTypes.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
template<class Derived>
Foam::CloudLagrangianFieldSource<Type>::CloudLagrangianFieldSource
(
    const Derived& field
)
:
    field_(static_cast<const LagrangianFieldSource<Type>&>(field)),
    cloud_(field_.db().template lookupType<Foam::cloud>())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
template<class Cloud, class ... Clouds>
bool Foam::CloudLagrangianFieldSource<Type>::isCloud() const
{
    return CloudTypes<Cloud, Clouds ...>::isA(cloud_);
}


template<class Type>
template<class Cloud, class ... Clouds>
void Foam::CloudLagrangianFieldSource<Type>::assertCloud
(
    const LagrangianModel& model,
    const LagrangianSubMesh& subMesh
) const
{
    if (!isCloud<Cloud, Clouds ...>())
    {
        FatalErrorInFunction
            << "The '" << field_.type() << "' source of field '"
            << (field_.db().dbDir()/field_.internalField().name()).c_str()
            << "' for the '" << model.type() << "' Lagrangian model '"
            << model.name() << "' of cloud '" << cloud_.mesh().name()
            << "' requires a cloud of type "
            << CloudTypes<Cloud, Clouds ...>::typesString("or").c_str()
            << " (or a derivation thereof), rather than '" << cloud_.type()
            << "'" << exit(FatalError);
    }
}


template<class LagrangianFieldSource>
template<class Cloud>
const Cloud& Foam::CloudLagrangianFieldSource<LagrangianFieldSource>::cloud
(
    const LagrangianModel& model,
    const LagrangianSubMesh& subMesh
) const
{
    assertCloud<Cloud>(model, subMesh);
    return refCast<const Cloud>(cloud_);
}


// ************************************************************************* //
