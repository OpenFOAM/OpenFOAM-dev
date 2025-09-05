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

#include "LagrangianFieldSourceBase.H"
#include "LagrangianModel.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class OtherFieldSourceType>
const OtherFieldSourceType&
Foam::LagrangianFieldSourceBase::fieldSourceCast
(
    const LagrangianModel& model
) const
{
    if (!isA<OtherFieldSourceType>(*this))
    {
        FatalErrorInFunction
            << "The '" << type() << "' source of field '"
            << (db().dbDir()/internalIo_.name()).c_str()
            << "' for the '" << model.type() << "' Lagrangian model '"
            << model.name() << "' requires the corresponding source of field '"
            << (db().dbDir()/internalName()).c_str()
            << "' to be of type '" << OtherFieldSourceType::typeName
            << "' (or a derivation thereof), rather than '" << type()
            << "'" << exit(FatalError);
    }

    return refCast<const OtherFieldSourceType>(*this);
}


template<class OtherType>
const Foam::LagrangianFieldSource<OtherType>&
Foam::LagrangianFieldSourceBase::fieldSource
(
    const word& name,
    const LagrangianModel& model
) const
{
    const LagrangianDynamicField<OtherType>& lf =
        db().template lookupObject<LagrangianDynamicField<OtherType>>(name);

    return lf.sources()[model.name()];
}


template<class OtherType, class OtherFieldSourceType>
const OtherFieldSourceType& Foam::LagrangianFieldSourceBase::fieldSourceCast
(
    const word& name,
    const LagrangianModel& model
) const
{
    const LagrangianFieldSource<OtherType>& lfs =
        fieldSource<OtherType>(name, model);

    if (!isA<OtherFieldSourceType>(lfs))
    {
        FatalErrorInFunction
            << "The '" << type() << "' source of field '"
            << (db().dbDir()/internalIo_.name()).c_str()
            << "' for the '" << model.type() << "' Lagrangian model '"
            << model.name() << "' requires the corresponding source of field '"
            << (db().dbDir()/name).c_str()
            << "' to be of type '" << OtherFieldSourceType::typeName
            << "' (or a derivation thereof), rather than '" << lfs.type()
            << "'" << exit(FatalError);
    }

    return refCast<const OtherFieldSourceType>(lfs);
}


template<class OtherModelType>
const OtherModelType& Foam::LagrangianFieldSourceBase::modelCast
(
    const LagrangianModel& model
) const
{
    if (!isA<OtherModelType>(model))
    {
        FatalErrorInFunction
            << "The '" << type() << "' source of field '"
            << (db().dbDir()/internalIo_.name()).c_str()
            << "' for the Lagrangian model '" << model.name()
            << "' requires a model of type '" << OtherModelType::typeName
            << "' (or a derivation thereof), rather than '" << model.type()
            << "'" << exit(FatalError);
    }

    return refCast<const OtherModelType>(model);
}


template<class OtherType>
Foam::tmp<Foam::LagrangianSubField<OtherType>>
Foam::LagrangianFieldSourceBase::value
(
    const word& name,
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    return fieldSource<OtherType>(name, source).value(source, subMesh);
}


template<class OtherType>
Foam::tmp<Foam::LagrangianSubField<OtherType>>
Foam::LagrangianFieldSourceBase::value
(
    const word& name,
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return fieldSource<OtherType>(name, injection).value(injection, subMesh);
}


// ************************************************************************* //
