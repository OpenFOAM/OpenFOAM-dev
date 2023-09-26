/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "fvModels.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class ... AlphaRhoFieldTypes>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::sourceTerm
(
    const VolField<Type>& eqnField,
    const dimensionSet& ds,
    const AlphaRhoFieldTypes& ... alphaRhoFields
) const
{
    checkApplied();

    tmp<fvMatrix<Type>> tmtx
    (
        new fvMatrix<Type>
        (
            eqnField,
            fvModel::sourceDims(ds, alphaRhoFields ...)
        )
    );
    fvMatrix<Type>& mtx = tmtx.ref();

    const PtrListDictionary<fvModel>& modelList(*this);

    const word& fieldName = fvModel::fieldName(alphaRhoFields ...);

    forAll(modelList, i)
    {
        const fvModel& model = modelList[i];

        if (model.addsSupToField(fieldName))
        {
            addSupFields_[i].insert(fieldName);

            if (debug)
            {
                Info<< "Applying model " << model.name() << " to field "
                    << fieldName << endl;
            }

            model.addSup(alphaRhoFields ..., mtx);
        }
    }

    return tmtx;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::sourceProxy
(
    const VolField<Type>& eqnField
) const
{
    return sourceTerm(eqnField, dimVolume/dimTime);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::source
(
    const VolField<Type>& field
) const
{
    return sourceTerm(field, dimVolume/dimTime, field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::sourceProxy
(
    const VolField<Type>& field,
    const VolField<Type>& eqnField
) const
{
    return sourceTerm(eqnField, dimVolume/dimTime, field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::source
(
    const volScalarField& rho,
    const VolField<Type>& field
) const
{
    return sourceTerm(field, dimVolume/dimTime, rho, field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::sourceProxy
(
    const volScalarField& rho,
    const VolField<Type>& field,
    const VolField<Type>& eqnField
) const
{
    return sourceTerm(eqnField, dimVolume/dimTime, rho, field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::source
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field
) const
{
    return sourceTerm(field, dimVolume/dimTime, alpha, rho, field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::sourceProxy
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    const VolField<Type>& eqnField
) const
{
    return sourceTerm(eqnField, dimVolume/dimTime, alpha, rho, field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::source
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const VolField<Type>& field
) const
{
    return this->source(field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::source
(
    const volScalarField& alpha,
    const geometricOneField& rho,
    const VolField<Type>& field
) const
{
    return this->source(alpha, field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::source
(
    const geometricOneField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field
) const
{
    return this->source(rho, field);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModels::d2dt2
(
    const VolField<Type>& field
) const
{
    return sourceTerm(field, dimVolume/sqr(dimTime), field);
}


// ************************************************************************* //
