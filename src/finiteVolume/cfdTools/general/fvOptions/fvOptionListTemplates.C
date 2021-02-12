/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class ... AlphaRhoFieldTypes>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::source
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName,
    const dimensionSet& ds,
    const AlphaRhoFieldTypes& ... alphaRhos
) const
{
    checkApplied();

    tmp<fvMatrix<Type>> tmtx
    (
        new fvMatrix<Type>
        (
            field,
            option::sourceDims(field, ds, alphaRhos ...)
        )
    );
    fvMatrix<Type>& mtx = tmtx.ref();

    forAll(*this, i)
    {
        const option& source = this->operator[](i);

        if (source.addsSupToField(fieldName))
        {
            addSupFields_[i].insert(fieldName);

            if (debug)
            {
                Info<< "Applying source " << source.name() << " to field "
                    << fieldName << endl;
            }

            source.addSup(alphaRhos ..., mtx, fieldName);
        }
    }

    return tmtx;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->operator()(rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime, rho);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->operator()(alpha, rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime, alpha, rho);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const geometricOneField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    volScalarField one
    (
        IOobject
        (
            "one",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    );

    return this->operator()(alpha, one, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const geometricOneField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->operator()(rho, field, field.name());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::d2dt2
(
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->d2dt2(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::d2dt2
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/sqr(dimTime));
}


template<class Type>
void Foam::fv::optionList::constrain(fvMatrix<Type>& eqn) const
{
    checkApplied();

    forAll(*this, i)
    {
        const option& constraint = this->operator[](i);

        if (constraint.constrainsField(eqn.psi().name()))
        {
            constrainedFields_[i].insert(eqn.psi().name());

            if (debug)
            {
                Info<< "Applying constraint " << constraint.name()
                    << " to field " << eqn.psi().name() << endl;
            }

            constraint.constrain(eqn, eqn.psi().name());
        }
    }
}


template<class Type>
void Foam::fv::optionList::correct
(
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    const word& fieldName = field.name();

    forAll(*this, i)
    {
        const option& correction = this->operator[](i);

        if (correction.correctsField(fieldName))
        {
            correctedFields_[i].insert(fieldName);

            if (debug)
            {
                Info<< "Applying correction " << correction.name()
                    << " for field " << fieldName << endl;
            }

            correction.correct(field);
        }
    }
}


// ************************************************************************* //
