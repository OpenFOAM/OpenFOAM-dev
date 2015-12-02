/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds = field.dimensions()/dimTime*dimVolume;

    tmp<fvMatrix<Type> > tmtx(new fvMatrix<Type>(field, ds));
    fvMatrix<Type>& mtx = tmtx();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << fieldName << endl;
                }

                source.addSup(mtx, fieldI);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds
    (
        rho.dimensions()*field.dimensions()/dimTime*dimVolume
    );

    tmp<fvMatrix<Type> > tmtx(new fvMatrix<Type>(field, ds));
    fvMatrix<Type>& mtx = tmtx();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << fieldName << endl;
                }

                source.addSup(rho, mtx, fieldI);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(alpha, rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds
    (
        alpha.dimensions()*rho.dimensions()*field.dimensions()
       /dimTime*dimVolume
    );

    tmp<fvMatrix<Type> > tmtx(new fvMatrix<Type>(field, ds));
    fvMatrix<Type>& mtx = tmtx();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << fieldName << endl;
                }

                source.addSup(alpha, rho, mtx, fieldI);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const geometricOneField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
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
        dimensionedScalar("one", dimless, 1.0)
    );

    return this->operator()(alpha, one, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::fv::optionList::operator()
(
    const geometricOneField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(rho, field, field.name());
}


template<class Type>
void Foam::fv::optionList::constrain(fvMatrix<Type>& eqn)
{
    checkApplied();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldI = source.applyToField(eqn.psi().name());

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying constraint " << source.name()
                        << " to field " << eqn.psi().name() << endl;
                }

                source.constrain(eqn, fieldI);
            }
        }
    }
}


template<class Type>
void Foam::fv::optionList::correct
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    const word& fieldName = field.name();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Correcting source " << source.name()
                        << " for field " << fieldName << endl;
                }

                source.correct(field);
            }
        }
    }
}


// ************************************************************************* //
