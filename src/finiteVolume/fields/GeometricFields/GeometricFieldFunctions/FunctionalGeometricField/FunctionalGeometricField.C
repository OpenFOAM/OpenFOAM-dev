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

#include "FunctionalGeometricField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class GeoMesh>
void Foam::FunctionalGeometricField<Type, GeoMesh>::setPatchFields()
{
    forAll(this->mesh().boundary(), patchi)
    {
        if (!this->mesh().boundary()[patchi].constraint())
        {
            patchFieldPtrs_.set
            (
                patchi,
                new SlicedDimensionedField<Type, GeoPatch>
                (
                    IOobject
                    (
                        this->name(),
                        this->mesh().time().name(),
                        this->mesh().db(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    this->mesh().boundary()[patchi],
                    this->dimensions(),
                    this->boundaryFieldRef()[patchi]
                )
            );
        }
    }
}


template<class Type, class GeoMesh>
bool Foam::FunctionalGeometricField<Type, GeoMesh>::readFuncs
(
    const dictionary& dict
)
{
    if (!dict.isDict(funcName_)) return false;

    internalFuncPtr_.set
    (
        DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>::New
        (
            dict.subDict(funcName_),
            this->internalFieldRef()
        ).ptr()
    );

    setPatchFields();

    forAll(this->mesh().boundary(), patchi)
    {
        if (!this->mesh().boundary()[patchi].constraint())
        {
            patchFuncPtrs_.set
            (
                patchi,
                DimensionedFieldFunction<DimensionedField<Type, GeoPatch>>::New
                (
                    dict.subDict(funcName_),
                    patchFieldPtrs_[patchi]
                ).ptr()
            );
        }
    }

    return true;
}


template<class Type, class GeoMesh>
void Foam::FunctionalGeometricField<Type, GeoMesh>::cloneFuncs
(
    const FunctionalGeometricField<Type, GeoMesh>& udff
)
{
    if (!udff.internalFuncPtr_.valid()) return;

    internalFuncPtr_.set
    (
        udff.internalFuncPtr_.clone
        (
            this->internalFieldRef()
        ).ptr()
    );

    setPatchFields();

    forAll(this->mesh().boundary(), patchi)
    {
        if (!this->mesh().boundary()[patchi].constraint())
        {
            patchFuncPtrs_.set
            (
                patchi,
                udff.boundaryFuncPtrs_[patchi].clone
                (
                    patchFieldPtrs_[patchi]
                ).ptr()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::FunctionalGeometricField<Type, GeoMesh>::FunctionalGeometricField
(
    const word& name,
    const word& funcName,
    const GeoMesh& mesh,
    const dimensionSet& dimensions,
    const dictionary& dict
)
:
    GeometricField<Type, GeoMesh>
    (
        IOobject
        (
            name,
            mesh.time().name(),
            mesh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensions
    ),
    funcName_(funcName),
    internalFuncPtr_(),
    patchFieldPtrs_(mesh.boundary().size()),
    patchFuncPtrs_(mesh.boundary().size()),
    defaultValue_
    (
        readFuncs(dict)
      ? pTraits<Type>::nan
      : dimensioned<Type>(funcName, dimensions, dict).value()
    )
{
    if (internalFuncPtr_.valid())
    {
        if (mesh.time().completeCase())
        {
            internalFuncPtr_->evaluate();

            forAll(this->mesh().boundary(), patchi)
            {
                if (!this->mesh().boundary()[patchi].constraint())
                {
                    patchFuncPtrs_[patchi].evaluate();
                }
            }

            this->correctBoundaryConditions();
        }
    }
    else
    {
        *this == dimensioned<Type>(dimensions, defaultValue_);
    }
}


template<class Type, class GeoMesh>
Foam::FunctionalGeometricField<Type, GeoMesh>::FunctionalGeometricField
(
    const word& name,
    const word& funcName,
    const GeoMesh& mesh,
    const dimensionSet& dimensions,
    const dictionary& dict,
    const Type& defaultValue
)
:
    GeometricField<Type, GeoMesh>
    (
        IOobject
        (
            name,
            mesh.time().name(),
            mesh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensions,
        defaultValue
    ),
    funcName_(funcName),
    internalFuncPtr_(),
    patchFieldPtrs_(mesh.boundary().size()),
    patchFuncPtrs_(mesh.boundary().size()),
    defaultValue_
    (
        readFuncs(dict)
      ? pTraits<Type>::nan
      : dict.found(funcName)
      ? dimensioned<Type>(funcName, dimensions, dict).value()
      : defaultValue_
    )
{
    if (internalFuncPtr_.valid())
    {
        if (mesh.time().completeCase())
        {
            internalFuncPtr_->evaluate();

            forAll(this->mesh().boundary(), patchi)
            {
                if (!this->mesh().boundary()[patchi].constraint())
                {
                    patchFuncPtrs_[patchi].evaluate();
                }
            }

            this->correctBoundaryConditions();
        }
    }
    else
    {
        *this == dimensioned<Type>(dimensions, defaultValue_);
    }
}


template<class Type, class GeoMesh>
Foam::FunctionalGeometricField<Type, GeoMesh>::FunctionalGeometricField
(
    const FunctionalGeometricField<Type, GeoMesh>& udff,
    const GeoMesh& mesh
)
:
    GeometricField<Type, GeoMesh>(udff, mesh, udff.dimensions()),
    funcName_(udff.funcName_),
    internalFuncPtr_(),
    patchFieldPtrs_(mesh.boundary().size()),
    patchFuncPtrs_(mesh.boundary().size()),
    defaultValue_(udff.defaultValue_)
{
    cloneFuncs(udff);
}


template<class Type, class GeoMesh>
Foam::FunctionalGeometricField<Type, GeoMesh>::FunctionalGeometricField
(
    const FunctionalGeometricField<Type, GeoMesh>& udff
)
:
    GeometricField<Type, GeoMesh>(udff),
    funcName_(udff.funcName_),
    internalFuncPtr_(),
    patchFieldPtrs_(this->mesh().boundary().size()),
    patchFuncPtrs_(this->mesh().boundary().size()),
    defaultValue_(udff.defaultValue_)
{
    cloneFuncs(udff);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
void Foam::FunctionalGeometricField<Type, GeoMesh>::reset()
{
    GeometricField<Type, GeoMesh>::reset
    (
        GeometricField<Type, GeoMesh>
        (
            static_cast<IOobject>(*this),
            this->mesh(),
            dimensioned<Type>("NaN", this->dimensions(), pTraits<Type>::nan)
        )
    );

    setPatchFields();

    if (internalFuncPtr_.valid())
    {
        internalFuncPtr_->reset();

        forAll(this->mesh().boundary(), patchi)
        {
            if (!this->mesh().boundary()[patchi].constraint())
            {
                patchFuncPtrs_[patchi].reset();
            }
        }

        this->correctBoundaryConditions();
    }
    else
    {
        *this == dimensioned<Type>(this->dimensions(), defaultValue_);
    }
}


template<class Type, class GeoMesh>
bool Foam::FunctionalGeometricField<Type, GeoMesh>::update()
{
    if (internalFuncPtr_.valid())
    {
        bool result = false;

        result = internalFuncPtr_->update() || result;

        forAll(this->mesh().boundary(), patchi)
        {
            if (!this->mesh().boundary()[patchi].constraint())
            {
                result = patchFuncPtrs_[patchi].update() || result;
            }
        }

        if (result) this->correctBoundaryConditions();

        return result;
    }
    else
    {
        return false;
    }
}


template<class Type, class GeoMesh>
void Foam::FunctionalGeometricField<Type, GeoMesh>::write(Ostream& os) const
{
    if (internalFuncPtr_.valid())
    {
        internalFuncPtr_->write(os);
    }
}


template<class Type, class GeoMesh>
void Foam::writeEntry
(
    Ostream& os,
    const FunctionalGeometricField<Type, GeoMesh>& udff
)
{
    if (udff.internalFuncPtr_.valid())
    {
        writeEntry(os, udff.funcName_, udff.internalFuncPtr_());
    }
    else
    {
        writeEntry(os, udff.funcName_, udff.defaultValue_);
    }
}


// ************************************************************************* //
