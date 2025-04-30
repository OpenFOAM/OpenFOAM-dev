/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "DimensionedField.H"
#include "uniformInletOutletFvFieldSource.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformInletOutletFvFieldSource<Type>::uniformInletOutletFvFieldSource
(
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvFieldSource<Type>(iF, dict),
    uniformInletValue_
    (
        Function1<Type>::New
        (
            "uniformInletValue",
            this->db().time().userUnits(),
            iF.dimensions(),
            dict
        )
    )
{}


template<class Type>
Foam::uniformInletOutletFvFieldSource<Type>::uniformInletOutletFvFieldSource
(
    const uniformInletOutletFvFieldSource<Type>& field,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvFieldSource<Type>(field, iF),
    uniformInletValue_(field.uniformInletValue_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformInletOutletFvFieldSource<Type>::~uniformInletOutletFvFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::uniformInletOutletFvFieldSource<Type>::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return
        DimensionedField<Type, Foam::volMesh>::New
        (
            model.name() + ":" + this->internalField().name() + "SourceValue",
            this->internalField().mesh(),
            dimensioned<Type>
            (
                this->internalField().dimensions(),
                uniformInletValue_->value(this->db().time().value())
            )
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::uniformInletOutletFvFieldSource<Type>::sourceValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return
        tmp<Field<Type>>
        (
            new Field<Type>
            (
                source.size(),
                uniformInletValue_->value(this->db().time().value())
            )
        );
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::uniformInletOutletFvFieldSource<Type>::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return neg0(source);
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::uniformInletOutletFvFieldSource<Type>::internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return neg0(source);
}


template<class Type>
void Foam::uniformInletOutletFvFieldSource<Type>::write(Ostream& os) const
{
    fvFieldSource<Type>::write(os);
    writeEntry
    (
        os,
        this->db().time().userUnits(),
        this->internalField().dimensions(),
        uniformInletValue_()
    );
}


// ************************************************************************* //
