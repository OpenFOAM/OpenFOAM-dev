/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
Foam::tmp<Foam::Field<Type>>
Foam::uniformInletOutletFvFieldSource<Type>::sourceValue
(
    const fvSource& source
) const
{
    const scalar t = this->db().time().value();
    const Type v = uniformInletValue_->value(t);
    return tmp<Field<Type>>(new Field<Type>(source.nCells(), v));
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::uniformInletOutletFvFieldSource<Type>::internalCoeff
(
    const fvSource& source
) const
{
    return
        neg0(source.source(this->internalField().name()))
       *scalarField(source.nCells(), scalar(1));
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
