/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "UniformTable1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::UniformTable<Type>::UniformTable
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction1<Type, UniformTable<Type>>(name),
    dictName_(dict.name()),
    low_(dict.lookup<scalar>("low")),
    high_(dict.lookup<scalar>("high")),
    values_(dict.lookup("values"))
{
    if (values_.size() < 2)
    {
        FatalErrorInFunction
            << "Table " << nl
            << "    " << dictName_ << nl
            << "    has less than 2 entries."
            << exit(FatalError);
    }
    else
    {
        delta_ = (high_ - low_)/(values_.size() - 1);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1s::UniformTable<Type>::value(scalar x) const
{
    const scalar nd = (x - low_)/delta_;
    const label i = nd;

    if (nd < 0 || i > values_.size() - 2)
    {
        FatalErrorInFunction
            << x << " out of range "
            << low_ << " to " << high_ << nl
            << "    of table " << dictName_
            << exit(FatalError);
    }

    const scalar xi = low_ + i*delta_;
    const scalar lambda = (x - xi)/delta_;

    return values_[i] + lambda*(values_[i + 1] - values_[i]);
}


template<class Type>
Type Foam::Function1s::UniformTable<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return Zero;
}


template<class Type>
void Foam::Function1s::UniformTable<Type>::write(Ostream& os) const
{
    writeEntry(os, "low", low_);
    writeEntry(os, "high", high_);
    writeEntry(os, "values", values_);
}


// ************************************************************************* //
