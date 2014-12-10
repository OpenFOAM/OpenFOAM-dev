/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "nonEquilibrium.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TimeScaleModels
{
    defineTypeNameAndDebug(nonEquilibrium, 0);

    addToRunTimeSelectionTable
    (
        TimeScaleModel,
        nonEquilibrium,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TimeScaleModels::nonEquilibrium::nonEquilibrium
(
    const dictionary& dict
)
:
    TimeScaleModel(dict)
{}


Foam::TimeScaleModels::nonEquilibrium::nonEquilibrium
(
    const nonEquilibrium& hc
)
:
    TimeScaleModel(hc)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TimeScaleModels::nonEquilibrium::~nonEquilibrium()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar> >
Foam::TimeScaleModels::nonEquilibrium::oneByTau
(
    const FieldField<Field, scalar>& alpha,
    const FieldField<Field, scalar>& r32,
    const FieldField<Field, scalar>& uSqr,
    const FieldField<Field, scalar>& f
) const
{
    static const scalar a =
        8.0*sqrt(2.0)/(3.0*constant::mathematical::pi)
       *0.25*(1.0 - e_*e_);

    return a*f*alphaPacked_/max(alphaPacked_ - alpha, SMALL);
}


// ************************************************************************* //
