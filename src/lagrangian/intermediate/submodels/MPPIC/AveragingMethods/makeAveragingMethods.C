/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "Field.H"
#include "fvcGrad.H"
#include "polyMeshTetDecomposition.H"

#include "Basic.H"
#include "Dual.H"
#include "Moment.H"

// Scalar interpolation
defineNamedTemplateTypeNameAndDebug(Foam::AveragingMethod<Foam::scalar>, 0);
namespace Foam
{
    defineTemplateRunTimeSelectionTable
    (
        AveragingMethod<Foam::scalar>,
        dictionary
    );
}

// Vector interpolation
defineNamedTemplateTypeNameAndDebug(Foam::AveragingMethod<Foam::vector>, 0);
namespace Foam
{
    defineTemplateRunTimeSelectionTable
    (
        Foam::AveragingMethod<Foam::vector>,
        dictionary
    );
}


// Basic interpolation
defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Basic<Foam::scalar>,
    0
);
Foam::AveragingMethod<Foam::scalar>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Basic<Foam::scalar>>
    addBasicscalarConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Basic<Foam::vector>,
    0
);
Foam::AveragingMethod<Foam::vector>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Basic<Foam::vector>>
    addBasicvectorConstructorToTable_;


// Dual interpolation
defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Dual<Foam::scalar>,
    0
);
Foam::AveragingMethod<Foam::scalar>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Dual<Foam::scalar>>
    addDualscalarConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Dual<Foam::vector>,
    0
);
Foam::AveragingMethod<Foam::vector>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Dual<Foam::vector>>
    addDualvectorConstructorToTable_;


// Moment interpolation
defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Moment<Foam::scalar>,
    0
);
Foam::AveragingMethod<Foam::scalar>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Moment<Foam::scalar>>
    addMomentscalarConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    Foam::AveragingMethods::Moment<Foam::vector>,
    0
);
Foam::AveragingMethod<Foam::vector>::
adddictionaryConstructorToTable<Foam::AveragingMethods::Moment<Foam::vector>>
    addMomentvectorConstructorToTable_;


// ************************************************************************* //
