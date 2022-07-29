/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "psiuMulticomponentThermo.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedUnburntEnthalpyFvPatchScalarField.H"
#include "gradientUnburntEnthalpyFvPatchScalarField.H"
#include "mixedUnburntEnthalpyFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiuMulticomponentThermo, 0);
    defineRunTimeSelectionTable(psiuMulticomponentThermo, fvMesh);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::psiuMulticomponentThermo::heuBoundaryTypes()
{
    const volScalarField::Boundary& tbf =
        this->Tu().boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedUnburntEnthalpyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientUnburntEnthalpyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedUnburntEnthalpyFvPatchScalarField::typeName;
        }
    }

    return hbt;
}

void Foam::psiuMulticomponentThermo::heuBoundaryCorrection(volScalarField& heu)
{
    volScalarField::Boundary& hbf = heu.boundaryFieldRef();

    forAll(hbf, patchi)
    {
        if
        (
            isA<gradientUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
        )
        {
            refCast<gradientUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
                .gradient() = hbf[patchi].fvPatchField::snGrad();
        }
        else if
        (
            isA<mixedUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
        )
        {
            refCast<mixedUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
                .refGrad() = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiuMulticomponentThermo::implementation::implementation
(
    const fvMesh& mesh,
    const word& phaseName
)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiuMulticomponentThermo>
Foam::psiuMulticomponentThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<psiuMulticomponentThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiuMulticomponentThermo::~psiuMulticomponentThermo()
{}


Foam::psiuMulticomponentThermo::implementation::~implementation()
{}


// ************************************************************************* //
