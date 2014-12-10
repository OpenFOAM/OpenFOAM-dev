/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "basic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDRDragModels
{
    defineTypeNameAndDebug(basic, 0);
    addToRunTimeSelectionTable(PDRDragModel, basic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRDragModels::basic::basic
(
    const dictionary& PDRProperties,
    const compressible::RASModel& turbulence,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    PDRDragModel(PDRProperties, turbulence, rho, U, phi),
    Csu("Csu", dimless, PDRDragModelCoeffs_.lookup("Csu")),
    Csk("Csk", dimless, PDRDragModelCoeffs_.lookup("Csk")),

    Aw_
    (
        IOobject
        (
            "Aw",
            U_.mesh().facesInstance(),
            U_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh()
    ),

    CR_
    (
        IOobject
        (
            "CR",
            U_.mesh().facesInstance(),
            U_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh()
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::PDRDragModels::basic::~basic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField> Foam::PDRDragModels::basic::Dcu() const
{
    tmp<volSymmTensorField> tDragDcu
    (
        new volSymmTensorField
        (
            IOobject
            (
                "tDragDcu",
                U_.mesh().time().constant(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedSymmTensor
            (
                "zero",
                dimMass/dimTime/pow(dimLength, 3),
                pTraits<symmTensor>::zero
            ),
            zeroGradientFvPatchSymmTensorField::typeName
        )
    );

    volSymmTensorField& DragDcu = tDragDcu();

    if (on_)
    {
        const volScalarField& betav =
            U_.db().lookupObject<volScalarField>("betav");

        DragDcu =
            (0.5*rho_)*CR_*mag(U_) + (Csu*I)*betav*turbulence_.muEff()*sqr(Aw_);
    }

    return tDragDcu;
}


Foam::tmp<Foam::volScalarField> Foam::PDRDragModels::basic::Gk() const
{
    tmp<volScalarField> tGk
    (
        new volScalarField
        (
            IOobject
            (
                "tGk",
                U_.mesh().time().constant(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimMass/dimLength/pow(dimTime, 3), 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    volScalarField& Gk = tGk();

    if (on_)
    {
        const volScalarField& betav =
            U_.db().lookupObject<volScalarField>("betav");

        const volSymmTensorField& CT =
            U_.db().lookupObject<volSymmTensorField>("CT");

        Gk =
            (0.5*rho_)*mag(U_)*(U_ & CT & U_)
          + Csk*betav*turbulence_.muEff()*sqr(Aw_)*magSqr(U_);
    }

    return tGk;
}


bool Foam::PDRDragModels::basic::read(const dictionary& PDRProperties)
{
    PDRDragModel::read(PDRProperties);

    PDRDragModelCoeffs_.lookup("Csu") >> Csu.value();
    PDRDragModelCoeffs_.lookup("Csk") >> Csk.value();

    return true;
}


void Foam::PDRDragModels::basic::writeFields() const
{
    Aw_.write();
    CR_.write();
}

// ************************************************************************* //
