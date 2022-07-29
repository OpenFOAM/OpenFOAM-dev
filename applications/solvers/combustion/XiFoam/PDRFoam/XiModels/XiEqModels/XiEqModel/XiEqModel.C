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

#include "XiEqModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XiEqModel, 0);
    defineRunTimeSelectionTable(XiEqModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModel::XiEqModel
(
    const dictionary& XiEqProperties,
    const psiuMulticomponentThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModelCoeffs_
    (
        XiEqProperties.subDict
        (
            word(XiEqProperties.lookup("XiEqModel")) + "Coeffs"
        )
    ),
    thermo_(thermo),
    turbulence_(turbulence),
    Su_(Su)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiEqModel::~XiEqModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::XiEqModel::read(const dictionary& XiEqProperties)
{
    XiEqModelCoeffs_ = XiEqProperties.optionalSubDict(type() + "Coeffs");

    return true;
}


void Foam::XiEqModel::writeFields() const
{
    //***HGW It is not clear why B is written here
    if (Su_.mesh().foundObject<volSymmTensorField>("B"))
    {
        const volSymmTensorField& B =
            Su_.mesh().lookupObject<volSymmTensorField>("B");
        B.write();
    }
}


Foam::tmp<Foam::volScalarField>
Foam::XiEqModel::calculateSchelkinEffect(const scalar uPrimeCoef) const
{
    const fvMesh& mesh = Su_.mesh();

    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const volSymmTensorField& CT = mesh.lookupObject<volSymmTensorField>("CT");
    const volScalarField& Nv = mesh.lookupObject<volScalarField>("Nv");
    const volSymmTensorField& nsv =
        mesh.lookupObject<volSymmTensorField>("nsv");

    tmp<volScalarField> tN
    (
        volScalarField::New
        (
            "tN",
            mesh,
            dimensionedScalar(Nv.dimensions(), 0)
        )
    );
    volScalarField& N = tN.ref();
    N.primitiveFieldRef() = Nv.primitiveField()*pow(mesh.V(), 2.0/3.0);

    volSymmTensorField ns
    (
        IOobject
        (
            "tns",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor
        (
            "zero",
            nsv.dimensions(),
            Zero
        )
    );
    ns.primitiveFieldRef() = nsv.primitiveField()*pow(mesh.V(), 2.0/3.0);

    const volVectorField Uhat
    (
        U/(mag(U) + dimensionedScalar(U.dimensions(), 1e-4))
    );

    const volScalarField nr(sqrt(max(N - (Uhat & ns & Uhat), scalar(1e-4))));

    const scalarField cellWidth(pow(mesh.V(), 1.0/3.0));

    const scalarField upLocal(uPrimeCoef*sqrt((U & CT & U)*cellWidth));

    const scalarField deltaUp(upLocal*(max(scalar(1), pow(nr, 0.5)) - 1.0));

    // Re use tN
    N.primitiveFieldRef() = upLocal*(max(scalar(1), pow(nr, 0.5)) - 1.0);

    return tN;
}


// ************************************************************************* //
