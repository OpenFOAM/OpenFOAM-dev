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

#include "basicXiSubXiEq.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(basicSubGrid, 0);
    addToRunTimeSelectionTable(XiEqModel, basicSubGrid, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::basicSubGrid::basicSubGrid
(
    const dictionary& XiEqProperties,
    const psiuMulticomponentThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, thermo, turbulence, Su),

    B_
    (
        IOobject
        (
            "B",
            Su.mesh().facesInstance(),
            Su.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        Su.mesh()
    ),

    XiEqModel_(XiEqModel::New(XiEqModelCoeffs_, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiEqModels::basicSubGrid::~basicSubGrid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::basicSubGrid::XiEq() const
{
    const fvMesh& mesh = Su_.mesh();
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    const volScalarField& Nv = mesh.lookupObject<volScalarField>("Nv");
    const volSymmTensorField& nsv =
        mesh.lookupObject<volSymmTensorField>("nsv");

    volScalarField magU(mag(U));
    volVectorField Uhat
    (
        U/(mag(U) + dimensionedScalar(U.dimensions(), 1e-4))
    );

    const scalarField Cw = pow(mesh.V(), 2.0/3.0);

    volScalarField N
    (
        IOobject
        (
            "N",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(Nv.dimensions(), 0)
    );
    N.primitiveFieldRef() = Nv.primitiveField()*Cw;

    volSymmTensorField ns
    (
        IOobject
        (
            "ns",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
            "zero",
            nsv.dimensions(),
            Zero
        )
    );
    ns.primitiveFieldRef() = nsv.primitiveField()*Cw;

    volScalarField n(max(N - (Uhat & ns & Uhat), scalar(1e-4)));
    volScalarField b((Uhat & B_ & Uhat)/sqrt(n));
    volScalarField up(sqrt((2.0/3.0)*turbulence_.k()));

    volScalarField XiSubEq
    (
        scalar(1)
      + max(2.2*sqrt(b), min(0.34*magU/up*sqrt(b), scalar(1.6)))
      * min(n, scalar(1))
    );

    return (XiSubEq*XiEqModel_->XiEq());
}


bool Foam::XiEqModels::basicSubGrid::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    return XiEqModel_->read(XiEqModelCoeffs_);
}


// ************************************************************************* //
