/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "zeroDimensionalMassSourceBase.H"
#include "fvCellSet.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(zeroDimensionalMassSourceBase, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::zeroDimensionalMassSourceBase::calcM0D() const
{
    tmp<volScalarField> tm =
        volScalarField::New
        (
            typedName("m0D"),
            mesh(),
            dimensionedScalar(dimMass, 0)
        );

    HashTable<const basicThermo*> thermos(mesh().lookupClass<basicThermo>());

    forAllConstIter(HashTable<const basicThermo*>, thermos, thermoIter)
    {
        const basicThermo& thermo = *thermoIter();

        tmp<volScalarField> tRho = thermo.rho();
        const volScalarField& rho = tRho();

        const word phaseName = thermo.phaseName();

        if (thermo.phaseName() != word::null)
        {
            const volScalarField& alpha =
                mesh().lookupObject<volScalarField>
                (
                    IOobject::groupName("alpha", phaseName)
                );

            tm.ref().ref() += alpha()*rho()*mesh().V();
        }
        else
        {
            tm.ref().ref() += rho()*mesh().V();
        }
    }

    return tm;
}


Foam::volScalarField& Foam::fv::zeroDimensionalMassSourceBase::initM0D() const
{
    if (!mesh().foundObject<volScalarField>(typedName("m0D")))
    {
        volScalarField* mPtr =
            new volScalarField
            (
                calcM0D()
            );

        mPtr->store();
    }

    return mesh().lookupObjectRef<volScalarField>(typedName("m0D"));
}


const Foam::volScalarField& Foam::fv::zeroDimensionalMassSourceBase::m() const
{
    // If not registered, then read or create the mass field
    if (!mesh().foundObject<volScalarField>(typedName("m")))
    {
        typeIOobject<volScalarField> mIo
        (
            typedName("m"),
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        volScalarField* mPtr =
            new volScalarField
            (
                mIo,
                mesh(),
                dimensionedScalar(dimMass, 0)
            );

        mPtr->store();

        if (!mIo.headerOk())
        {
            *mPtr = m0D_;
        }

        volScalarField* factorPtr =
            new volScalarField
            (
                IOobject
                (
                    typedName("factor"),
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                *mPtr/m0D_
            );

        factorPtr->store();
    }

    volScalarField& m =
        mesh().lookupObjectRef<volScalarField>(typedName("m"));

    volScalarField& factor =
        mesh().lookupObjectRef<volScalarField>(typedName("factor"));

    // Update the mass if changes are available
    if (mesh().foundObject<volScalarField>(typedName("deltaM")))
    {
        volScalarField& deltaM =
            mesh().lookupObjectRef<volScalarField>(typedName("deltaM"));

        m = m.oldTime() + deltaM;

        factor = m/m0D_;

        deltaM.checkOut();
    }

    return m;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::zeroDimensionalMassSourceBase::zeroDimensionalMassSourceBase
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    massSourceBase(name, modelType, mesh, dict),
    m0D_(initM0D())
{
    if (mesh.nGeometricD() != 0)
    {
        FatalIOErrorInFunction(dict)
            << "Zero-dimensional fvModel applied to a "
            << mesh.nGeometricD() << "-dimensional mesh"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelUList Foam::fv::zeroDimensionalMassSourceBase::cells() const
{
    static labelList zero(1, 0);
    return labelUList(zero);
}


Foam::label Foam::fv::zeroDimensionalMassSourceBase::nCells() const
{
    return 1;
}


Foam::scalar Foam::fv::zeroDimensionalMassSourceBase::V() const
{
    return mesh().V()[0];
}


Foam::dimensionedScalar Foam::fv::zeroDimensionalMassSourceBase::S() const
{
    return
        dimensionedScalar
        (
            dimMass/dimTime,
            massFlowRate()*m0D_[0]/m()[0]
        );
}


bool Foam::fv::zeroDimensionalMassSourceBase::movePoints()
{
    return true;
}


void Foam::fv::zeroDimensionalMassSourceBase::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fv::zeroDimensionalMassSourceBase::mapMesh
(
    const polyMeshMap& map
)
{}


void Foam::fv::zeroDimensionalMassSourceBase::distribute
(
    const polyDistributionMap& map
)
{}


void Foam::fv::zeroDimensionalMassSourceBase::correct()
{
    // Correct the zero-dimensional mass
    m0D_ = calcM0D();

    // Create the mass change
    if (!mesh().foundObject<volScalarField>(typedName("deltaM")))
    {
        volScalarField* dMPtr =
            new volScalarField
            (
                IOobject
                (
                    typedName("deltaM"),
                    mesh().time().timeName(),
                    mesh()
                ),
                mesh(),
                dimensionedScalar(dimMass, 0)
            );

        dMPtr->store();
    }

    volScalarField& deltaM =
        mesh().lookupObjectRef<volScalarField>(typedName("deltaM"));

    deltaM +=
        mesh().time().deltaT()
       *dimensionedScalar(dimMass/dimTime, massFlowRate());
}


// ************************************************************************* //
