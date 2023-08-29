/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

#include "phaseChange.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace driftModels
{
    defineTypeNameAndDebug(phaseChange, 0);
    addToRunTimeSelectionTable(driftModel, phaseChange, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::driftModels::phaseChange::phaseChange
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    driftModel(popBal, dict),
    interfaces_
    (
        dict.lookup("interfaces"),
        phaseInterface::iNew(popBal_.fluid())
    ),
    numberWeighted_(dict.lookupOrDefault<Switch>("numberWeighted", false)),
    W_(interfaces_.size()),
    dmdtfName_(dict.lookup("dmdtf")),
    specieName_(dict.lookupOrDefault("specie", word()))
{
    forAll(interfaces_, i)
    {
        W_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(typedName("W"), interfaces_[i].name()),
                    popBal_.mesh().time().name(),
                    popBal_.mesh()
                ),
                popBal_.mesh(),
                dimensionedScalar
                (
                    inv(numberWeighted_ ? dimVolume : dimLength),
                    Zero
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::driftModels::phaseChange::precompute()
{
    forAll(interfaces_, i)
    {
        W_[i] = Zero;
    }

    forAll(interfaces_, k)
    {
        forAllConstIter(phaseInterface, interfaces_[k], iter)
        {
            const phaseModel& phase = iter();

            if (!isA<velocityGroup>(phase.diameter())) continue;

            const velocityGroup& velGroup =
                refCast<const velocityGroup>(phase.diameter());

            forAll(velGroup.sizeGroups(), i)
            {
                const sizeGroup& fi = velGroup.sizeGroups()[i];

                if (numberWeighted_)
                {
                    W_[k] += fi*max(fi.phase(), small)/fi.x();
                }
                else
                {
                    W_[k] += fi*max(fi.phase(), small)/fi.x()*fi.a();
                }
            }
        }
    }
}


void Foam::diameterModels::driftModels::phaseChange::addToDriftRate
(
    volScalarField& driftRate,
    const label i
)
{
    const velocityGroup& velGrp = popBal_.sizeGroups()[i].group();

    forAll(interfaces_, k)
    {
        if (interfaces_[k].contains(velGrp.phase()))
        {
            const volScalarField& dmidtf =
                popBal_.mesh().lookupObject<volScalarField>
                (
                    IOobject::groupName
                    (
                        IOobject::groupName
                        (
                            dmdtfName_,
                            specieName_
                        ),
                        interfaces_[k].name()
                    )
                );

            const scalar dmidtfSign =
                interfaces_[k].index(velGrp.phase()) == 0 ? +1 : -1;

            const sizeGroup& fi = popBal_.sizeGroups()[i];

            tmp<volScalarField> dDriftRate
            (
                dmidtfSign*dmidtf/(fi.phase().rho()*W_[k])
            );

            if (!numberWeighted_)
            {
                dDriftRate.ref() *= fi.a();
            }

            driftRate += velGrp.phase()*dDriftRate;
        }
    }
}


// ************************************************************************* //
