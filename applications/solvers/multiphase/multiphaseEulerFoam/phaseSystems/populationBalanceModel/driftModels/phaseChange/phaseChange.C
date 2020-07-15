/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
#include "phasePairKey.H"

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
    pairKeys_(dict.lookup("pairs")),
    numberWeighted_(dict.lookupOrDefault<Switch>("numberWeighted", false)),
    W_(pairKeys_.size()),
    dmdtfName_(dict.lookup("dmdtf")),
    specieName_(dict.lookupOrDefault("specie", word()))
{
    const phaseSystem& fluid = popBal_.fluid();

    forAll(pairKeys_, i)
    {
        const phasePair& pair = fluid.phasePairs()[pairKeys_[i]];

        W_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(type() + ":W", pair.name()),
                    popBal_.mesh().time().timeName(),
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

void Foam::diameterModels::driftModels::phaseChange::correct()
{
    const phaseSystem& fluid = popBal_.fluid();

    forAll(pairKeys_, i)
    {
        W_[i] = Zero;
    }

    forAll(pairKeys_, k)
    {
        if (fluid.phasePairs().found(pairKeys_[k]))
        {
            const phasePair& pair = fluid.phasePairs()[pairKeys_[k]];

            forAll(popBal_.velocityGroups(), j)
            {
                const velocityGroup& vgj = popBal_.velocityGroups()[j];
                if (pair.contains(vgj.phase()))
                {
                    forAll(vgj.sizeGroups(), i)
                    {
                        const sizeGroup& fi = vgj.sizeGroups()[i];

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
    }
}


void Foam::diameterModels::driftModels::phaseChange::addToDriftRate
(
    volScalarField& driftRate,
    const label i
)
{
    const velocityGroup& vg = popBal_.sizeGroups()[i].VelocityGroup();

    forAll(pairKeys_, k)
    {
        const phasePair& pair =
                popBal_.fluid().phasePairs()[pairKeys_[k]];

        if (pair.contains(vg.phase()))
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
                        pair.name()
                    )
                );

            const scalar dmidtfSign =
                vg.phase().name() == pair.first() ? +1 : -1;

            const sizeGroup& fi = popBal_.sizeGroups()[i];

            tmp<volScalarField> dDriftRate
            (
                dmidtfSign*dmidtf/(fi.phase().rho()*W_[k])
            );

            if (!numberWeighted_)
            {
                dDriftRate.ref() *= fi.a();
            }

            driftRate += dDriftRate;
        }
    }
}


// ************************************************************************* //
