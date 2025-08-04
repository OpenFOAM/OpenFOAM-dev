/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(populationBalance, 0);
    addToRunTimeSelectionTable(diameterModel, populationBalance, dictionary);

    // Backwards compatible lookup as "velocityGroup"
    addNamedToRunTimeSelectionTable
    (
        diameterModel,
        populationBalance,
        dictionary,
        velocityGroup
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalance::populationBalance
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    popBalName_(diameterProperties.lookup("populationBalance")),
    popBalPtr_(nullptr),
    nGroups_(diameterProperties.lookup<label>("nGroups")),
    iFirst_(-1),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase.time().name(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar("d", dimLength, NaN)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalance::~populationBalance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::populationBalanceModel&
Foam::diameterModels::populationBalance::popBal() const
{
    if (popBalPtr_ == nullptr)
    {
        popBalPtr_ =
            &phase().mesh().lookupObject<populationBalanceModel>(popBalName_);
    }

    return *popBalPtr_;
}


Foam::label Foam::diameterModels::populationBalance::iFirst() const
{
    if (iFirst_ == -1)
    {
        const UPtrList<const populationBalance>& uniqueDiameters =
            popBal().uniqueDiameters();

        // Determine the unique phase index within the population balance
        label uniquePhasei = -1;
        forAll(uniqueDiameters, uniquePhasej)
        {
            if (&uniqueDiameters[uniquePhasej] == this)
            {
                uniquePhasei = uniquePhasej;
                break;
            }
        }

        // If this is the first unique phase, then the first group of this
        // phase is the first of the entire population balance. Otherwise it is
        // the group immediately after those in the previous unique phase.
        iFirst_ =
            uniquePhasei == 0
          ? 0
          : uniqueDiameters[uniquePhasei - 1].iLast() + 1;
    }

    return iFirst_;
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalance::d() const
{
    return d_;
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalance::Av() const
{
    tmp<volScalarField> tsumFiAbyV
    (
        volScalarField::New
        (
            "sumFiAbyV",
            phase().mesh(),
            dimensionedScalar(inv(dimLength), Zero)
        )
    );
    volScalarField& sumFiAbyV = tsumFiAbyV.ref();

    const populationBalanceModel& popBal = this->popBal();

    for (label i = iFirst(); i <= iLast(); ++ i)
    {
        sumFiAbyV += popBal.f(i)*popBal.a(i)/popBal.v(i);
    }

    return phase()*tsumFiAbyV;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::diameterModels::populationBalance::fSum() const
{
    tmp<volScalarField::Internal> tsumFi =
        volScalarField::Internal::New
        (
            "sumFi",
            phase().mesh(),
            dimensionedScalar(dimless, 0)
        );
    volScalarField::Internal& sumFi = tsumFi.ref();

    const populationBalanceModel& popBal = this->popBal();

    for (label i = iFirst(); i <= iLast(); ++ i)
    {
        sumFi += popBal.f(i)();
    }

    return tsumFi;
}


void Foam::diameterModels::populationBalance::correct()
{
    const populationBalanceModel& popBal = this->popBal();

    if (popBal.solveOnFinalIterOnly() && !popBal.fluid().pimple().finalIter())
    {
        return;
    }

    tmp<volScalarField> tsumFi
    (
        volScalarField::New
        (
            "sumFi",
            phase().mesh(),
            dimensionedScalar(dimless, Zero)
        )
    );
    tmp<volScalarField> tsumFiAbyV
    (
        volScalarField::New
        (
            "sumFiAbyV",
            phase().mesh(),
            dimensionedScalar(inv(dimLength), Zero)
        )
    );
    volScalarField& sumFi = tsumFi.ref();
    volScalarField& sumFiAbyV = tsumFiAbyV.ref();

    for (label i = iFirst(); i <= iLast(); ++ i)
    {
        sumFi += max(popBal.f(i), rootVSmall);
        sumFiAbyV += max(popBal.f(i), rootVSmall)*popBal.a(i)/popBal.v(i);
    }

    d_ = 6*sumFi/tsumFiAbyV;

    Info<< phase().name() << " Sauter mean diameter, min, max = "
        << d_.weightedAverage(d_.mesh().V()).value()
        << ' ' << min(d_).value() << ' ' << max(d_).value() << endl;
}


bool Foam::diameterModels::populationBalance::read
(
    const dictionary& phaseProperties
)
{
    diameterModel::read(phaseProperties);

    return true;
}


// ************************************************************************* //
