/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "populationBalanceSetSizeDistribution.H"
#include "distribution.H"
#include "addToRunTimeSelectionTable.H"
#include "populationBalance.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(populationBalanceSetSizeDistribution, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        populationBalanceSetSizeDistribution,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::populationBalanceSetSizeDistribution::
populationBalanceSetSizeDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    popBalName_(dict.lookupOrDefault("populationBalance", word::null)),
    phaseName_(dict.lookupOrDefault("phase", word::null)),
    distribution_
    (
        distribution::New(dimLength, dict.subDict("distribution"), 3, -1)
    )
{
    const bool havePopBal = popBalName_ != word::null;
    const bool havePhase = phaseName_ != word::null;
    if (havePopBal == havePhase)
    {
        FatalIOErrorInFunction(dict)
            << (havePopBal ? "both" : "neither") << " of keywords "
            << "populationBalance " << (havePopBal ? "and" : "or")
            << " phase defined in dictionary " << dict.name()
            << exit(FatalIOError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::populationBalanceSetSizeDistribution::
~populationBalanceSetSizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::populationBalanceSetSizeDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::populationBalanceSetSizeDistribution::write()
{
    if (!functionObject::postProcess)
    {
        WarningInFunction
            << "The " << typeName << " function cannot be executed at run-time"
            << endl;

        return false;
    }

    Info<< nl << indent << name() << token::COLON << endl << incrIndent;

    if (popBalName_ != word::null)
    {
        populationBalanceModel& popBal =
            obr_.lookupObjectRef<populationBalanceModel>(popBalName_);

        HashTable<const diameterModels::populationBalance*> diameterPtrs;

        // Set the group fractions for the population balance
        forAll(popBal.fs(), i)
        {
            volScalarField& fi = popBal.f(i);

            fi == popBal.etaV(i, distribution_());

            Info<< indent << "Writing " << fi.name() << endl;

            fi.write();

            diameterPtrs.insert
            (
                popBal.phases()[i].name(),
                &popBal.diameters()[i]
            );
        }

        // Set the volume fractions if there are multiple phases in this
        // population balance
        if (diameterPtrs.size() > 1)
        {
            tmp<volScalarField> talphas =
                volScalarField::New
                (
                    IOobject::groupName("alpha", popBal.name()),
                    mesh(),
                    dimensionedScalar(dimless, scalar(0))
                );
            volScalarField& alphas = talphas.ref();

            forAllConstIter
            (
                HashTable<const diameterModels::populationBalance*>,
                diameterPtrs,
                iter
            )
            {
                alphas += iter()->phase();
            }

            forAllConstIter
            (
                HashTable<const diameterModels::populationBalance*>,
                diameterPtrs,
                iter
            )
            {
                volScalarField& alpha =
                    const_cast<phaseModel&>(iter()->phase());

                const label i0 = iter()->iFirst();
                const label i1 = iter()->iLast();

                alpha == popBal.etaV(labelPair(i0, i1), distribution_)*alphas;

                Info<< indent << "Writing " << alpha.name() << endl;

                alpha.write();
            }
        }
    }
    else
    {
        diameterModels::populationBalance& diameter =
            refCast<diameterModels::populationBalance>
            (
                obr_.lookupObjectRef<phaseSystem>
                (
                    phaseSystem::propertiesName
                ).phases()[phaseName_].diameter()
            );

        populationBalanceModel& popBal =
            obr_.lookupObjectRef<populationBalanceModel>(diameter.popBalName());

        // Set the group fractions for this phase
        for (label i = diameter.iFirst(); i <= diameter.iLast(); ++ i)
        {
            volScalarField& fi = popBal.f(i);

            fi == popBal.etaV(i, distribution_());

            Info<< indent << "Writing " << fi.name() << endl;

            fi.write();
        }
    }

    Info<< decrIndent;

    return true;
}


// ************************************************************************* //
