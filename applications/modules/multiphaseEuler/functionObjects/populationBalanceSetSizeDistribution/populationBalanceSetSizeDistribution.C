/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2026 OpenFOAM Foundation
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
#include "fvCellZoneIndexBoundaryField.H"
#include "distribution.H"
#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"

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


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class EtaVArg>
void Foam::functionObjects::populationBalanceSetSizeDistribution::setField
(
    const populationBalanceModel& popBal,
    const EtaVArg& etaVArg,
    volScalarField& f
)
{
    const fvCellZoneIndexBoundaryField& zoneIndexBf =
        fvCellZoneIndexBoundaryField::New(mesh(), zoneIndices_);

    f == popBal.etaV(etaVArg, defaultDistribution_());

    PtrList<dimensionedScalar> zoneEtaVs(zoneIndices_.size());
    forAll(zoneIndices_, zonei)
    {
        zoneEtaVs.set
        (
            zonei,
            new dimensionedScalar
            (
                popBal.etaV(etaVArg, zoneDistributions_[zonei])
            )
        );
    }

    forAll(zoneIndices_, zonei)
    {
        const cellZone& cells = mesh().cellZones()[zoneIndices_[zonei]];

        UIndirectList<scalar>(f.primitiveFieldRef(), cells) =
            zoneEtaVs[zonei].value();
    }

    volScalarField::BoundaryField& fBf = f.boundaryFieldRef();

    forAll(mesh().boundary(), patchi)
    {
        forAll(mesh().boundary()[patchi], patchFacei)
        {
            const label zonei = zoneIndexBf[patchi][patchFacei];

            if (zonei != -1)
            {
                fBf[patchi][patchFacei] = zoneEtaVs[zonei].value();
            }
        }
    }
}


void Foam::functionObjects::populationBalanceSetSizeDistribution::setFields()
{
    if (popBalName_ != word::null)
    {
        populationBalanceModel& popBal =
            obr_.lookupObjectRef<populationBalanceModel>(popBalName_);

        // Set the group fractions for the population balance
        forAll(popBal.fs(), i)
        {
            setField(popBal, i, popBal.f(i));
        }

        // Set the volume fractions if there are multiple phases in this
        // population balance. But retain the same total amount of phase across
        // the population balance as a whole.
        if (popBal.uniquePhases().size() > 1)
        {
            tmp<volScalarField> talphas =
                volScalarField::New
                (
                    IOobject::groupName("alpha", popBal.name()),
                    mesh(),
                    dimensionedScalar(dimless, scalar(0))
                );
            volScalarField& alphas = talphas.ref();
            forAll(popBal.uniquePhases(), uniquePhasei)
            {
                alphas += popBal.uniquePhases()[uniquePhasei];
            }

            forAll(popBal.uniquePhases(), uniquePhasei)
            {
                const phaseModel& phase = popBal.uniquePhases()[uniquePhasei];

                volScalarField& alpha = const_cast<phaseModel&>(phase);

                setField
                (
                    popBal,
                    labelPair
                    (
                        popBal.uniqueDiameters()[uniquePhasei].iFirst(),
                        popBal.uniqueDiameters()[uniquePhasei].iLast()
                    ),
                    alpha
                );

                alpha *= alphas;
            }
        }

        popBal.info();
    }
    else
    {
        const diameterModels::populationBalance& diameter =
            refCast<const diameterModels::populationBalance>
            (
                obr_.lookupObject<phaseSystem>
                (
                    phaseSystem::propertiesName
                ).phases()[phaseName_].diameter()
            );

        populationBalanceModel& popBal =
            const_cast<populationBalanceModel&>(diameter.popBal());

        // Set the group fractions for this phase
        for (label i = diameter.iFirst(); i <= diameter.iLast(); ++ i)
        {
            setField(popBal, i, popBal.f(i));
        }

        popBal.info();
    }
}


void Foam::functionObjects::populationBalanceSetSizeDistribution::writeFields()
{
    if (popBalName_ != word::null)
    {
        const populationBalanceModel& popBal =
            obr_.lookupObject<populationBalanceModel>(popBalName_);

        // Write the group fractions for the population balance
        forAll(popBal.fs(), i)
        {
            const volScalarField& fi = popBal.f(i);

            Info<< indent << "Writing " << fi.name() << endl;

            fi.write();
        }

        // Write the volume fractions if there are multiple phases in this
        // population balance
        if (popBal.uniquePhases().size() > 1)
        {
            forAll(popBal.uniquePhases(), uniquePhasei)
            {
                const volScalarField& alpha =
                    popBal.uniquePhases()[uniquePhasei];

                Info<< indent << "Writing " << alpha.name() << endl;

                alpha.write();
            }
        }
    }
    else
    {
        const diameterModels::populationBalance& diameter =
            refCast<const diameterModels::populationBalance>
            (
                obr_.lookupObject<phaseSystem>
                (
                    phaseSystem::propertiesName
                ).phases()[phaseName_].diameter()
            );

        const populationBalanceModel& popBal = diameter.popBal();

        // Write the group fractions for this phase
        for (label i = diameter.iFirst(); i <= diameter.iLast(); ++ i)
        {
            const volScalarField& fi = popBal.f(i);

            Info<< indent << "Writing " << fi.name() << endl;

            fi.write();
        }
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
    phaseName_(dict.lookupOrDefault("phase", word::null))
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

    const bool haveDistribution = dict.found("distribution");
    const bool haveDefaultOrZones =
        dict.found("default") || dict.found("zones");
    if (haveDistribution == haveDefaultOrZones)
    {
        FatalIOErrorInFunction(dict)
            << (haveDistribution ? "both" : "neither") << " of keywords "
            << "distribution " << (haveDistribution ? "and" : "or")
            << "default/zones defined in dictionary " << dict.name()
            << exit(FatalIOError);
    }

    defaultDistribution_.set
    (
        distribution::New
        (
            dimensions::length,
            dict.subDict(haveDistribution ? "distribution" : "default"),
            3,
            -1
        ).ptr()
    );

    if (haveDefaultOrZones)
    {
        const dictionary& zonesDict = dict.subDict("zones");

        forAllConstIter(dictionary, zonesDict, iter)
        {
            zoneIndices_.append(mesh().cellZones().findIndex(iter().keyword()));
            zoneDistributions_.append
            (
                distribution::New
                (
                    dimensions::length,
                    iter().dict(),
                    3,
                    -1
                ).ptr()
            );
        }
    }

    if (!functionObject::postProcess && mesh().time().timeIndex() == 0)
    {
        setFields();
    }
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
    if (functionObject::postProcess)
    {
        Info<< nl << indent << name() << token::COLON << endl << incrIndent;

        setFields();

        writeFields();

        Info<< decrIndent;
    }

    return true;
}


// ************************************************************************* //
