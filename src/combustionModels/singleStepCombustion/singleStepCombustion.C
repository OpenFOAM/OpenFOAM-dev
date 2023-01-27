/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "fvSolution.H"
#include "singleStepCombustion.H"
#include "fvmSup.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::combustionModels::singleStepCombustion::calculateqFuel()
{
    const scalar Wu = mixture_.Wi(fuelIndex_);

    forAll(reaction_.lhs(), i)
    {
        const label speciei = reaction_.lhs()[i].index;
        const scalar stoichCoeff = reaction_.lhs()[i].stoichCoeff;
        specieStoichCoeffs_[speciei] = -stoichCoeff;
        qFuel_.value() +=
            mixture_.Hf(speciei)*mixture_.Wi(speciei)*stoichCoeff/Wu;
    }

    forAll(reaction_.rhs(), i)
    {
        const label speciei = reaction_.rhs()[i].index;
        const scalar stoichCoeff = reaction_.rhs()[i].stoichCoeff;
        specieStoichCoeffs_[speciei] = stoichCoeff;
        qFuel_.value() -=
            mixture_.Hf(speciei)*mixture_.Wi(speciei)*stoichCoeff/Wu;
        specieProd_[speciei] = -1;
    }

    Info << "Fuel heat of combustion: " << qFuel_.value() << endl;
}


void Foam::combustionModels::singleStepCombustion::massAndAirStoichRatios()
{
    const label O2Index = mixture_.species()["O2"];
    const scalar Wu = mixture_.Wi(fuelIndex_);

    stoicRatio_ =
        (
            mixture_.Wi(mixture_.defaultSpecie())
           *specieStoichCoeffs_[mixture_.defaultSpecie()]
          + mixture_.Wi(O2Index)*mag(specieStoichCoeffs_[O2Index])
        )/(Wu*mag(specieStoichCoeffs_[fuelIndex_]));

    s_ = mixture_.Wi(O2Index)*mag(specieStoichCoeffs_[O2Index])
        /(Wu*mag(specieStoichCoeffs_[fuelIndex_]));

    Info << "stoichiometric air-fuel ratio: " << stoicRatio_.value() << endl;
    Info << "stoichiometric oxygen-fuel ratio: " << s_.value() << endl;
}


void Foam::combustionModels::singleStepCombustion::calculateMaxProducts()
{
    scalar Wm = 0.0;
    scalar totalMol = 0.0;
    forAll(reaction_.rhs(), i)
    {
        label speciei = reaction_.rhs()[i].index;
        totalMol += mag(specieStoichCoeffs_[speciei]);
    }

    scalarList Xi(reaction_.rhs().size());

    forAll(reaction_.rhs(), i)
    {
        const label speciei = reaction_.rhs()[i].index;
        Xi[i] = mag(specieStoichCoeffs_[speciei])/totalMol;
        Wm += Xi[i]*mixture_.Wi(speciei);
    }

    forAll(reaction_.rhs(), i)
    {
        const label speciei = reaction_.rhs()[i].index;
        Yprod0_[speciei] =  mixture_.Wi(speciei)/Wm*Xi[i];
    }

    Info << "Maximum products mass concentrations: " << nl;
    forAll(Yprod0_, i)
    {
        if (Yprod0_[i] > 0)
        {
            Info<< "    " << mixture_.species()[i] << ": " << Yprod0_[i] << nl;
        }
    }

    // Normalise the stoichiometric coeff to mass
    forAll(specieStoichCoeffs_, i)
    {
        specieStoichCoeffs_[i] =
            specieStoichCoeffs_[i]*mixture_.Wi(i)
           /(mixture_.Wi(fuelIndex_)*mag(specieStoichCoeffs_[fuelIndex_]));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::singleStepCombustion::singleStepCombustion
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    combustionModel(modelType, thermo, turb, combustionProperties),
    mixture_(dynamic_cast<const basicSpecieMixture&>(this->thermo())),
    reaction_(mixture_.species(), this->subDict("reaction")),
    stoicRatio_(dimensionedScalar("stoicRatio", dimless, 0)),
    s_(dimensionedScalar("s", dimless, 0)),
    qFuel_(dimensionedScalar("qFuel", sqr(dimVelocity), 0)),
    specieStoichCoeffs_(mixture_.species().size(), 0.0),
    Yprod0_(mixture_.species().size(), 0.0),
    fres_(Yprod0_.size()),
    fuelIndex_(mixture_.species()[thermo.properties().lookup("fuel")]),
    specieProd_(Yprod0_.size(), 1),
    wFuel_
    (
        IOobject
        (
            this->thermo().phasePropertyName("wFuel"),
            this->mesh().time().name(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    ),
    semiImplicit_(readBool(this->coeffs_.lookup("semiImplicit")))
{
    forAll(fres_, fresI)
    {
        IOobject header
        (
            "fres_" + mixture_.species()[fresI],
            this->mesh().time().name(),
            this->mesh()
        );

        fres_.set
        (
            fresI,
            new volScalarField
            (
                header,
                this->mesh(),
                dimensionedScalar("fres" + Foam::name(fresI), dimless, 0)
            )
        );
    }

    calculateqFuel();

    massAndAirStoichRatios();

    calculateMaxProducts();

    if (semiImplicit_)
    {
        Info<< "Combustion mode: semi-implicit" << endl;
    }
    else
    {
        Info<< "Combustion mode: explicit" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::singleStepCombustion::~singleStepCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::combustionModels::singleStepCombustion::R(const label speciei) const
{
    return wFuel_()*specieStoichCoeffs()[speciei];
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::singleStepCombustion::R(volScalarField& Y) const
{
    const label specieI = mixture_.species()[Y.member()];

    volScalarField wSpecie
    (
        wFuel_*specieStoichCoeffs()[specieI]
    );

    if (semiImplicit_)
    {
        const label fNorm = specieProd()[specieI];
        const volScalarField fres(this->fres(specieI));
        wSpecie /= max(fNorm*(Y - fres), scalar(1e-2));

        return -fNorm*wSpecie*fres + fNorm*fvm::Sp(wSpecie, Y);
    }
    else
    {
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::singleStepCombustion::Qdot() const
{
    const label fuelI = fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermo().composition().Y(fuelI));

    return -qFuel()*(R(YFuel) & YFuel);
}


bool Foam::combustionModels::singleStepCombustion::read()
{
    if (combustionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::combustionModels::singleStepCombustion::fresCorrect()
{
    const label O2Index = mixture_.species()["O2"];
    const volScalarField& YFuel = mixture_.Y()[fuelIndex_];
    const volScalarField& YO2 = mixture_.Y()[O2Index];

    // reactants
    forAll(reaction_.lhs(), i)
    {
        const label speciei = reaction_.lhs()[i].index;
        if (speciei == fuelIndex_)
        {
            fres_[speciei] = max(YFuel - YO2/s_, scalar(0));
        }
        else if (speciei == O2Index)
        {
            fres_[speciei] = max(YO2 - YFuel*s_, scalar(0));
        }
    }

    // products
    forAll(reaction_.rhs(), i)
    {
        const label speciei = reaction_.rhs()[i].index;
        if (speciei != mixture_.defaultSpecie())
        {
            forAll(fres_[speciei], celli)
            {
                if (fres_[fuelIndex_][celli] > 0.0)
                {
                    // rich mixture
                    fres_[speciei][celli] =
                        Yprod0_[speciei]
                      * (1.0 + YO2[celli]/s_.value() - YFuel[celli]);
                }
                else
                {
                    // lean mixture
                    fres_[speciei][celli] =
                        Yprod0_[speciei]
                      * (
                            1.0
                          - YO2[celli]/s_.value()*stoicRatio_.value()
                          + YFuel[celli]*stoicRatio_.value()
                        );
                }
            }
        }
    }
}


// ************************************************************************* //
