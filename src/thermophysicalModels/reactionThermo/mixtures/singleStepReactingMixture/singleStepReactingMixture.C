/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "singleStepReactingMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void Foam::singleStepReactingMixture<ThermoType>::calculateqFuel()
{
    const Reaction<ThermoType>& reaction = this->operator[](0);
    const  scalar Wu = this->speciesData()[fuelIndex_].W();

    forAll(reaction.lhs(), i)
    {
        const label specieI = reaction.lhs()[i].index;
        const scalar stoichCoeff = reaction.lhs()[i].stoichCoeff;
        specieStoichCoeffs_[specieI] = -stoichCoeff;
        qFuel_.value() += this->speciesData()[specieI].hc()*stoichCoeff/Wu;
    }

    forAll(reaction.rhs(), i)
    {
        const label specieI = reaction.rhs()[i].index;
        const scalar stoichCoeff = reaction.rhs()[i].stoichCoeff;
        specieStoichCoeffs_[specieI] = stoichCoeff;
        qFuel_.value() -= this->speciesData()[specieI].hc()*stoichCoeff/Wu;
        specieProd_[specieI] = -1;
    }

    Info << "Fuel heat of combustion :" << qFuel_.value() << endl;
}


template<class ThermoType>
void Foam::singleStepReactingMixture<ThermoType>::massAndAirStoichRatios()
{
    const label O2Index = this->species()["O2"];
    const scalar Wu = this->speciesData()[fuelIndex_].W();

    stoicRatio_ =
       (this->speciesData()[inertIndex_].W()
      * specieStoichCoeffs_[inertIndex_]
      + this->speciesData()[O2Index].W()
      * mag(specieStoichCoeffs_[O2Index]))
      / (Wu*mag(specieStoichCoeffs_[fuelIndex_]));

    s_ =
        (this->speciesData()[O2Index].W()
      * mag(specieStoichCoeffs_[O2Index]))
      / (Wu*mag(specieStoichCoeffs_[fuelIndex_]));

    Info << "stoichiometric air-fuel ratio :" << stoicRatio_.value() << endl;

    Info << "stoichiometric oxygen-fuel ratio :" << s_.value() << endl;
}


template<class ThermoType>
void Foam::singleStepReactingMixture<ThermoType>::calculateMaxProducts()
{
    const Reaction<ThermoType>& reaction = this->operator[](0);

    scalar Wm = 0.0;
    scalar totalMol = 0.0;
    forAll(reaction.rhs(), i)
    {
        label specieI = reaction.rhs()[i].index;
        totalMol += mag(specieStoichCoeffs_[specieI]);
    }

    scalarList Xi(reaction.rhs().size());

    forAll(reaction.rhs(), i)
    {
        const label specieI = reaction.rhs()[i].index;
        Xi[i] = mag(specieStoichCoeffs_[specieI])/totalMol;

        Wm += Xi[i]*this->speciesData()[specieI].W();
    }

    forAll(reaction.rhs(), i)
    {
        const label specieI = reaction.rhs()[i].index;
        Yprod0_[specieI] =  this->speciesData()[specieI].W()/Wm*Xi[i];
    }

    Info << "Maximum products mass concentrations:" << nl;
    forAll(Yprod0_, i)
    {
        if (Yprod0_[i] > 0)
        {
            Info<< "    " << this->species()[i] << ": " << Yprod0_[i] << nl;
        }
    }

    // Normalize the stoichiometric coeff to mass
    forAll(specieStoichCoeffs_, i)
    {
        specieStoichCoeffs_[i] =
            specieStoichCoeffs_[i]
          * this->speciesData()[i].W()
          / (this->speciesData()[fuelIndex_].W()
          * mag(specieStoichCoeffs_[fuelIndex_]));
    }
}


template<class ThermoType>
void Foam::singleStepReactingMixture<ThermoType>::fresCorrect()
{
    const Reaction<ThermoType>& reaction = this->operator[](0);

    label O2Index = this->species()["O2"];
    const volScalarField& YFuel = this->Y()[fuelIndex_];
    const volScalarField& YO2 = this->Y()[O2Index];

    // reactants
    forAll(reaction.lhs(), i)
    {
        const label specieI = reaction.lhs()[i].index;
        if (specieI == fuelIndex_)
        {
            fres_[specieI] =  max(YFuel - YO2/s_, scalar(0.0));
        }
        else if (specieI == O2Index)
        {
            fres_[specieI] =  max(YO2 - YFuel*s_, scalar(0.0));
        }
    }


    // products
    forAll(reaction.rhs(), i)
    {
        const label specieI = reaction.rhs()[i].index;
        if (specieI != inertIndex_)
        {
            forAll(fres_[specieI], cellI)
            {
                if (fres_[fuelIndex_][cellI] > 0.0)
                {
                    // rich mixture
                    fres_[specieI][cellI] =
                        Yprod0_[specieI]
                      * (1.0 + YO2[cellI]/s_.value() - YFuel[cellI]);
                }
                else
                {
                    // lean mixture
                    fres_[specieI][cellI] =
                        Yprod0_[specieI]
                      * (
                            1.0
                          - YO2[cellI]/s_.value()*stoicRatio_.value()
                          + YFuel[cellI]*stoicRatio_.value()
                        );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::singleStepReactingMixture<ThermoType>::singleStepReactingMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
:
    reactingMixture<ThermoType>(thermoDict, mesh),
    stoicRatio_(dimensionedScalar("stoicRatio", dimless, 0.0)),
    s_(dimensionedScalar("s", dimless, 0.0)),
    qFuel_(dimensionedScalar("qFuel", sqr(dimVelocity), 0.0)),
    specieStoichCoeffs_(this->species_.size(), 0.0),
    Yprod0_(this->species_.size(), 0.0),
    fres_(Yprod0_.size()),
    inertIndex_(this->species()[thermoDict.lookup("inertSpecie")]),
    fuelIndex_(this->species()[thermoDict.lookup("fuel")]),
    specieProd_(Yprod0_.size(), 1)
{
    if (this->size() == 1)
    {
        forAll(fres_, fresI)
        {
            IOobject header
            (
                "fres_" + this->species()[fresI],
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            );

            fres_.set
            (
                fresI,
                new volScalarField
                (
                    header,
                    mesh,
                    dimensionedScalar("fres" + name(fresI), dimless, 0.0)
                )
            );
        }

        calculateqFuel();

        massAndAirStoichRatios();

        calculateMaxProducts();

        autoPtr<chemistryReader<ThermoType> >::clear();
    }
    else
    {
        FatalErrorIn
        (
            "singleStepReactingMixture::<ThermoType>::singleStepReactingMixture"
            "("
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Only one reaction required for single step reaction"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::singleStepReactingMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{}


// ************************************************************************* //
