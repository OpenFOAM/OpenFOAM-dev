/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "turbulentDispersion.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToIncompressibleFluid.H"
#include "coupledToFluid.H"
#include "standardNormal.H"
#include "maxLagrangianFieldSources.H"
#include "NaNLagrangianFieldSources.H"
#include "internalLagrangianFieldSources.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(turbulentDispersion, 0);
    addToRunTimeSelectionTable
    (
        LagrangianModel,
        turbulentDispersion,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::Tuple2<bool, Foam::CloudStateField<Type>&>
Foam::Lagrangian::turbulentDispersion::initialiseTurbField
(
    const word& name,
    const dimensionSet& dims,
    const Type& value
)
{
    CloudStateField<Type>& ref =
        cloud().stateField<Type>
        (
            IOobject
            (
                name,
                mesh().time().name(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensioned<Type>(name, dims, value)
        );

    return Tuple2<bool, CloudStateField<Type>&>(ref.headerOk(), ref);
}


template<class InjectionFieldSourceType, class Type>
void Foam::Lagrangian::turbulentDispersion::completeTurbField
(
    Tuple2<bool, CloudStateField<Type>&>& turbField
)
{
    if (turbField.first()) return;

    turbField.first() = true;

    LagrangianModels& modelList = cloud().LagrangianModels();

    turbField.second().sourcesRef().table().transfer
    (
        typename LagrangianDynamicField<Type>::Sources
        (
            turbField.second(),
            modelList.modelTypeFieldSourceTypes
            <
                LagrangianInjection,
                InjectionFieldSourceType,
                LagrangianSource,
                internalLagrangianFieldSource<Type>
            >()
        ).table()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::turbulentDispersion::turbulentDispersion
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianModel(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    dragPtr_(nullptr),
    momentumTransportModel_(mesh.mesh().lookupType<momentumTransportModel>()),
    kc_
    (
        cloud<clouds::coupled>().carrierField<scalar>
        (
            clouds::coupled::carrierName
            (
                momentumTransportModel_.k()().name()
            ),
            [&]()
            {
                return momentumTransportModel_.k();
            }
        )
    ),
    epsilonc_
    (
        cloud<clouds::coupled>().carrierField<scalar>
        (
            clouds::coupled::carrierName
            (
                momentumTransportModel_.epsilon()().name()
            ),
            [&]()
            {
                return momentumTransportModel_.epsilon();
            }
        )
    ),
    fractionTurb_(initialiseTurbField("fractionTurb", dimless, vGreat)),
    tTurb_(initialiseTurbField("tTurb", dimTime, NaN)),
    Uturb_(initialiseTurbField("Uturb", dimVelocity, vector::uniform(NaN))),
    rndGen_("rndGen", stateDict, name, false),
    avgUturbPtr_(nullptr)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::Lagrangian::turbulentDispersion::addSupFields() const
{
    return wordList(1, cloud().U.name());
}


void Foam::Lagrangian::turbulentDispersion::postConstruct()
{
    LagrangianModels& modelList = cloud().LagrangianModels();

    dragPtr_ = nullptr;

    forAll(modelList, i)
    {
        if (!isA<drag>(modelList[i])) continue;

        if (dragPtr_ != nullptr)
        {
            FatalErrorInFunction
                << "Multiple drag models found. Turbulent dispersion "
                << "requires exactly one drag model."
                << exit(FatalError);
        }

        dragPtr_ = &refCast<const drag>(modelList[i]);
    }

    if (dragPtr_ == nullptr)
    {
        FatalErrorInFunction
            << "No drag models found. Turbulent dispersion "
            << "requires exactly one drag model."
            << exit(FatalError);
    }

    completeTurbField<maxLagrangianScalarFieldSource>(fractionTurb_);
    completeTurbField<NaNLagrangianScalarFieldSource>(tTurb_);
    completeTurbField<NaNLagrangianVectorFieldSource>(Uturb_);
}


void Foam::Lagrangian::turbulentDispersion::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    // References to the evolving fields
    LagrangianSubScalarSubField& fractionTurb =
        fractionTurb_.second().ref(subMesh);
    LagrangianSubScalarSubField& tTurb = tTurb_.second().ref(subMesh);
    LagrangianSubVectorSubField& Uturb = Uturb_.second().ref(subMesh);
    const LagrangianSubScalarField& fractionTurb0 = fractionTurb.oldTime();
    const LagrangianSubVectorField& Uturb0 = Uturb.oldTime();

    // Update the eddy time-scale
    const dimensionedScalar cps(dimless, 0.16432); // (model constant ?)
    const LagrangianSubScalarField magUrel
    (
        mag(cloud().U(subMesh) - cloud<clouds::coupled>().Uc(subMesh))
    );
    static const dimensionedScalar rootVSmallEpsilon
    (
        dimEnergy/dimMass/dimTime,
        rootVSmall
    );
    static const dimensionedScalar rootVSmallEpsilonU
    (
        dimEnergy/dimMass/dimTime*dimVelocity,
        rootVSmall
    );
    tTurb =
        max
        (
            kc_(subMesh)/max(epsilonc_(subMesh), rootVSmallEpsilon),
            cps*kc_(subMesh)*sqrt(kc_(subMesh))
           /max(epsilonc_(subMesh)*magUrel, rootVSmallEpsilonU)
        );

    // Velocity fluctuation magnitude
    const LagrangianSubScalarField magUturb(sqrt(2.0/3.0*kc_(subMesh)));

    // Initialise the average turbulent velocities
    avgUturbPtr_.reset
    (
        LagrangianSubVectorField::New
        (
            "avgUturb",
            subMesh,
            dimensionedVector(dimVelocity, Zero)
        ).ptr()
    );

    // Consider each particle in turn
    forAll(subMesh, subi)
    {
        // Create independent generators for each particle. That way if the
        // sequence grows or shrinks across iterations, that effect doesn't
        // propagate and affect other particles' sequences
        randomGenerator rndGen
        (
            rndGen_.sampleAB(labelMin, labelMax),
            false
        );
        distributions::standardNormal stdNormal
        (
            rndGen_.sampleAB(labelMin, labelMax),
            false
        );

        // Set up sub-stepping
        const scalar Dt = max(deltaT[subi], rootVSmall);
        scalar dt = 0;

        // Continue/complete the previous eddy
        if (fractionTurb0[subi] < 1)
        {
            dt += tTurb[subi]*(1 - fractionTurb0[subi]);

            avgUturbPtr_()[subi] += min(dt/Dt, 1)*Uturb0[subi];
        }

        // Add new eddies across the time-step
        while (dt < Dt)
        {
            const scalar dtPrev = dt;
            dt += tTurb[subi];

            // Create a random direction with normally distributed magnitude
            const scalar theta =
                rndGen.scalar01()*constant::mathematical::twoPi;
            const scalar u = 2*rndGen.scalar01() - 1;
            const scalar a = sqrt(1 - sqr(u));
            const vector dir(a*cos(theta), a*sin(theta), u);

            // Set the new turbulent fluctuation velocity
            Uturb[subi] = dir*stdNormal.sample()*magUturb[subi];

            // Add it to the average
            avgUturbPtr_()[subi] += (min(dt/Dt, 1) - dtPrev/Dt)*Uturb[subi];
        }

        // Set the fraction to where we got to in the current eddy
        fractionTurb[subi] = 1 - (dt - Dt)/tTurb[subi];
    }

    // If this is not the final iteration then rewind the generator so that the
    // same numbers are generated in the next iteration
    rndGen_.start(!final);
}


void Foam::Lagrangian::turbulentDispersion::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    assertCloud<clouds::coupledToIncompressibleFluid>();

    eqn.Su += dragPtr_->D(deltaT.mesh())*avgUturbPtr_();
}


void Foam::Lagrangian::turbulentDispersion::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    assertCloud<clouds::coupledToIncompressibleFluid, clouds::coupledToFluid>();

    eqn.Su += dragPtr_->D(deltaT.mesh())*avgUturbPtr_();
}


void Foam::Lagrangian::turbulentDispersion::writeProcessorState
(
    Ostream& os
) const
{
    LagrangianModel::writeProcessorState(os);

    writeEntry(os, "rndGen", rndGen_);
}


// ************************************************************************* //
