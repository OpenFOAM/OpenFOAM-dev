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

#include "collisionPhaseTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "carried.H"
#include "grouped.H"
#include "massive.H"
#include "standardNormal.H"
#include "LagrangianmSp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(collisionPhaseTransfer, 0);
    addToRunTimeSelectionTable
    (
        LagrangianModel,
        collisionPhaseTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Lagrangian::collisionPhaseTransfer::addSupType
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    LagrangianEqn<scalar>& eqn
) const
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();
    const SubList<LagrangianState> subStates = subMesh.sub(mesh().states());

    if (eqn.isPsi(vOrM))
    {
        // Add the implicit coefficient to the models' total
        if (sumDeltaTSp_.valid()) sumDeltaTSp_.ref() += PhitPtr_();

        // Add an implicit source to the cloud
        eqn.deltaTSp -= PhitPtr_();
    }
    else
    {
        // Calculate the rate of volume/mass consumption
        tmp<LagrangianSubScalarField> deltaTSu = PhitPtr_()*vOrM;

        // Modify the rate of consumption to fully consume removed particles
        forAll(subMesh, subi)
        {
            if (subStates[subi] == LagrangianState::toBeRemoved)
            {
                deltaTSu.ref()[subi] *= 1 + 1/max(sumDeltaTSp_()[subi], small);
            }
        }

        // Add an explicit source to the carrier
        eqn.deltaTSu -= deltaTSu;
    }
}


template<class Type>
void Foam::Lagrangian::collisionPhaseTransfer::addSupType
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubSubField<Type>& field,
    LagrangianEqn<Type>& eqn
) const
{
    FatalErrorInFunction
        << "The Lagrangian model '" << name() << "' of cloud '"
        << mesh().name() << "' cannot add a transfer for the field '"
        << field.name() << "' because the equation is not volume or "
        << "mass-weighted" << exit(FatalError);
}


template<class Type>
void Foam::Lagrangian::collisionPhaseTransfer::addSupType
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubSubField<Type>& field,
    LagrangianEqn<Type>& eqn
) const
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();
    const SubList<LagrangianState> subStates = subMesh.sub(mesh().states());

    if (eqn.isPsi(field))
    {
        // Add an implicit source to the cloud
        eqn.deltaTSp -= PhitPtr_()*vOrM;
    }
    else
    {
        // Calculate the rate of consumption
        tmp<LagrangianSubField<Type>> deltaTSu = PhitPtr_()*vOrM*field;

        // Modify the rate of consumption to fully consume removed particles
        forAll(subMesh, subi)
        {
            if (subStates[subi] == LagrangianState::toBeRemoved)
            {
                deltaTSu.ref()[subi] *= 1 + 1/max(sumDeltaTSp_()[subi], small);
            }
        }

        // Add an explicit source to the carrier
        eqn.deltaTSu -= deltaTSu;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::collisionPhaseTransfer::collisionPhaseTransfer
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianSource(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    alphac_
    (
        cloud<clouds::carried>().carrierField
        (
            mesh.mesh().lookupObject<volScalarField>
            (
                IOobject::groupName
                (
                    "alpha",
                    cloud<clouds::carried>().carrierPhaseName()
                )
            )
        )
    ),
    rndGen_("rndGen", stateDict, name, false),
    PhitPtr_(nullptr),
    sumDeltaTSp_()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList
Foam::Lagrangian::collisionPhaseTransfer::addSupFields() const
{
    return
        isCloud<clouds::massive>()
      ? wordList({word::null, clouds::massive::mName})
      : wordList({word::null, clouds::shaped::vName});
}


bool Foam::Lagrangian::collisionPhaseTransfer::addsSupToField
(
    const word& fieldName,
    const word& eqnFieldName
) const
{
    const word eqnPhaseName = IOobject::group(eqnFieldName);
    const word eqnMember = IOobject::member(eqnFieldName);

    const bool isCarrierEqn = eqnMember[eqnMember.size() - 1] == 'c';

    return
        !isCarrierEqn
     || (
            cloud<clouds::carried>().hasPhase()
         && (
                eqnPhaseName == word::null
             || eqnPhaseName == cloud<clouds::carried>().phaseName()
            )
        );
}


void Foam::Lagrangian::collisionPhaseTransfer::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();
    SubList<LagrangianState> subStates = subMesh.sub(mesh().states());

    auto sampleUniform01 = [&]()
    {
        return
            LagrangianSubScalarField::New
            (
                "sampleUniform01",
                subMesh,
                dimless,
                rndGen_.sample01<scalar>(subMesh.size())
            );
    };

    auto sampleStandardNormal = [&]()
    {
        return
            LagrangianSubScalarField::New
            (
                "sampleStandardNormal",
                subMesh,
                dimless,
                distributions::standardNormal
                (
                    rndGen_.sampleAB(labelMin, labelMax),
                    false
                ).sample(subMesh.size())
            );
    };

    // Length scale for the containing cells
    tmp<LagrangianSubScalarField> L =
        LagrangianSubScalarField::New
        (
            "L",
            subMesh,
            dimLength,
            cbrt
            (
                scalarField
                (
                    subMesh.mesh().mesh().cellVolumes(),
                    subMesh.sub(subMesh.mesh().celli())
                )
            )
        );

    // This is a simple model for the probability that a particle hits the
    // interface during this step
    PhitPtr_.set
    (
        (
            max(mag(alphac_.grad(subMesh)), rootVSmall/L)
           *mag(cloud().U(subMesh))
           /max(1 - alphac_(subMesh), rootVSmall)
           *deltaT
        ).ptr()
    );
    LagrangianSubScalarField& Phit = PhitPtr_();

    // Sample to randomise the collisions
    if (!isCloud<clouds::grouped>())
    {
        Phit = pos(Phit - sampleUniform01());
    }
    else
    {
        const LagrangianSubScalarSubField number
        (
            cloud<clouds::grouped>().number(subMesh)
        );

        Phit.maxMin(scalar(0), scalar(1));

        Phit += sqrt(Phit*(1 - Phit)/number)*sampleStandardNormal();

        Phit.maxMin(scalar(0), scalar(1));
    }

    // If this is not the final iteration then rewind the generator so that the
    // same numbers are generated in the next iteration
    rndGen_.start(!final);

    if (!final) return;

    // Mark any particles with high probabilities as to be removed
    forAll(subMesh, subi)
    {
        if (Phit[subi] == 1)
        {
            subStates[subi] = LagrangianState::toBeRemoved;
        }
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::Lagrangian::collisionPhaseTransfer::source
(
    const word& fieldName,
    const LagrangianSubMesh& subMesh
)
const
{
    NotImplemented;
    return tmp<LagrangianSubScalarField>(nullptr);
}


void Foam::Lagrangian::collisionPhaseTransfer::preAddSup
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    if (!final) return;

    const LagrangianSubMesh& subMesh = deltaT.mesh();
    const SubList<LagrangianState> subStates = subMesh.sub(mesh().states());

    // If particles are being removed then create/lookup the sum of the
    // implicit coefficients
    if (findIndex(subStates, LagrangianState::toBeRemoved) != -1)
    {
        sumDeltaTSp_.set
        (
            IOobject
            (
                "sumDeltaTSpvOrM",
                mesh().time().name(),
                mesh()
            ),
            subMesh,
            dimensionedScalar(dimless, scalar(0))
        );
    }
}


void Foam::Lagrangian::collisionPhaseTransfer::addSup
(
    const LagrangianSubScalarField& deltaT,
    LagrangianEqn<scalar>& eqn
) const
{
    eqn.deltaTSp -= PhitPtr_();
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_LAGRANGIAN_MODEL_ADD_FIELD_SUP,
    Lagrangian::collisionPhaseTransfer
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_LAGRANGIAN_MODEL_ADD_V_OR_M_FIELD_SUP,
    Lagrangian::collisionPhaseTransfer
)


void Foam::Lagrangian::collisionPhaseTransfer::postAddSup
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    PhitPtr_.clear();

    if (final) sumDeltaTSp_.clear();
}


void Foam::Lagrangian::collisionPhaseTransfer::writeProcessorState
(
    Ostream& os
) const
{
    LagrangianModel::writeProcessorState(os);

    writeEntry(os, "rndGen", rndGen_);
}


// ************************************************************************* //
