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

#include "LagrangianSubFieldsFwd.H"
#include "constantFluxCarrierTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToConstantDensityFluid.H"
#include "massive.H"
#include "shaped.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(constantFluxCarrierTransfer, 0);
    addToRunTimeSelectionTable
    (
        LagrangianModel,
        constantFluxCarrierTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Lagrangian::constantFluxCarrierTransfer::readCoeffs
(
    const dictionary& modelDict
)
{
    const bool haveVolumeFlux = modelDict.found("volumeFlux");
    const bool haveMassFlux = modelDict.found("massFlux");

    if (haveVolumeFlux == haveMassFlux)
    {
        FatalIOErrorInFunction(modelDict)
            << "keywords volumeFlux and massFlux are both "
            << (haveVolumeFlux ? "" : "un" ) << "defined in "
            << "dictionary " << modelDict.name()
            << exit(FatalIOError);
    }

    flux_.read(modelDict);

    if (flux_.value() < 0)
    {
        FatalIOErrorInFunction(modelDict)
            << flux_.name() << " must be positive"
            << exit(FatalError);
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::Lagrangian::constantFluxCarrierTransfer::Sp
(
    const LagrangianSubMesh& subMesh
) const
{
    tmp<LagrangianSubScalarField> tSByV =
        cloud<clouds::shaped>().aByV(subMesh)*flux_;

    if (flux_.dimensions() == dimMass/dimArea/dimTime)
    {
        assertCloud
        <
            clouds::coupledToConstantDensityFluid,
            clouds::massive
        >();

        if (isCloud<clouds::coupledToConstantDensityFluid>())
        {
            return
                tSByV
               /cloud<clouds::coupledToConstantDensityFluid>().rho();
        }
        else
        {
            return
                tSByV
               /cloud<clouds::massive>().rho(subMesh);
        }
    }
    else
    {
        return tSByV;
    }
}


void Foam::Lagrangian::constantFluxCarrierTransfer::addSupType
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    LagrangianEqn<scalar>& eqn
) const
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();
    const SubList<LagrangianState> subStates = subMesh.sub(mesh().states());

    const LagrangianSubScalarField Sp(this->Sp(subMesh));

    if (eqn.isPsi(vOrM))
    {
        // Add the implicit coefficient to the models' total
        if (sumDeltaTSp_.valid()) sumDeltaTSp_.ref() += deltaT*Sp;

        // Add an implicit source to the cloud
        eqn.Sp -= Sp;
    }
    else
    {
        // Calculate the rate of consumption
        tmp<LagrangianSubScalarField> Su = Sp*vOrM;

        // Modify the rate of consumption to fully consume removed particles
        forAll(subMesh, subi)
        {
            if (subStates[subi] == LagrangianState::toBeRemoved)
            {
                Su.ref()[subi] *= 1 + 1/max(sumDeltaTSp_()[subi], small);
            }
        }

        // Add an explicit source to the carrier
        eqn.Su -= Su;
    }
}


template<class Type>
void Foam::Lagrangian::constantFluxCarrierTransfer::addSupType
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
void Foam::Lagrangian::constantFluxCarrierTransfer::addSupType
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubSubField<Type>& field,
    LagrangianEqn<Type>& eqn
) const
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();
    const SubList<LagrangianState> subStates = subMesh.sub(mesh().states());

    const LagrangianSubScalarField Sp(this->Sp(subMesh));

    if (eqn.isPsi(field))
    {
        // Add an implicit source to the cloud
        eqn.Sp -= Sp*vOrM;
    }
    else
    {
        // Calculate the rate of consumption
        tmp<LagrangianSubField<Type>> Su = Sp*vOrM*field;

        // Modify the rate of consumption to fully consume removed particles
        forAll(subMesh, subi)
        {
            if (subStates[subi] == LagrangianState::toBeRemoved)
            {
                Su.ref()[subi] *= 1 + 1/max(sumDeltaTSp_()[subi], small);
            }
        }

        // Add an explicit source to the carrier
        eqn.Su -= Su;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::constantFluxCarrierTransfer::constantFluxCarrierTransfer
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianSource(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    flux_
    (
        modelDict.found("volumeFlux")
      ? dimensionedScalar("volumeFlux", dimVolume/dimArea/dimTime, modelDict)
      : modelDict.found("massFlux")
      ? dimensionedScalar("massFlux", dimMass/dimArea/dimTime, modelDict)
      : dimensionedScalar()
    ),
    sumDeltaTSp_()
{
    readCoeffs(modelDict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList
Foam::Lagrangian::constantFluxCarrierTransfer::addSupFields() const
{
    return
        isCloud<clouds::massive>()
      ? wordList({clouds::massive::mName})
      : wordList({clouds::shaped::vName});
}


bool Foam::Lagrangian::constantFluxCarrierTransfer::addsSupToField
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
     || eqnPhaseName == word::null
     || eqnPhaseName == cloud<clouds::carried>().carrierPhaseName();
}


void Foam::Lagrangian::constantFluxCarrierTransfer::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    if (!final) return;

    const LagrangianSubMesh& subMesh = deltaT.mesh();
    SubList<LagrangianState> subStates = subMesh.sub(mesh().states());

    // Surface area divided by volume to the power of two-thirds. Could/should
    // be made a virtual function of the shaped cloud if non-spherical
    // representations are added.
    static const scalar aByVPowTwoThirds = cbrt(36*constant::mathematical::pi);

    // Calculate the particle's remaining lifetime
    tmp<LagrangianSubScalarField> tlifetime =
        3
       *cbrt(cloud<clouds::shaped>().v(subMesh))
       /aByVPowTwoThirds
       /flux_;
    if (flux_.dimensions() == dimMass/dimArea/dimTime)
    {
        assertCloud
        <
            clouds::coupledToConstantDensityFluid,
            clouds::massive
        >();

        if (isCloud<clouds::coupledToConstantDensityFluid>())
        {
            tlifetime =
                tlifetime/cloud<clouds::coupledToConstantDensityFluid>().rho();
        }
        else
        {
            tlifetime =
                tlifetime/cloud<clouds::massive>().rho(subMesh);
        }
    }

    // Remove any particles with a lifetime shorter than the time-step
    forAll(subMesh, subi)
    {
        if (tlifetime()[subi] < deltaT[subi])
        {
            subStates[subi] = LagrangianState::toBeRemoved;
        }
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::Lagrangian::constantFluxCarrierTransfer::source
(
    const word& fieldName,
    const LagrangianSubMesh& subMesh
)
const
{
    return
      - this->Sp(subMesh)
       *(
            isCloud<clouds::massive>()
          ? cloud<clouds::massive>().m(subMesh)
          : cloud<clouds::shaped>().v(subMesh)
        );
}


void Foam::Lagrangian::constantFluxCarrierTransfer::preAddSup
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


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_LAGRANGIAN_MODEL_ADD_FIELD_SUP,
    Lagrangian::constantFluxCarrierTransfer
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_LAGRANGIAN_MODEL_ADD_V_OR_M_FIELD_SUP,
    Lagrangian::constantFluxCarrierTransfer
)


void Foam::Lagrangian::constantFluxCarrierTransfer::postAddSup
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    if (final) sumDeltaTSp_.clear();
}


bool Foam::Lagrangian::constantFluxCarrierTransfer::read
(
    const dictionary& modelDict
)
{
    if (LagrangianModel::read(modelDict))
    {
        readCoeffs(modelDict);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
