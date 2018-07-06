/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "VoFSolidificationMeltingSource.H"
#include "twoPhaseMixtureThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(VoFSolidificationMeltingSource, 0);

        addToRunTimeSelectionTable
        (
            option,
            VoFSolidificationMeltingSource,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::VoFSolidificationMeltingSource::update()
{
    if (curTimeIndex_ == mesh_.time().timeIndex())
    {
        return;
    }

    if (debug)
    {
        Info<< type() << ": " << name_
            << " - updating solid phase fraction" << endl;
    }

    alphaSolid_.oldTime();

    const twoPhaseMixtureThermo& thermo
    (
        mesh_.lookupObject<twoPhaseMixtureThermo>
        (
            twoPhaseMixtureThermo::dictName
        )
    );

    const volScalarField& TVoF = thermo.thermo1().T();
    const volScalarField CpVoF(thermo.thermo1().Cp());
    const volScalarField& alphaVoF = thermo.alpha1();

    forAll(cells_, i)
    {
        const label celli = cells_[i];

        alphaSolid_[celli] = min
        (
            relax_*alphaVoF[celli]*alphaSolidT_->value(TVoF[celli])
          + (1 - relax_)*alphaSolid_[celli],
            alphaVoF[celli]
        );
    }

    alphaSolid_.correctBoundaryConditions();

    curTimeIndex_ = mesh_.time().timeIndex();
}


Foam::word Foam::fv::VoFSolidificationMeltingSource::alphaSolidName() const
{
    const twoPhaseMixtureThermo& thermo
    (
        mesh_.lookupObject<twoPhaseMixtureThermo>
        (
            twoPhaseMixtureThermo::dictName
        )
    );

    const volScalarField& alphaVoF = thermo.alpha1();

    return IOobject::groupName(alphaVoF.name(), "solid");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFSolidificationMeltingSource::VoFSolidificationMeltingSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    alphaSolidT_(Function1<scalar>::New("alphaSolidT", coeffs_)),
    L_("L", dimEnergy/dimMass, coeffs_),
    relax_(coeffs_.lookupOrDefault("relax", 0.9)),
    Cu_(coeffs_.lookupOrDefault<scalar>("Cu", 100000)),
    q_(coeffs_.lookupOrDefault("q", 0.001)),
    alphaSolid_
    (
        IOobject
        (
            alphaSolidName(),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha1", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    curTimeIndex_(-1)
{
    fieldNames_.setSize(2);
    fieldNames_[0] = "U";
    fieldNames_[1] = "T";
    applied_.setSize(fieldNames_.size(), false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::VoFSolidificationMeltingSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    apply(geometricOneField(), eqn);
}


void Foam::fv::VoFSolidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    apply(rho, eqn);
}


void Foam::fv::VoFSolidificationMeltingSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    update();

    scalarField& Sp = eqn.diag();
    const scalarField& V = mesh_.V();

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        const scalar Vc = V[celli];
        const scalar alphaFluid = 1 - alphaSolid_[celli];

        const scalar S = Cu_*sqr(1 - alphaFluid)/(pow3(alphaFluid) + q_);

        Sp[celli] -= Vc*S;
    }
}


void Foam::fv::VoFSolidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    // Momentum source uses a Boussinesq approximation - redirect
    addSup(eqn, fieldi);
}


// ************************************************************************* //
