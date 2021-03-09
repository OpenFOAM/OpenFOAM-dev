/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
#include "fvcDdt.H"
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
            fvModel,
            VoFSolidificationMeltingSource,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::VoFSolidificationMeltingSource::readCoeffs()
{
    alphaSolidT_.reset(Function1<scalar>::New("alphaSolidT", coeffs()).ptr());
    L_ = dimensionedScalar("L", dimEnergy/dimMass, coeffs());
    relax_ = coeffs().lookupOrDefault<scalar>("relax", 0.9);
    Cu_ = coeffs().lookupOrDefault<scalar>("Cu", 100000);
    q_ = coeffs().lookupOrDefault<scalar>("q", 0.001);
}


void Foam::fv::VoFSolidificationMeltingSource::update() const
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    if (debug)
    {
        Info<< type() << ": " << name()
            << " - updating solid phase fraction" << endl;
    }

    alphaSolid_.oldTime();

    const twoPhaseMixtureThermo& thermo
    (
        mesh().lookupObject<twoPhaseMixtureThermo>
        (
            twoPhaseMixtureThermo::dictName
        )
    );

    const volScalarField& TVoF = thermo.thermo1().T();
    const volScalarField CpVoF(thermo.thermo1().Cp());
    const volScalarField& alphaVoF = thermo.alpha1();

    const labelList& cells = set_.cells();

    forAll(cells, i)
    {
        const label celli = cells[i];

        alphaSolid_[celli] = min
        (
            relax_*alphaVoF[celli]*alphaSolidT_->value(TVoF[celli])
          + (1 - relax_)*alphaSolid_[celli],
            alphaVoF[celli]
        );
    }

    alphaSolid_.correctBoundaryConditions();

    curTimeIndex_ = mesh().time().timeIndex();
}


Foam::word Foam::fv::VoFSolidificationMeltingSource::alphaSolidName() const
{
    const twoPhaseMixtureThermo& thermo
    (
        mesh().lookupObject<twoPhaseMixtureThermo>
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
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    alphaSolidT_(),
    L_("L", dimEnergy/dimMass, NaN),
    relax_(NaN),
    Cu_(NaN),
    q_(NaN),
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
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    curTimeIndex_(-1)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::VoFSolidificationMeltingSource::addSupFields() const
{
    return wordList({"U", "T"});
}


void Foam::fv::VoFSolidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    update();

    const twoPhaseMixtureThermo& thermo
    (
        mesh().lookupObject<twoPhaseMixtureThermo>
        (
            twoPhaseMixtureThermo::dictName
        )
    );

    const volScalarField CpVoF(thermo.thermo1().Cp());

    if (eqn.psi().dimensions() == dimTemperature)
    {
        eqn += L_/CpVoF*(fvc::ddt(rho, alphaSolid_));
    }
    else
    {
        eqn += L_*(fvc::ddt(rho, alphaSolid_));
    }
}


void Foam::fv::VoFSolidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    update();

    scalarField& Sp = eqn.diag();
    const scalarField& V = mesh().V();

    const labelList& cells = set_.cells();

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar Vc = V[celli];
        const scalar alphaFluid = 1 - alphaSolid_[celli];

        const scalar S = Cu_*sqr(1 - alphaFluid)/(pow3(alphaFluid) + q_);

        Sp[celli] -= Vc*rho[celli]*S;
    }
}


void Foam::fv::VoFSolidificationMeltingSource::updateMesh
(
    const mapPolyMesh& mpm
)
{
    set_.updateMesh(mpm);
}


bool Foam::fv::VoFSolidificationMeltingSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
