/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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

#include "solidificationMeltingSource.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "uniformDimensionedFields.H"
#include "zeroGradientFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        fv::solidificationMeltingSource::thermoMode,
        2
    >::names[] =
    {
        "thermo",
        "lookup"
    };

    namespace fv
    {
        defineTypeNameAndDebug(solidificationMeltingSource, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            solidificationMeltingSource,
            dictionary
        );
    }
}

const Foam::NamedEnum<Foam::fv::solidificationMeltingSource::thermoMode, 2>
    Foam::fv::solidificationMeltingSource::thermoModeTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::solidificationMeltingSource::readCoeffs()
{
    Tsol_ = coeffs().lookup<scalar>("Tsol");
    Tliq_ = coeffs().lookupOrDefault<scalar>("Tliq", Tsol_);
    alpha1e_ = coeffs().lookupOrDefault<scalar>("alpha1e", 0.0);
    L_ = coeffs().lookup<scalar>("L");

    relax_ = coeffs().lookupOrDefault("relax", 0.9);

    mode_ = thermoModeTypeNames_.read(coeffs().lookup("thermoMode"));

    rhoRef_ = coeffs().lookup<scalar>("rhoRef");
    TName_ = coeffs().lookupOrDefault<word>("T", "T");
    CpName_ = coeffs().lookupOrDefault<word>("Cp", "Cp");
    UName_ = coeffs().lookupOrDefault<word>("U", "U");
    phiName_ = coeffs().lookupOrDefault<word>("phi", "phi");

    Cu_ = coeffs().lookupOrDefault<scalar>("Cu", 100000);
    q_ = coeffs().lookupOrDefault("q", 0.001);

    beta_ = coeffs().lookup<scalar>("beta");
}


Foam::tmp<Foam::volScalarField>
Foam::fv::solidificationMeltingSource::Cp() const
{
    switch (mode_)
    {
        case thermoMode::thermo:
        {
            const basicThermo& thermo =
                mesh().lookupObject<basicThermo>(physicalProperties::typeName);

            return thermo.Cp();
            break;
        }
        case thermoMode::lookup:
        {
            if (CpName_ == "CpRef")
            {
                scalar CpRef = coeffs().lookup<scalar>("CpRef");

                return volScalarField::New
                (
                    name() + ":Cp",
                    mesh(),
                    dimensionedScalar
                    (
                        dimEnergy/dimMass/dimTemperature,
                        CpRef
                    ),
                    extrapolatedCalculatedFvPatchScalarField::typeName
                );
            }
            else
            {
                return mesh().lookupObject<volScalarField>(CpName_);
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled thermo mode: " << thermoModeTypeNames_[mode_]
                << abort(FatalError);
        }
    }

    return tmp<volScalarField>(nullptr);
}


Foam::vector Foam::fv::solidificationMeltingSource::g() const
{
    if (mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        const uniformDimensionedVectorField& value =
            mesh().lookupObject<uniformDimensionedVectorField>("g");
        return value.value();
    }
    else
    {
        return coeffs().lookup("g");
    }
}


void Foam::fv::solidificationMeltingSource::update
(
    const volScalarField& Cp
) const
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    if (debug)
    {
        Info<< type() << ": " << name()
            << " - updating phase indicator" << endl;
    }

    // update old time alpha1 field
    alpha1_.oldTime();

    const volScalarField& T = mesh().lookupObject<volScalarField>(TName_);

    const labelList& cells = set_.cells();

    forAll(cells, i)
    {
        const label celli = cells[i];

        const scalar Tc = T[celli];
        const scalar Cpc = Cp[celli];
        const scalar alpha1New =
            alpha1_[celli]
          + relax_*Cpc
           *(
                Tc
              - max
                (
                    Tsol_,
                    Tsol_
                  + (Tliq_ - Tsol_)*(alpha1_[celli] - alpha1e_)/(1 - alpha1e_)
                )
            )/L_;

        alpha1_[celli] = max(0, min(alpha1New, 1));
        deltaT_[i] =
            Tc
          - max
            (
                Tsol_,
                Tsol_
              + (Tliq_ - Tsol_)*(alpha1_[celli] - alpha1e_)/(1 - alpha1e_)
            );
    }

    alpha1_.correctBoundaryConditions();

    curTimeIndex_ = mesh().time().timeIndex();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::solidificationMeltingSource::solidificationMeltingSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    set_(mesh, coeffs()),
    Tsol_(NaN),
    Tliq_(NaN),
    alpha1e_(NaN),
    L_(NaN),
    relax_(NaN),
    mode_(thermoMode::thermo),
    rhoRef_(NaN),
    TName_(word::null),
    CpName_(word::null),
    UName_(word::null),
    phiName_(word::null),
    Cu_(NaN),
    q_(NaN),
    beta_(NaN),
    alpha1_
    (
        IOobject
        (
            this->name() + ":alpha1",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    curTimeIndex_(-1),
    deltaT_(set_.cells().size(), 0)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::solidificationMeltingSource::addSupFields() const
{
    switch (mode_)
    {
        case thermoMode::thermo:
        {
            const basicThermo& thermo =
                mesh().lookupObject<basicThermo>(physicalProperties::typeName);

            return wordList({UName_, thermo.he().name()});
        }
        case thermoMode::lookup:
        {
            return wordList({UName_, TName_});
        }
    }

    return wordList::null();
}


void Foam::fv::solidificationMeltingSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    apply(geometricOneField(), eqn);
}


void Foam::fv::solidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    apply(rho, eqn);
}


void Foam::fv::solidificationMeltingSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField Cp(this->Cp());

    update(Cp);

    vector g = this->g();

    scalarField& Sp = eqn.diag();
    vectorField& Su = eqn.source();
    const scalarField& V = mesh().V();

    const labelList& cells = set_.cells();

    forAll(cells, i)
    {
        const label celli = cells[i];

        const scalar Vc = V[celli];
        const scalar alpha1c = alpha1_[celli];

        const scalar S = -Cu_*sqr(1.0 - alpha1c)/(pow3(alpha1c) + q_);
        const vector Sb = rhoRef_*g*beta_*deltaT_[i];

        Sp[celli] += Vc*S;
        Su[celli] += Vc*Sb;
    }
}


void Foam::fv::solidificationMeltingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    // Momentum source uses a Boussinesq approximation - redirect
    addSup(eqn, fieldName);
}


bool Foam::fv::solidificationMeltingSource::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::solidificationMeltingSource::topoChange
(
    const polyTopoChangeMap& map
)
{
    set_.topoChange(map);
}


void Foam::fv::solidificationMeltingSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::solidificationMeltingSource::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


bool Foam::fv::solidificationMeltingSource::read(const dictionary& dict)
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
