/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "phaseChange.H"
#include "fluidMulticomponentThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(phaseChange, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::phaseChange::readCoeffs(const dictionary& dict)
{
    energySemiImplicit_ = dict.lookup<bool>("energySemiImplicit");
}


const Foam::List<Foam::labelPair> Foam::fv::phaseChange::initSpecieis() const
{
    const ThermoRefPair<multicomponentThermo>& mcThermos =
        thermos().thermos<multicomponentThermo>();

    List<labelPair> result(species().size(), labelPair(-1, -1));

    forAll(phaseNames(), i)
    {
        if (mcThermos.valid()[i])
        {
            forAll(species(), mDoti)
            {
                result[mDoti][i] = mcThermos[i].species()[species()[mDoti]];
            }
        }
    }

    return result;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::wordList Foam::fv::phaseChange::readSpecie
(
    const dictionary& dict,
    const bool required
) const
{
    const bool haveSpecie = dict.found("specie");

    if (required && !haveSpecie)
    {
        dict.lookup<word>("specie");
    }

    const wordList result =
        haveSpecie
      ? wordList(1, dict.lookup<word>("specie"))
      : wordList();

    return result;
}


Foam::wordList Foam::fv::phaseChange::readSpecies
(
    const dictionary& dict,
    const bool required
) const
{
    const bool haveSpecie = dict.found("specie");
    const bool haveSpecies = dict.found("species");

    if (haveSpecie && haveSpecies)
    {
        FatalIOErrorInFunction(dict)
            << "Both keywords specie and species "
            << " are defined in dictionary " << dict.name()
            << exit(FatalError);
    }

    if (required && !haveSpecie && !haveSpecies)
    {
        FatalIOErrorInFunction(dict)
            << "Neither keywords specie or species "
            << " are defined in dictionary " << dict.name()
            << exit(FatalError);
    }

    const wordList result =
        haveSpecie ? wordList(1, dict.lookup<word>("specie"))
      : haveSpecies ? dict.lookup<wordList>("species")
      : wordList();

    return result;
}


void Foam::fv::phaseChange::reReadSpecie(const dictionary& dict) const
{
    if (species() != readSpecie(dict, false))
    {
        FatalIOErrorInFunction(dict)
            << "Cannot change the specie of a " << type() << " model "
            << "at run time" << exit(FatalIOError);
    }
}


void Foam::fv::phaseChange::reReadSpecies(const dictionary& dict) const
{
    if (species() != readSpecies(dict, false))
    {
        FatalIOErrorInFunction(dict)
            << "Cannot change the species of a " << type() << " model "
            << "at run time" << exit(FatalIOError);
    }
}


void Foam::fv::phaseChange::setSpecies
(
    const word& name,
    const word& modelType,
    const wordList& species
)
{
    const ThermoRefPair<multicomponentThermo>& mcThermos =
        thermos().thermos<multicomponentThermo>();

    // Set the names
    species_ = species;

    // Set the indices
    specieis_ = List<labelPair>(species.size(), labelPair(-1, -1));
    forAll(phaseNames(), i)
    {
        if (mcThermos.valid()[i])
        {
            forAll(species, mDoti)
            {
                specieis_[mDoti][i] = mcThermos[i].species()[species[mDoti]];
            }
        }
    }

    // Checks ...

    // If either phase is multicomponent then species should
    // have been specified
    if (mcThermos.either() && species_.empty())
    {
        FatalErrorInFunction
            << "Mixture transfer specified by model " << name << " of type "
            << modelType << " for two phases " << phaseNames().first()
            << " and " << phaseNames().second() << " but ";

        mcThermos.both()
      ? FatalErrorInFunction
            << "both phases have"
      : FatalErrorInFunction
            << "phase " << phaseNames()[mcThermos.valid().second()] << " has";

        FatalErrorInFunction
            << " multiple species" << exit(FatalError);
    }

    // If neither phase is multicomponent then species should
    // not have been specified
    if (!mcThermos.either() && species_.size())
    {
        FatalErrorInFunction
            << "Specie transfer specified by model " << name
            << " of type " << modelType << " for two pure phases "
            << phaseNames().first() << " and " << phaseNames().second()
            << exit(FatalError);
    }

    // If either phase is pure then there can be at most one specie
    if (!mcThermos.both() && species_.size() > 1)
    {
        FatalErrorInFunction
            << "Multi-specie transfer specified by model " << name
            << " of type " << modelType << " for phases "
            << phaseNames().first() << " and " << phaseNames().second()
            << " but phase " << phaseNames()[mcThermos.valid().first()]
            << " is pure " << exit(FatalError);
    }
}


void Foam::fv::phaseChange::setSpecies(const wordList& species)
{
    setSpecies(name(), type(), species);
}


void Foam::fv::phaseChange::reSetSpecies(const wordList& species)
{
    if (this->species() != species)
    {
        FatalErrorInFunction
            << "Cannot change the species of a " << type() << " model "
            << "at run time" << exit(FatalError);
    }
}


const Foam::volScalarField& Foam::fv::phaseChange::p() const
{
    for (label i = 0; i < 2; ++ i)
    {
        if (isA<fluidThermo>(thermos()[i]))
        {
            return refCast<const fluidThermo>(thermos()[i]).p();
        }
    }

    return mesh().lookupObject<volScalarField>("p");
}


Foam::tmp<Foam::volScalarField> Foam::fv::phaseChange::vifToVf
(
    const tmp<volScalarField::Internal>& tvif
)
{
    tmp<volScalarField> tvf =
        volScalarField::New
        (
            tvif().name(),
            tvif().mesh(),
            tvif().dimensions(),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        );

    tvf->internalFieldRef() = tvif();
    tvf->correctBoundaryConditions();

    tvif.clear();

    return tvf;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::phaseChange::vfToVif
(
    const tmp<volScalarField>& tvf
)
{
    tmp<volScalarField::Internal> tvif(tvf.ptr());

    tvf.clear();

    return tvif;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseChange::phaseChange
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const wordList& species
)
:
    massTransfer(name, modelType, mesh, dict),
    thermos_(mesh, phaseNames()),
    heNames_(thermos_.first().he().name(), thermos_.second().he().name()),
    species_(),
    specieis_(),
    energySemiImplicit_(false)
{
    readCoeffs(coeffs(dict));

    if (notNull(species)) setSpecies(name, modelType, species);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::ThermoRefPair<Foam::fluidThermo>
Foam::fv::phaseChange::fluidThermos
(
    const bool firstRequired,
    const bool secondRequired
) const
{
    return
        thermos_.thermos<fluidThermo>
        (
            {firstRequired, secondRequired},
            *this,
            "fluid"
        );
}


const Foam::ThermoRefPair<Foam::multicomponentThermo>
Foam::fv::phaseChange::multicomponentThermos
(
    const bool firstRequired,
    const bool secondRequired
) const
{
    return
        thermos_.thermos<multicomponentThermo>
        (
            {firstRequired, secondRequired},
            *this,
            "multicomponent"
        );
}


const Foam::ThermoRefPair<Foam::fluidMulticomponentThermo>
Foam::fv::phaseChange::fluidMulticomponentThermos
(
    const bool firstRequired,
    const bool secondRequired
) const
{
    return
        thermos_.thermos<fluidMulticomponentThermo>
        (
            {firstRequired, secondRequired},
            *this,
            "fluid-multicomponent"
        );
}


const Foam::labelPair& Foam::fv::phaseChange::specieis(const label mDoti) const
{
    static const labelPair noSpecieis(-1, -1);

    if (mDoti == -1)
    {
        if (species().empty())
        {
            return noSpecieis;
        }
        else if (species().size() == 1)
        {
            return specieis_[0];
        }
        else
        {
            FatalErrorInFunction
                << "Requested mixture/single-specie indices from multi-specie "
                << "model of type" << type() << exit(FatalError);
        }
    }

    return specieis_[mDoti];
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseChange::Tchange() const
{
    tmp<volScalarField::Internal> tTchange =
        volScalarField::Internal::New
        (
            name() + ":Tchange",
            mesh(),
            dimTemperature
        );
    volScalarField::Internal& Tchange = tTchange.ref();

    tmp<volScalarField::Internal> tmDot = this->mDot();
    const volScalarField::Internal& mDot = tmDot();

    const volScalarField::Internal& T1 = thermos().first().T();
    const volScalarField::Internal& T2 = thermos().second().T();

    forAll(Tchange, i)
    {
        Tchange[i] =
            mDot[i] > 0 ? T1[i]
          : mDot[i] < 0 ? T2[i]
          : (T1[i] + T2[i])/2;
    }

    return tTchange;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseChange::Lfraction() const
{
    const volScalarField::Internal& kappa1 = thermos().first().kappa();
    const volScalarField::Internal& kappa2 = thermos().second().kappa();

    return kappa2/(kappa1 + kappa2);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::phaseChange::L
(
    const label mDoti
) const
{
    return L(this->Tchange(), mDoti);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::phaseChange::L
(
    const volScalarField::Internal& Tchange,
    const label mDoti
) const
{
    const ThermoRefPair<multicomponentThermo>& mcThermos =
        thermos().thermos<multicomponentThermo>();

    const labelPair specieis = this->specieis(mDoti);

    const volScalarField::Internal& p = this->p();

    // Absolute enthalpies at the interface
    Pair<tmp<volScalarField::Internal>> has;
    for (label j = 0; j < 2; ++ j)
    {
        has[j] =
            specieis[j] == -1
          ? thermos()[j].ha(p, Tchange)
          : mcThermos[j].hai(specieis[j], p, Tchange);
    }

    // Latent heat of phase change
    return has.second() - has.first();
}


Foam::tmp<Foam::scalarField> Foam::fv::phaseChange::L
(
    const label patchi,
    const scalarField& Tchange,
    const label mDoti
) const
{
    const ThermoRefPair<multicomponentThermo>& mcThermos =
        thermos().thermos<multicomponentThermo>();

    const labelPair specieis = this->specieis(mDoti);

    const fvPatchScalarField& p = this->p().boundaryField()[patchi];

    // Absolute enthalpies at the interface
    Pair<tmp<scalarField>> has;
    for (label j = 0; j < 2; ++ j)
    {
        has[j] =
            specieis[j] == -1
          ? thermos()[j].ha(Tchange, patchi)
          : mcThermos[j].hai(specieis[j], p, Tchange);
    }

    // Latent heat of phase change
    return has.second() - has.first();
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::phaseChange::mDot() const
{
    if (species().empty())
    {
        FatalErrorInFunction
            << "Mixture phase change rate not defined by model of type "
            << type() << exit(FatalError);
    }

    tmp<volScalarField::Internal> tmDot =
        volScalarField::Internal::New
        (
            "mDot",
            mesh(),
            dimensionedScalar(dimDensity/dimTime, Zero)
        );

    forAll(species(), mDoti)
    {
        tmDot.ref() += mDot(mDoti);
    }

    return tmDot;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::phaseChange::mDot
(
    const label mDoti
) const
{
    if (mDoti == -1)
    {
        return mDot();
    }

    if (mDoti == 0 && species().size() == 1)
    {
        return mDot();
    }

    if (species().size() > 1)
    {
        FatalErrorInFunction
            << "Specie phase change rate not defined by model of type "
            << type() << exit(FatalError);
    }

    return tmp<volScalarField::Internal>(nullptr);
}


void Foam::fv::phaseChange::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& heOrYi,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(phaseNames(), alpha.group());
    const label s = sign(phaseNames(), alpha.group());

    const ThermoRefPair<multicomponentThermo>& mcThermos =
        thermos().thermos<multicomponentThermo>();

    // Energy equation
    if (index(heNames(), heOrYi.name()) != -1)
    {
        const volScalarField::Internal& p = this->p();
        tmp<volScalarField::Internal> tTchange = this->Tchange();

        for
        (
            label mDoti = species().empty() ? -1 : 0;
            mDoti < species().size();
            mDoti ++
        )
        {
            const labelPair specieis = this->specieis(mDoti);
            tmp<volScalarField::Internal> tmDot = this->mDot(mDoti);

            // Direct transfer of energy due to mass transfer
            eqn +=
                s
               *tmDot()
               *(
                    specieis[i] == -1
                  ? thermos()[i].hs(p, tTchange())
                  : mcThermos[i].hsi(specieis[i], p, tTchange())
               );

            // Optional linearisation
            if (energySemiImplicit_)
            {
                eqn += -fvm::SuSp(-s*tmDot(), heOrYi) - s*tmDot()*heOrYi;
            }

            // Latent heat of phase change
            eqn -=
                (i == 0 ? 1 - Lfraction() : Lfraction())
               *tmDot
               *L(tTchange(), mDoti);
        }

        return;
    }

    // Mass fraction equation
    const word specieName = heOrYi.member();
    if (mcThermos.valid()[i] && mcThermos[i].containsSpecie(specieName))
    {
        // A non-transferring specie. Do not add a source.
        if (!species().found(specieName)) return;

        // A transferring specie. Add a source.
        eqn += s*mDot(species()[specieName]);

        return;
    }

    // Something else. Fall back.
    massTransfer::addSup(alpha, rho, heOrYi, eqn);
}


bool Foam::fv::phaseChange::read(const dictionary& dict)
{
    if (massTransfer::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
