/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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

#include "multicomponentThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multicomponentThermo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multicomponentThermo::implementation::scaleMassFractions()
{
    if (Y_.size())
    {
        tmp<volScalarField> tYt
        (
            volScalarField::New
            (
                IOobject::groupName("Yt", Y_[0].group()),
                Y_[0],
                calculatedFvPatchScalarField::typeName
            )
        );
        volScalarField& Yt = tYt.ref();

        for (label i=1; i<Y_.size(); i++)
        {
            Y_[i].bound(0, 1);
            Yt += Y_[i];
        }

        const scalar minYt = min(Yt).value();
        const scalar maxYt = max(Yt).value();

        if (minYt == 0 && maxYt == 0)
        {
            FatalErrorInFunction
                << "Sum of specie mass fractions = 0"
                << exit(FatalError);
        }

        if (minYt < 1 - normalisationLimit_)
        {
            FatalErrorInFunction
                << "Min sum of specie mass fractions " << minYt
                << " < " << 1 - normalisationLimit_
                << exit(FatalError);
        }

        if (maxYt > 1 + normalisationLimit_)
        {
            FatalErrorInFunction
                << "Max sum of specie mass fractions " << maxYt
                << " > " << 1 + normalisationLimit_
                << exit(FatalError);
        }

        forAll(Y_, i)
        {
            Y_[i] /= Yt;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentThermo::implementation::implementation
(
    const dictionary& dict,
    const speciesTable& species,
    const fvMesh& mesh,
    const word& phaseName,
    const bool normal
)
:
    defaultSpecieName_
    (
        normal && species.size()
      ? dict.lookupOrDefault<word>("defaultSpecie", "undefined")
      : word("undefined")
    ),
    defaultSpeciei_
    (
        normal
      ? species.size() && dict.found("defaultSpecie")
          ? species[defaultSpecieName_]
          : -1
      : -2
    ),
    normalisationLimit_
    (
        dict.lookupOrDefault<scalar>("normalisationLimit", 1e-3)
    ),
    Y_(species.size())
{
    bool Yset = false;

    // Read the species' mass fractions
    forAll(species, i)
    {
        typeIOobject<volScalarField> header
        (
            IOobject::groupName(species[i], phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ
        );

        if (header.headerOk())
        {
            // Read the mass fraction field
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species[i], phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimless
                )
            );

            Yset = true;
        }
        else
        {
            // Read Ydefault if not already read
            if (!Ydefault_.valid())
            {
                Ydefault_ = new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Ydefault", phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimless
                );
            }

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species[i], phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault_()
                )
            );
        }
    }

    // If the mass fractions have been specified
    // check and normalise if normal
    if (Yset && defaultSpeciei_ != -2)
    {
        scaleMassFractions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multicomponentThermo::~multicomponentThermo()
{}


Foam::multicomponentThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::multicomponentThermo::implementation::defaultSpecie() const
{
    return defaultSpeciei_;
}


void Foam::multicomponentThermo::implementation::syncActiveSpecies() const
{
    if (Pstream::parRun())
    {
        boolList& activeSpecies =
            const_cast<boolList&>(this->activeSpecies());

        Pstream::listCombineGather(activeSpecies, orEqOp());
        Pstream::listCombineScatter(activeSpecies);

        PtrList<volScalarField>& Y =
            const_cast<PtrList<volScalarField>&>(this->Y());

        forAll(Y, speciei)
        {
            if (activeSpecies[speciei])
            {
                Y[speciei].writeOpt() = IOobject::AUTO_WRITE;
            }
        }
    }

    if (Ydefault_.valid())
    {
        Ydefault_->writeOpt() = IOobject::AUTO_WRITE;
    }
}


Foam::PtrList<Foam::volScalarField>&
Foam::multicomponentThermo::implementation::Y()
{
    return Y_;
}


const Foam::PtrList<Foam::volScalarField>&
Foam::multicomponentThermo::implementation::Y() const
{
    return Y_;
}


void Foam::multicomponentThermo::implementation::normaliseY()
{
    if (defaultSpeciei_ >= 0)
    {
        if (species().size())
        {
            tmp<volScalarField> tYt
            (
                volScalarField::New
                (
                    IOobject::groupName("Yt", phaseName()),
                    Y_[0].mesh(),
                    dimensionedScalar(dimless, 0)
                )
            );
            volScalarField& Yt = tYt.ref();

            forAll(Y_, i)
            {
                if (solveSpecie(i))
                {
                    Y_[i].bound(0, 1);
                    Yt += Y_[i];
                }
            }

            Y_[defaultSpeciei_] = scalar(1) - Yt;
            Y_[defaultSpeciei_].boundLower(0);
        }

        if (Ydefault_.valid() && Ydefault_->writeOpt() == IOobject::NO_WRITE)
        {
            Ydefault_.clear();
        }
    }
    else if (defaultSpeciei_ == -1)
    {
        scaleMassFractions();
    }
    else
    {
        forAll(Y_, i)
        {
            Y_[i].bound(0, 1);
        }
    }
}


// ************************************************************************* //
