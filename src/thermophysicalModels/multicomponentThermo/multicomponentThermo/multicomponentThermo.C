/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::multicomponentThermo::implementation::correctMassFractions()
{
    if (species_.size())
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
            Yt += Y_[i];
        }

        if (mag(min(Yt).value()) < rootVSmall)
        {
            FatalErrorInFunction
                << "Sum of mass fractions is zero for species " << species()
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
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    species_(specieNames),
    defaultSpecieName_
    (
        species_.size()
      ? dict.lookupBackwardsCompatible<word>
        (
            {"defaultSpecie", "inertSpecie"}
        )
      : word("undefined")
    ),
    defaultSpeciei_
    (
        species_.size()
     && species_.found(defaultSpecieName_)
      ? species_[defaultSpecieName_]
      : -1
    ),
    active_(species_.size(), true),
    Y_(species_.size())
{
    if (species_.size() && defaultSpeciei_ == -1)
    {
        FatalIOErrorInFunction(dict)
            << "default specie " << defaultSpecieName_
            << " not found in available species " << species_
            << exit(FatalIOError);
    }

    // Read the species' mass fractions
    tmp<volScalarField> tYdefault;
    forAll(species_, i)
    {
        typeIOobject<volScalarField> header
        (
            IOobject::groupName(species_[i], phaseName),
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
                        IOobject::groupName(species_[i], phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            // Read Ydefault if not already read
            if (!tYdefault.valid())
            {
                const word YdefaultName
                (
                    IOobject::groupName("Ydefault", phaseName)
                );

                typeIOobject<volScalarField> timeIO
                (
                    YdefaultName,
                    mesh.time().name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                typeIOobject<volScalarField> constantIO
                (
                    YdefaultName,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                typeIOobject<volScalarField> time0IO
                (
                    YdefaultName,
                    Time::timeName(0),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (timeIO.headerOk())
                {
                    tYdefault = new volScalarField(timeIO, mesh);
                }
                else if (constantIO.headerOk())
                {
                    tYdefault = new volScalarField(constantIO, mesh);
                }
                else
                {
                    tYdefault = new volScalarField(time0IO, mesh);
                }
            }

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tYdefault()
                )
            );
        }
    }

    correctMassFractions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multicomponentThermo::~multicomponentThermo()
{}


Foam::multicomponentThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::speciesTable&
Foam::multicomponentThermo::implementation::species() const
{
    return species_;
}


Foam::label Foam::multicomponentThermo::implementation::defaultSpecie() const
{
    return defaultSpeciei_;
}


const Foam::List<bool>&
Foam::multicomponentThermo::implementation::speciesActive() const
{
    return active_;
}


void Foam::multicomponentThermo::syncSpeciesActive() const
{
    if (Pstream::parRun())
    {
        List<bool>& speciesActive =
            const_cast<List<bool>&>(this->speciesActive());

        Pstream::listCombineGather(speciesActive, orEqOp<bool>());
        Pstream::listCombineScatter(speciesActive);

        PtrList<volScalarField>& Y =
            const_cast<PtrList<volScalarField>&>(this->Y());

        forAll(Y, speciei)
        {
            if (speciesActive[speciei])
            {
                Y[speciei].writeOpt() = IOobject::AUTO_WRITE;
            }
        }
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


void Foam::multicomponentThermo::normaliseY()
{
    if (species().size())
    {
        tmp<volScalarField> tYt
        (
            volScalarField::New
            (
                IOobject::groupName("Yt", phaseName()),
                Y()[0].mesh(),
                dimensionedScalar(dimless, 0)
            )
        );
        volScalarField& Yt = tYt.ref();

        forAll(Y(), i)
        {
            if (solveSpecie(i))
            {
                Y()[i].max(scalar(0));
                Yt += Y()[i];
            }
        }

        Y()[defaultSpecie()] = scalar(1) - Yt;
        Y()[defaultSpecie()].max(scalar(0));
    }
}


// ************************************************************************* //
