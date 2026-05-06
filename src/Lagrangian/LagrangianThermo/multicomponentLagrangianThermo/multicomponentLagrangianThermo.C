/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "LagrangianFieldsFwd.H"
#include "multicomponentLagrangianThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multicomponentLagrangianThermo, 0);
    defineRunTimeSelectionTable(multicomponentLagrangianThermo, LagrangianMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentLagrangianThermo::implementation::implementation
(
    const dictionary& dict,
    const speciesTable& species,
    const LagrangianMesh& mesh,
    const word& phaseName
)
:
    defaultSpecieName_
    (
        species.size()
      ? dict.lookup<word>("defaultSpecie")
      : "undefined"
    ),
    defaultSpeciei_
    (
        species.size()
      ? species[defaultSpecieName_]
      : -1
    ),
    Y_(species.size())
{
    tmp<LagrangianScalarField> Ydefault;

    bool Yset = false;

    // Read the species' mass fractions
    forAll(species, i)
    {
        typeIOobject<LagrangianScalarDynamicField> io
        (
            IOobject::groupName(species[i], phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ
        );

        if (io.headerOk())
        {
            Y_.set
            (
                i,
                new LagrangianScalarDynamicField
                (
                    IOobject
                    (
                        IOobject::groupName(species[i], phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );

            Yset = true;
        }
        else
        {
            if (!Ydefault.valid())
            {
                Ydefault = new LagrangianScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("Ydefault", phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );
            }

            Y_.set
            (
                i,
                new LagrangianScalarDynamicField
                (
                    IOobject
                    (
                        IOobject::groupName(species[i], phaseName),
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault()
                )
            );
        }
    }

    // If the mass fractions have been specified check and normalise
    if (Yset)
    {
        // Scale the mass fractions so that they sum to one in case the input is
        // not exact
        tmp<LagrangianScalarField> tYt
        (
            new LagrangianScalarField
            (
                IOobject::groupName("Yt", phaseName),
                Y_[0]
            )
        );
        LagrangianScalarField& Yt = tYt.ref();

        for (label i = 1; i < Y_.size(); ++ i)
        {
            Yt += Y_[i];
        }

        if (min(Yt.primitiveField()) == 0 && max(Yt.primitiveField()) == 0)
        {
            FatalErrorInFunction
                << "Sum of specie mass fractions = 0"
                << exit(FatalError);
        }

        if (min(Yt.primitiveField()) < 0.999)
        {
            FatalErrorInFunction
                << "Min sum of specie mass fractions "
                << min(Yt.primitiveField())
                << " < 0.999"
                << exit(FatalError);
        }

        if (max(Yt.primitiveField()) > 1.001)
        {
            FatalErrorInFunction
                << "Max sum of specie mass fractions "
                << max(Yt.primitiveField())
                << " > 1.001"
                << exit(FatalError);
        }

        forAll(Y_, i)
        {
            Y_[i] /= Yt;
        }
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multicomponentLagrangianThermo>
Foam::multicomponentLagrangianThermo::New
(
    const LagrangianMesh& mesh,
    const word& phaseName
)
{
    return
        basicLagrangianThermo::New<multicomponentLagrangianThermo>
        (
            mesh,
            phaseName
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multicomponentLagrangianThermo::~multicomponentLagrangianThermo()
{}


Foam::multicomponentLagrangianThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label
Foam::multicomponentLagrangianThermo::implementation::defaultSpecie() const
{
    return defaultSpeciei_;
}


Foam::PtrList<Foam::LagrangianScalarDynamicField>&
Foam::multicomponentLagrangianThermo::implementation::Y()
{
    return Y_;
}


const Foam::PtrList<Foam::LagrangianScalarDynamicField>&
Foam::multicomponentLagrangianThermo::implementation::Y() const
{
    return Y_;
}


void Foam::multicomponentLagrangianThermo::implementation::normaliseY
(
    const LagrangianSubMesh& subMesh
)
{
    if (defaultSpeciei_ == -1 || !species().size()) return;

    tmp<LagrangianSubScalarField> tYt
    (
        LagrangianSubScalarField::New
        (
            IOobject::groupName("Yt", phaseName()),
            subMesh,
            dimensionedScalar(dimless, 0)
        )
    );
    LagrangianSubScalarField& Yt = tYt.ref();

    forAll(Y(), i)
    {
        if (i == defaultSpeciei_) continue;

        LagrangianSubScalarSubField subYi(subMesh.sub(Y_[i]));
        subYi = max(subYi, scalar(0));
        Yt += subYi;
    }

    LagrangianSubScalarSubField subYdefault(subMesh.sub(Y_[defaultSpeciei_]));
    subYdefault = max(scalar(1) - Yt, scalar(0));
}


// ************************************************************************* //
