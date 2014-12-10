/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "vorticity.H"
#include "volFields.H"
#include "dictionary.H"
#include "fvcCurl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vorticity, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vorticity::vorticity
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    UName_("U"),
    outputName_(typeName)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "vorticity::vorticity"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    read(dict);

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volVectorField* vorticityPtr
        (
            new volVectorField
            (
                IOobject
                (
                    outputName_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("0", dimless/dimTime, vector::zero)
            )
        );

        mesh.objectRegistry::store(vorticityPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vorticity::~vorticity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vorticity::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
        if (UName_ != "U")
        {
            outputName_ = typeName + "(" + UName_ + ")";
        }
    }
}


void Foam::vorticity::execute()
{
    if (active_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        volVectorField& vorticity =
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>(outputName_)
            );

        vorticity = fvc::curl(U);
    }
}


void Foam::vorticity::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::vorticity::timeSet()
{
    // Do nothing
}


void Foam::vorticity::write()
{
    if (active_)
    {
        const volVectorField& vorticity =
            obr_.lookupObject<volVectorField>(outputName_);

        Info<< type() << " " << name_ << " output:" << nl
            << "    writing field " << vorticity.name() << nl
            << endl;

        vorticity.write();
    }
}


// ************************************************************************* //
