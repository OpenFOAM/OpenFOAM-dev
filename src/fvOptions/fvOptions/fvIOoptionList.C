/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "fvIOoptionList.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::fv::IOoptionList::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "fvOptions",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        Info<< "Creating finite volume options from "
            << io.instance()/io.name() << nl
            << endl;

        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        // Check if the fvOptions file is in system
        io.instance() = mesh.time().system();

        if (io.headerOk())
        {
            Info<< "Creating finite volume options from "
                << io.instance()/io.name() << nl
                << endl;

            io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            return io;
        }
        else
        {
            Info<< "No finite volume options present" << nl << endl;

            io.readOpt() = IOobject::NO_READ;
            return io;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::IOoptionList::IOoptionList
(
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(mesh)),
    optionList(mesh, *this)
{}


bool Foam::fv::IOoptionList::read()
{
    if (regIOobject::read())
    {
        optionList::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
