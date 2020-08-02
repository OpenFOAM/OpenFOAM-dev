/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "fvOptions.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(options, 0);
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::fv::options::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
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

        if (io.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Creating finite volume options from "
                << io.instance()/io.name() << nl
                << endl;

            io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
            return io;
        }
        else
        {
            io.readOpt() = IOobject::NO_READ;
            return io;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::options::options
(
    const fvMesh& mesh
)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, options>
    (
        mesh,
        createIOobject(mesh)
    ),
    optionList(mesh)
{
    dictionary::name() = IOobject::objectPath();
    readHeaderOk(IOstream::ASCII, typeName);

    // Add file watch on the fvOptions dictionary for MUST_READ_IF_MODIFIED
    addWatch();

    optionList::reset(*this);
}


bool Foam::fv::options::readData(Istream& is)
{
    is >> *this;
    return !is.bad();
}


bool Foam::fv::options::read()
{
    if (IOdictionary::regIOobject::read())
    {
        optionList::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::fv::options::writeData(Ostream& os) const
{
    dictionary::write(os, false);
    return os.good();
}


void Foam::fv::options::updateMesh(const mapPolyMesh& mpm)
{
    optionList::updateMesh(mpm);
}


bool Foam::fv::options::movePoints()
{
    return optionList::movePoints();
}


// ************************************************************************* //
