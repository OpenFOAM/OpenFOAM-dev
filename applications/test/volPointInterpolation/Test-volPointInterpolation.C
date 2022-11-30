/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Application
    volPointInterpolationTest

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool interpolate(const fvMesh& mesh, const word& name)
{
    typeIOobject<VolField<Type>> io
    (
        name,
        mesh.time().name(),
        mesh,
        IOobject::MUST_READ
    );

    if (!io.headerOk()) return false;

    Info<< "Reading field " << name << nl << endl;

    const VolField<Type> vf(io, mesh);
    const PointField<Type> pf(volPointInterpolation::New(mesh).interpolate(vf));

    Info<< "Writing field " << pf.name() << nl << endl;

    return pf.write();
}


int main(int argc, char *argv[])
{
    argList::validArgs.append("field");

    Foam::timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const pointMesh& pMesh = pointMesh::New(mesh);
    const pointBoundaryMesh& pbm = pMesh.boundary();

    Info<< "pointMesh boundary" << nl;
    forAll(pbm, patchi)
    {
        Info<< "patch=" << pbm[patchi].name()
            << ", type=" << pbm[patchi].type()
            << ", coupled=" << pbm[patchi].coupled()
            << endl;
    }

    const word name = args.argRead<word>(1);

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.userTimeName() << endl;

        mesh.readUpdate();

        if
        (
            !interpolate<scalar>(mesh, name)
         && !interpolate<vector>(mesh, name)
         && !interpolate<sphericalTensor>(mesh, name)
         && !interpolate<symmTensor>(mesh, name)
         && !interpolate<tensor>(mesh, name)
        )
        {
            WarningInFunction
                << "Could not find field " << name << nl << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
