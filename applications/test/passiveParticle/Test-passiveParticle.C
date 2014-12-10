/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    testPassiveParticle

Description
    Test cloud of passive particles.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "passiveParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("cloudName");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    runTime.functionObjects().off();

    const word cloudName = args[1];

    {
        // Start with empty cloud
        passiveParticleCloud particles
        (
            mesh,
            cloudName,
            IDLList<passiveParticle>()
        );
        Pout<< "Starting particles:" << particles.size() << endl;

        Pout<< "Adding a particle." << endl;
        particles.addParticle(new passiveParticle(mesh, vector::zero, -1));

        forAllConstIter(passiveParticleCloud, particles, iter)
        {
            Pout<< "    " << iter().position() << " cell:" << iter().cell()
                << " origProc:" << iter().origProc()
                << " origId:" << iter().origId()
                << endl;
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime++;
        Pout<< "Writing particles to time " << runTime.timeName() << endl;
        particles.write();
    }

    {
        Pout<< "Rereading particles from time " << runTime.timeName()
            << endl;
        passiveParticleCloud particles(mesh, cloudName);
        Pout<< "Reread particles:" << particles.size() << endl;

        forAllConstIter(passiveParticleCloud, particles, iter)
        {
            Pout<< "    " << iter().position() << " cell:" << iter().cell()
                << " origProc:" << iter().origProc()
                << " origId:" << iter().origId()
                << endl;
        }
    }


    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
