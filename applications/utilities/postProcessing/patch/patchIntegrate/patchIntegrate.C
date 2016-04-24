/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    patchIntegrate

Description
    Calculates the integral of the specified field over the specified patch.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FieldType>
void printIntegrate
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const label patchi,
    bool& done
)
{
    if (!done && fieldHeader.headerClassName() == FieldType::typeName)
    {
        Info<< "    Reading " << fieldHeader.headerClassName() << " "
            << fieldHeader.name() << endl;

        FieldType field(fieldHeader, mesh);

        Info<< "    Integral of " << fieldHeader.name()
            << " over vector area of patch "
            << mesh.boundary()[patchi].name() << '[' << patchi << ']' << " = "
            << gSum
               (
                   mesh.Sf().boundaryField()[patchi]
                  *field.boundaryField()[patchi]
               )
            << nl;

        Info<< "    Integral of " << fieldHeader.name()
            << " over area magnitude of patch "
            << mesh.boundary()[patchi].name() << '[' << patchi << ']' << " = "
            << gSum
               (
                   mesh.magSf().boundaryField()[patchi]
                  *field.boundaryField()[patchi]
               )
            << nl;

        done = true;
    }
}


template<class FieldType>
void printSum
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const label patchi,
    bool& done
)
{
    if (!done && fieldHeader.headerClassName() == FieldType::typeName)
    {
        Info<< "    Reading " << FieldType::typeName << " "
            << fieldHeader.name() << endl;

        FieldType field(fieldHeader, mesh);
        typename FieldType::value_type sumField = gSum
        (
            field.boundaryField()[patchi]
        );

        Info<< "    Integral of " << fieldHeader.name() << " over patch "
            << mesh.boundary()[patchi].name() << '[' << patchi << ']' << " = "
            << sumField << nl;

        done = true;
    }
}



int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const word fieldName = args[1];
    const word patchName = args[2];

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check field exists
        if (fieldHeader.headerOk())
        {
            mesh.readUpdate();

            const label patchi = mesh.boundaryMesh().findPatchID(patchName);
            if (patchi < 0)
            {
                FatalError
                    << "Unable to find patch " << patchName << nl
                    << exit(FatalError);
            }

            // Give patch area
            Info<< "    Area vector of patch "
                << patchName << '[' << patchi << ']' << " = "
                << gSum(mesh.Sf().boundaryField()[patchi]) << endl;
            Info<< "    Area magnitude of patch "
                << patchName << '[' << patchi << ']' << " = "
                << gSum(mesh.magSf().boundaryField()[patchi]) << endl;

            // Read field and calc integral
            bool done = false;
            printIntegrate<volScalarField>
            (
                mesh,
                fieldHeader,
                patchi,
                done
            );
            printIntegrate<volVectorField>
            (
                mesh,
                fieldHeader,
                patchi,
                done
            );

            //- No tensor integrations
            //printIntegrate<volSphericalTensorField>
            //(
            //    mesh,
            //    fieldHeader,
            //    patchi,
            //    done
            //);
            //printIntegrate<volSymmTensorField>
            //(
            //    mesh,
            //    fieldHeader,
            //    patchi,
            //    done
            //);
            //printIntegrate<volTensorField>
            //(
            //    mesh,
            //    fieldHeader,
            //    patchi,
            //    done
            //);

            printSum<surfaceScalarField>
            (
                mesh,
                fieldHeader,
                patchi,
                done
            );
            printSum<surfaceVectorField>
            (
                mesh,
                fieldHeader,
                patchi,
                done
            );
            printSum<volSphericalTensorField>
            (
                mesh,
                fieldHeader,
                patchi,
                done
            );
            printSum<volSymmTensorField>
            (
                mesh,
                fieldHeader,
                patchi,
                done
            );
            printSum<volTensorField>
            (
                mesh,
                fieldHeader,
                patchi,
                done
            );

            if (!done)
            {
                FatalError
                    << "Only possible to integrate "
                    << "volFields and surfaceFields."
                    << " Field " << fieldName << " is of type "
                    << fieldHeader.headerClassName()
                    << nl << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No field " << fieldName << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
