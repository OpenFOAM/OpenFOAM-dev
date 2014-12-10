/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    smapToFoam

Description
    Translates a STAR-CD SMAP data file into OpenFOAM field format.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("SMAP fileName");

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

#   include "createTime.H"

    fileNameList fieldNames = readDir(runTime.timePath(), fileName::FILE);
    dictionary fieldNameDict;
    forAll(fieldNames, i)
    {
        fieldNameDict.add(fieldNames[i], word(fieldNames[i]));
    }

    dictionary nameMap;
    if (fieldNameDict.found("U")) nameMap.add("SU", word("U"));
    if (fieldNameDict.found("p")) nameMap.add("P", word("p"));
    if (fieldNameDict.found("T")) nameMap.add("T", word("T"));
    if (fieldNameDict.found("rho")) nameMap.add("DENS", word("rho"));
    if (fieldNameDict.found("k")) nameMap.add("TE", word("k"));
    if (fieldNameDict.found("epsilon")) nameMap.add("ED", word("epsilon"));
    if (fieldNameDict.found("nuEff")) nameMap.add("VIS", word("nuEff"));

#   include "createMesh.H"

    IFstream smapFile(args[1]);

    if (!smapFile.good())
    {
        FatalErrorIn(args.executable())
            << "Cannot open SMAP file " << smapFile.name()
            << exit(FatalError);
    }

    while (!smapFile.eof())
    {
        wordList starFieldNames(10);

        token fieldName(smapFile);

        if (!smapFile.good())
        {
            break;
        }

        if
        (
            fieldName.type() != token::WORD
         && fieldName.wordToken() != "CELL"
        )
        {
            FatalErrorIn(args.executable())
                << "Expected first CELL, found "
                << fieldName
                << exit(FatalError);
        }

        label nCols = 0;
        smapFile >> fieldName;
        while (fieldName.type() == token::WORD)
        {
            starFieldNames[nCols++] = fieldName.wordToken();
            smapFile >> fieldName;
        }

        List<volScalarField*> sFields
        (
            nCols,
            reinterpret_cast<volScalarField*>(0)
        );

        List<volVectorField*> vFields
        (
            nCols,
            reinterpret_cast<volVectorField*>(0)
        );

        label i=0;
        while (i < nCols)
        {
            if (nameMap.found(starFieldNames[i]))
            {
                if (starFieldNames[i] == "SU")
                {
                    vFields[i] =
                    new volVectorField
                    (
                        IOobject
                        (
                            nameMap.lookup(starFieldNames[i]),
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh
                    );

                    i += 3;
                }
                else
                {
                    sFields[i] =
                    new volScalarField
                    (
                        IOobject
                        (
                            nameMap.lookup(starFieldNames[i]),
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh
                    );

                    i++;
                }
            }
            else
            {
                i++;
            }
        }


        label cell;
        scalar value;
        forAll(mesh.cells(), celli)
        {
            if (celli > 0)
            {
                smapFile >> cell;
            }

            label i=0;
            while (i < nCols)
            {
                if (sFields[i])
                {
                    smapFile >> (*sFields[i])[celli];
                    i++;
                }
                else if (vFields[i])
                {
                    smapFile >> (*vFields[i])[celli].x();
                    smapFile >> (*vFields[i])[celli].y();
                    smapFile >> (*vFields[i])[celli].z();
                    i += 3;
                }
                else
                {
                    smapFile >> value;
                    i++;
                }
            }
        }

        for (label i=0; i<nCols; i++)
        {
            if (sFields[i])
            {
                sFields[i]->correctBoundaryConditions();
                sFields[i]->write();
                delete sFields[i];
                sFields[i] = NULL;
            }
            else if (vFields[i])
            {
                vFields[i]->correctBoundaryConditions();
                vFields[i]->write();
                delete vFields[i];
                vFields[i] = NULL;
            }
        }

        // Read dummy entry and check the cell index
        smapFile >> cell;

        if (cell != 0)
        {
            FatalErrorIn(args.executable())
                << "Expected first SMAP dummy entry to be cell 0, found "
                << cell
                << exit(FatalError);
        }

        for (label i=0; i<nCols; i++)
        {
            smapFile >> value;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
