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

\*---------------------------------------------------------------------------*/

#include "randomise.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(randomise, 0);
        addToRunTimeSelectionTable(calcType, randomise, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::randomise::randomise()
:
    calcType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::randomise::~randomise()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::randomise::init()
{
    argList::validArgs.append("randomise");
    argList::validArgs.append("perturbation");
    argList::validArgs.append("fieldName");
}


void Foam::calcTypes::randomise::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{}


void Foam::calcTypes::randomise::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    const scalar pertMag = args.argRead<scalar>(2);
    const word fieldName = args[3];

    Random rand(1234567);

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
        bool processed = false;

        writeRandomField<vector>
        (
            fieldHeader,
            pertMag,
            rand,
            mesh,
            processed
        );
        writeRandomField<sphericalTensor>
        (
            fieldHeader,
            pertMag,
            rand,
            mesh,
            processed
        );
        writeRandomField<symmTensor>
        (
            fieldHeader,
            pertMag,
            rand,
            mesh,
            processed
        );
        writeRandomField<tensor>
        (
            fieldHeader,
            pertMag,
            rand,
            mesh,
            processed
        );

        if (!processed)
        {
            FatalError
                << "Unable to process " << fieldName << nl
                << "No call to randomise for fields of type "
                << fieldHeader.headerClassName() << nl << nl
                << exit(FatalError);
        }
    }
    else
    {
        Info<< "    No " << fieldName << endl;
    }
}


// ************************************************************************* //
