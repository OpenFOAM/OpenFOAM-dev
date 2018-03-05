/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "helpBoundary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace helpTypes
    {
        defineTypeNameAndDebug(helpBoundary, 0);
        addNamedToRunTimeSelectionTable
        (
            helpType,
            helpBoundary,
            dictionary,
            boundary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helpTypes::helpBoundary::helpBoundary()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helpTypes::helpBoundary::~helpBoundary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helpTypes::helpBoundary::init()
{
    helpType::init();

    argList::validArgs.append("boundary");

    argList::addOption
    (
        "field",
        "word",
        "list available conditions for field"
    );
    argList::addBoolOption
    (
        "constraint",
        "list constraint patches"
    );
    argList::addBoolOption
    (
        "fixedValue",
        "list fixed value patches (use with -field option)"
    );
}


void Foam::helpTypes::helpBoundary::execute
(
    const argList& args,
    const fvMesh& mesh
)
{
    bool abortVar(env("FOAM_ABORT"));
    if (abortVar)
    {
        FatalErrorInFunction
            << "Please unset FOAM_ABORT to use this utility"
            << exit(FatalError);
    }

    word condition(word::null);
    word fieldName(word::null);

    if (args.optionReadIfPresent("browse", condition))
    {
        // TODO: strip scoping info if present?
        // e.g. conditions with leading "compressible::" will not be found
        // ".*[fF]vPatchField.*" + className + ".*"
        displayDoc(condition, ".*[fF]vPatchField.*", false);
    }
    else if (args.optionFound("constraint"))
    {
        HashSet<word> constraintTypes(fvPatch::constraintTypes());
        Info<< "Constraint types:" << nl;
        forAllConstIter(HashSet<word>, constraintTypes, iter)
        {
            Info<< "    " << iter.key() << nl;
        }
        Info<< endl;
    }
    else if (args.optionReadIfPresent("field", fieldName))
    {
        IOobject fieldHeader
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check for any type of volField
        if (fieldHeader.typeHeaderOk<volScalarField>(false))
        {
            if (args.optionFound("fixedValue"))
            {
                fixedValueFieldConditions<scalar>(fieldHeader);
                fixedValueFieldConditions<vector>(fieldHeader);
                fixedValueFieldConditions<sphericalTensor>(fieldHeader);
                fixedValueFieldConditions<symmTensor>(fieldHeader);
                fixedValueFieldConditions<tensor>(fieldHeader);
            }
            else
            {
                (void)fieldConditions<scalar>(fieldHeader, true);
                (void)fieldConditions<vector>(fieldHeader, true);
                (void)fieldConditions<sphericalTensor>(fieldHeader, true);
                (void)fieldConditions<symmTensor>(fieldHeader, true);
                (void)fieldConditions<tensor>(fieldHeader, true);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unable to read field " << fieldName << exit(FatalError);
        }
    }
    else if (args.optionReadIfPresent("fixedValue", fieldName))
    {
        FatalErrorInFunction
            << "-field option must be specified when using the -fixedValue "
            << "option" << exit(FatalError);
    }
    else
    {
        // TODO: strip scoping info if present?
        // e.g. conditions with leading "compressible::" will not be found
        // ".*[fF]vPatchField.*" + className + ".*"
        displayDocOptions(".*[fF]vPatchField.*", false);
    }
}


// ************************************************************************* //
