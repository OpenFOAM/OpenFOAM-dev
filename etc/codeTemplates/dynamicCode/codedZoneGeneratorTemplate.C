/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "codedZoneGeneratorTemplate.H"
#include "volFields.H"
#include "read.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(${typeName}ZoneGenerator, 0);

addRemovableToRunTimeSelectionTable
(
    zoneGenerator,
    ${typeName}ZoneGenerator,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // Unique function name that can be checked
    // to ensure the correct library version has been loaded
    void ${uniqueFunctionName}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}ZoneGenerator::${typeName}ZoneGenerator
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}ZoneGenerator::~${typeName}ZoneGenerator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet ${typeName}ZoneGenerator::generate() const
{
    if (${verbose})
    {
        Info<<"fields ${typeName} sha1: ${SHA1sum}\n";
    }

    enabledLabelList pointIndices;
    enabledLabelList cellIndices;
    enabledLabelList faceIndices;
    boolList flipMap;

//{{{ begin code
    ${code}
//}}} end code

    return zoneSet
    (
        pointIndices.enabled()
      ? new pointZone
        (
            zoneName_,
            pointIndices,
            mesh_.pointZones(),
            moveUpdate_,
            true
        )
      : nullptr,
        cellIndices.enabled()
      ? new cellZone
        (
            zoneName_,
            cellIndices,
            mesh_.cellZones(),
            moveUpdate_,
            true
        )
      : nullptr,
        faceIndices.enabled()
      ? flipMap.size()
          ? new faceZone
            (
                zoneName_,
                faceIndices,
                flipMap,
                mesh_.faceZones(),
                moveUpdate_,
                true
            )
          : new faceZone
            (
                zoneName_,
                faceIndices,
                mesh_.faceZones(),
                moveUpdate_,
                true
            )
      : nullptr
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
