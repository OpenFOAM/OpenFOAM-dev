/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

Class
    Foam::zoneGenerators::face

Description
    A zoneGenerator which converts the point, cell and face zones
    from a list of zoneGenerators into a faceZone

    If a faceZone with the same name as the faceZone to be generated is
    returned by a zoneGenerator in the list it is ignored.

    The returned faceZone is oriented and provides the corresponding flipMap.

Usage
    \table
        Property     | Description             | Required | Default value
        type         | Type: face              | yes      |
        name         | Name of the zone        | no       | zoneGenerator name
        cellFaces    | Cell face selection option | no    | all
        moveUpdate   | Switch to update after mesh motion | no | false
    \endtable

    cellFaces face selection controls:
    \table
        Option        | Description
        all           | Select all faces of the cells in the cellZone
        inner         | Select the faces inside the cellZone
        outer         | Select the faces on the outer surface of the cellZone
        outerInternal | Select the internal faces on the outer surface \\
                        of the cellZone
    \endtable

    To convert cellZone named \c cone1 generated by
    \verbatim
        cone1
        {
            type        truncatedCone;
            zoneType    cell;

            point1      (-0.0075 0 -1);
            point2      (-0.003 0.0025 1);
            radius1     0.001;
            radius2     0.0005;
        }
    \endverbatim
    into a faceZone also named \c cone1, the following zoneGenerator named \c
    cone1Faces can be used
    \verbatim
        cone1Faces
        {
            type        face;
            name        cone1;

            cone1;
        }
    \endverbatim
    To create a faceZone which includes only the internal mesh faces on the
    outer surface of a pre-generated cellZone named \c leftFluid the \c
    outerInternal \c cellFaces face selection option can be used:
    \verbatim
        cycLeft
        {
            type        face;
            cellFaces   outerInternal;
            leftFluid;
        }
    \endverbatim

SourceFiles
    face_zoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef face_zoneGenerator_H
#define face_zoneGenerator_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class face Declaration
\*---------------------------------------------------------------------------*/

class face
:
    public zoneGenerator
{
public:

        //- Enumeration defining the valid options
        enum class cellFaces
        {
            all,
            inner,
            outer,
            outerInternal
        };

        //- Names of the valid options
        static const NamedEnum<cellFaces, 4> cellFacesNames;


private:

    // Private Data

        //- List of zoneGenerators
        //  which generate the zones to convert to the faceZone
        zoneGeneratorList zoneGenerators_;

        //- Cell face selection option
        cellFaces cellFaces_;


public:

    //- Runtime type information
    TypeName("face");


    // Constructors

        //- Construct from dictionary
        face
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~face();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
