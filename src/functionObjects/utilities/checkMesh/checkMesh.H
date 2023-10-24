/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::functionObjects::checkMesh

Description
    Executes primitiveMesh::checkMesh(true) every execute time for which the
    mesh changed, i.e. moved or changed topology.

    Useful to check the correctness of changing and morphing meshes.

Usage
    \table
        Property      | Description               | Required     | Default value
        type          | type name: checkMesh      | yes          |
        noTopology    | Skip checking the mesh topology | no     | false
        allTopology   | Check all addressing      | no           | false
        allGeometry   | Check all geometry        | no           | false
        writeSurfaces | Reconstruct and write problem faces | no | false
        surfaceFormat | Format for problem faceSets | no         | vtk
        writeSets     | Reconstruct and write problem points | no | false
        setFormat     | Format used to write the problem pointSets | no | vtk
        nonOrthThreshold | Threshold for non-orthogonality errors | no | 70 deg
        skewThreshold | Threshold for reporting skewness errors | no | 4
    \endtable

    Example of checkMesh specification:
    \verbatim
    checkMesh
    {
        type            checkMesh;
        libs            ("libutilityFunctionObjects.so");

        executeControl  timeStep;
        executeInterval 10;

        allGeometry     true;
        allTopology     true;

        writeSurfaces   true;
        surfaceFormat   vtk;

        writeSets       true;
        setFormat       vtk;
    }
    \endverbatim

    or using the standard configuration file:

    \verbatim
    #includeFunc checkMesh(executeInterval=10, allGeometry=true)
    \endverbatim

SourceFiles
    checkMesh.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_checkMesh_H
#define functionObjects_checkMesh_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class checkMesh Declaration
\*---------------------------------------------------------------------------*/

class checkMesh
:
    public fvMeshFunctionObject
{
    // Private Data

        bool noTopology_;
        bool allGeometry_;
        bool allTopology_;

        bool writeSurfaces_;
        word surfaceFormat_;

        bool writeSets_;
        word setFormat_;

        scalar nonOrthThreshold_;
        scalar skewThreshold_;


public:

    //- Runtime type information
    TypeName("checkMesh");


    // Constructors

        //- Construct from Time and dictionary
        checkMesh
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~checkMesh();


    // Member Functions

        //- Read the checkMesh controls
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute primitiveMesh::checkMesh(true)
        virtual bool execute();

        //- Do nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
