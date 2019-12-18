/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::extendedFeatureEdgeMesh

Description

    extendedEdgeMesh + IO.

SourceFiles
    extendedFeatureEdgeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef extendedFeatureEdgeMesh_H
#define extendedFeatureEdgeMesh_H

#include "extendedEdgeMesh.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;

/*---------------------------------------------------------------------------*\
                       Class extendedFeatureEdgeMesh Declaration
\*---------------------------------------------------------------------------*/

class extendedFeatureEdgeMesh
:
    public regIOobject,
    public extendedEdgeMesh
{

public:

    //- Runtime type information
    TypeName("extendedFeatureEdgeMesh");


    // Constructors

        //- Construct (read) given an IOobject
        extendedFeatureEdgeMesh(const IOobject&);

        //- Copy constructor
        extendedFeatureEdgeMesh
        (
            const IOobject&,
            const extendedEdgeMesh&
        );

        //- Construct given a surface with selected edges,point
        //  (surfaceFeatures), an objectRegistry and a
        //  fileName to write to.
        //  Extracts, classifies and reorders the data from surfaceFeatures.
        extendedFeatureEdgeMesh
        (
            const surfaceFeatures& sFeat,
            const objectRegistry& obr,
            const fileName& sFeatFileName,
            const boolList& surfBaffleRegions
        );

        //- Construct from PrimitivePatch
        extendedFeatureEdgeMesh
        (
            const IOobject&,
            const PrimitivePatch<faceList, pointField>& surf,
            const labelList& featureEdges,
            const labelList& regionFeatureEdges,
            const labelList& featurePoints
        );

        //- Construct from all components
        extendedFeatureEdgeMesh
        (
            const IOobject& io,
            const pointField& pts,
            const edgeList& eds,
            label concaveStart,
            label mixedStart,
            label nonFeatureStart,
            label internalStart,
            label flatStart,
            label openStart,
            label multipleStart,
            const vectorField& normals,
            const List<sideVolumeType>& normalVolumeTypes,
            const vectorField& edgeDirections,
            const labelListList& normalDirections,
            const labelListList& edgeNormals,
            const labelListList& featurePointNormals,
            const labelListList& featurePointEdges,
            const labelList& regionEdges
        );


    //- Destructor
    virtual ~extendedFeatureEdgeMesh();


    // IO

        //- Give precedence to the regIOobject write
        using regIOobject::write;

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<extendedFeatureEdgeMesh>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
