/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::faceCoupleInfo

Description
    Container for information needed to couple to meshes. When constructed
    from two meshes and a list of coupled faces returns the mapping between
    points.

SourceFiles
    faceCoupleInfo.C

\*---------------------------------------------------------------------------*/

#ifndef faceCoupleInfo_H
#define faceCoupleInfo_H

#include "pointField.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class faceCoupleInfo Declaration
\*---------------------------------------------------------------------------*/

class faceCoupleInfo
{
    // Private Data

        //- Master patch
        autoPtr<indirectPrimitivePatch> masterPatchPtr_;

        //- Slave patch
        autoPtr<indirectPrimitivePatch> slavePatchPtr_;

        //- The number of coupled points
        label nCouplePoints_;

        //- Master-to-couple point addressing
        labelList masterToCouplePoints_;

        //- Slave-to-couple point addressing
        labelList slaveToCouplePoints_;


    // Private Member Functions

        // Debugging

            //- Calculate face centres from (subset of) faces.
            template<template<class> class FaceList>
            static pointField calcFaceCentres
            (
                const FaceList<face>&,
                const pointField&,
                const label start,
                const label size
            );

            //- Calculate face point averages from (subset of) faces.
            template<template<class> class FaceList>
            static pointField calcFacePointAverages
            (
                const FaceList<face>&,
                const pointField&,
                const label start,
                const label size
            );

            //- Write edges
            static void writeOBJ
            (
                const fileName& fName,
                const edgeList& edges,
                const pointField& points
            );

            //- Write edges
            static void writeOBJ
            (
                const fileName& fName,
                const pointField& points0,
                const pointField& points1
            );

            //- Write connections between corresponding points and faces
            //  as .obj files.
            void writePointsFaces() const;


public:

    //- Runtime type information
    ClassName("faceCoupleInfo");


    // Constructors

        //- Construct from meshes and subset of mesh faces (i.e.,
        //  indirectPrimitivePatch addressing). Faces should be supplied in
        //  order and the points within them should be ordered (i.e., coupled
        //  faces should have the same zero-point and opposite orientations).
        faceCoupleInfo
        (
            const polyMesh& masterMesh,
            const labelList& masterAddressing,
            const polyMesh& slaveMesh,
            const labelList& slaveAddressing
        );


    //- Destructor
    ~faceCoupleInfo();


    // Member Functions

        // Access

            //- Addressing engine for coupled faces on mesh0
            const indirectPrimitivePatch& masterPatch() const
            {
                return masterPatchPtr_();
            }

            //- Addressing engine for coupled faces on mesh1
            const indirectPrimitivePatch& slavePatch() const
            {
                return slavePatchPtr_();
            }

            //- Access master-to-couple point addressing
            const labelList& masterToCouplePoints() const
            {
                return masterToCouplePoints_;
            }

            //- Access slave-to-couple point addressing
            const labelList& slaveToCouplePoints() const
            {
                return slaveToCouplePoints_;
            }

            //- Generate couple-to-master point addressing
            labelListList coupleToMasterPoints() const
            {
                return invertOneToMany(nCouplePoints_, masterToCouplePoints_);
            }

            //- Generate couple-to-slave point addressing
            labelListList coupleToSlavePoints() const
            {
                return invertOneToMany(nCouplePoints_, slaveToCouplePoints_);
            }

            //- Return a coupled face
            const face coupleFace(const label coupleFacei) const
            {
                return
                    renumber
                    (
                        masterToCouplePoints(),
                        masterPatch().localFaces()[coupleFacei]
                    );
            }
 };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faceCoupleInfoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
