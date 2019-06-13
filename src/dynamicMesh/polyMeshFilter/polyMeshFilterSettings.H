/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::polyMeshFilterSettings

Description
    Class to store the settings for the polyMeshFilter class.

SourceFiles
    polyMeshFilterSettings.C
    polyMeshFilterSettingsI.H

\*---------------------------------------------------------------------------*/

#ifndef polyMeshFilterSettings_H
#define polyMeshFilterSettings_H

#include "dictionary.H"
#include "Switch.H"
#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class polyMeshFilterSettings Declaration
\*---------------------------------------------------------------------------*/

class polyMeshFilterSettings
{
    // Private Data

        //- Dictionary containing the coefficient sub-dictionaries
        const dictionary dict_;

        //- After collapsing, check the mesh quality and redo the collapsing
        //  iteration if there are too many bad faces in the mesh
        Switch controlMeshQuality_;

        //- Coefficients for collapsing edges
        const dictionary& collapseEdgesCoeffDict_;

        //- Coefficients for collapsing faces
        dictionary collapseFacesCoeffDict_;

        //- Coefficients for controlling the mesh quality
        dictionary meshQualityCoeffDict_;

        //- Remove edges shorter than this length
        const scalar minLen_;

        //- Merge points that are only attached to two edges and have an angle
        //  between the edge greater than this value
        const scalar maxCos_;

        //- The amount that the local minimum edge length will be reduced by if
        //  the edge is part of a collapse string that generates poor quality
        //  faces
        const scalar edgeReductionFactor_;

        //- Maximum number of outer iterations
        const label maxIterations_;

        //- Maximum number of smoothing iterations of minEdgeLen_ and
        //  faceFilterFactor_
        const label maxSmoothIters_;

        //- Initialisation value of faceFilterFactor_
        const scalar initialFaceLengthFactor_;

        //- The amount that the local face size factor will be reduced by if
        //  the face is part of a collapse string that generates poor quality
        //  faces
        const scalar faceReductionFactor_;

        //- Maximum number of times a deleted point can be associated with the
        //  creation of a bad face it is forced to be kept.
        const label maxPointErrorCount_;


public:

    // Constructors

        //- Construct from dictionary
        explicit polyMeshFilterSettings(const dictionary& dict);


    //- Destructor
    ~polyMeshFilterSettings(){};


    // Member Functions

        // Access

            inline const dictionary& collapseEdgesCoeffDict() const;

            inline const dictionary& collapseFacesCoeffDict() const;

            inline const dictionary& meshQualityCoeffDict() const;


            inline const Switch& controlMeshQuality() const;


            inline const scalar& minLen() const;

            inline const scalar& maxCos() const;

            inline const scalar& edgeReductionFactor() const;

            inline const label& maxIterations() const;

            inline const label& maxSmoothIters() const;

            inline const scalar& initialFaceLengthFactor() const;

            inline const scalar& faceReductionFactor() const;

            inline const label& maxPointErrorCount() const;


        // Write

            //- Write the settings to a stream
            void writeSettings(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshFilterSettingsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
