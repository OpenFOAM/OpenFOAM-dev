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
    Foam::layerInfo

Description
    Class to be used with FaceCellWave which enumerates layers of cells

SourceFiles
    layerInfoI.H

\*---------------------------------------------------------------------------*/

#ifndef layerInfo_H
#define layerInfo_H

#include "pointField.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;

// Forward declaration of friend functions and operators
class layerInfo;
inline Ostream& operator<<(Ostream&, const layerInfo&);
inline Istream& operator>>(Istream&, layerInfo&);

/*---------------------------------------------------------------------------*\
                          Class layerInfo Declaration
\*---------------------------------------------------------------------------*/

class layerInfo
{
    // Private Data

        //- Layer index
        label layer_;

        //- Direction of propagation (+1 or -1) relative to the face normal.
        //  Only valid when on a face. Takes a value of 0 when in a cell.
        label direction_;

        //- The previous face index. Only valid when in a cell. Takes a value
        //  of -1 if on a face.
        label prevFace_;


    // Private Member Functions

        //- Generate an error in the event that opposing waves collided
        void collide() const;


public:

    // Constructors

        //- Construct null
        inline layerInfo();

        //- Construct given the face layer index and direction
        inline layerInfo(const label faceLayer, const label direction);


    // Member Functions

        // Access

            //- Return the face layer index
            inline label faceLayer() const;

            //- Return the cell layer index
            inline label cellLayer() const;


        // Needed by FaceCellWave

            //- Check whether the layerInfo has been changed at all or still
            //  contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for checking
            //  consistency across cyclics.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const layerInfo&,
                const scalar,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label nbrFacei,
                const layerInfo& nbrInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label nbrCelli,
                const layerInfo& nbrInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const layerInfo& nbrInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test equality
            template<class TrackingData>
            inline bool equal
            (
                const layerInfo&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const layerInfo&) const;
        inline bool operator!=(const layerInfo&) const;


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const layerInfo&);
        inline friend Istream& operator>>(Istream&, layerInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
inline bool contiguous<layerInfo>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "layerInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
