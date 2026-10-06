/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::topoDistance

Description
    For use with FaceCellWave. Determines topological distance to starting faces

SourceFiles
    topoDistanceI.H
    topoDistance.C

\*---------------------------------------------------------------------------*/

#ifndef topoDistance_H
#define topoDistance_H

#include "scalar.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;

// Forward declaration of friend functions and operators
class topoDistance;
inline Istream& operator>>(Istream&, topoDistance&);
inline Ostream& operator<<(Ostream&, const topoDistance&);

/*---------------------------------------------------------------------------*\
                      Class topoDistance Declaration
\*---------------------------------------------------------------------------*/

class topoDistance
{
    // Private Data

        //- Distance
        label distance_;


public:

    // Constructors

        //- Construct null
        inline topoDistance();

        //- Construct from distance
        inline topoDistance(const label distance);


    // Member Functions

        // Access

            //- Distance
            inline label distance() const
            {
                return distance_;
            }


        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const topoDistance&,
                const scalar,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class Patch, class TrackingData>
            inline void transform
            (
                const Patch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            template<class Mesh, class NeighbourFacei, class TrackingData>
            inline bool updateCell
            (
                const Mesh&,
                const label thisCelli,
                const NeighbourFacei& neighbourFacei,
                const topoDistance& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            template<class Mesh, class ThisFacei, class TrackingData>
            inline bool updateFace
            (
                const Mesh&,
                const ThisFacei& thisFacei,
                const label neighbourCelli,
                const topoDistance& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class Mesh, class ThisFacei, class TrackingData>
            inline bool updateFace
            (
                const Mesh&,
                const ThisFacei& thisFacei,
                const topoDistance& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const topoDistance&, TrackingData& td) const;


    // Member Operators

        inline bool operator==(const topoDistance&) const;
        inline bool operator!=(const topoDistance&) const;


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const topoDistance&);
        inline friend Istream& operator>>(Istream&, topoDistance&);
};


//- Data associated with topoDistance type are contiguous
template<>
inline bool contiguous<topoDistance>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "topoDistanceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
