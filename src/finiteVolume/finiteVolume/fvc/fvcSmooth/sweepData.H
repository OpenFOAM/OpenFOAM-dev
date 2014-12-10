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

Class
    Foam::sweepData

Description
    Helper class used by fvc::sweep function.

SourceFiles
    sweepData.H
    sweepDataI.H

\*---------------------------------------------------------------------------*/

#ifndef sweepData_H
#define sweepData_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sweepData Declaration
\*---------------------------------------------------------------------------*/

class sweepData
{
    scalar value_;
    point origin_;

    // Private Member Functions

        //- Update - gets information from neighbouring face/cell and
        //  uses this to update itself (if necessary) and return true
        template<class TrackingData>
        inline bool update
        (
            const sweepData& svf,
            const point& position,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline sweepData();

        //- Construct from component
        inline sweepData(const scalar value, const point& origin);


    // Member Functions

        // Access

            //- Return value
            scalar value() const
            {
                return value_;
            }

            //- Return origin
            const point& origin() const
            {
                return origin_;
            }


        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data
            //  Used for cyclics checking
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const sweepData&,
                const scalar,
                TrackingData& td
            ) const;

            //- Convert any absolute coordinates into relative to
            //  (patch)face centre
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFaceI,
                const point& faceCentre,
                TrackingData& td
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFaceI,
                const point& faceCentre,
                TrackingData& td
            );

            //- Apply rotation matrix to any coordinates
            template<class TrackingData>
            inline void transform
            (
                const polyMesh&,
                const tensor&,
                TrackingData& td
            );

            //- Influence of neighbouring face
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCellI,
                const label neighbourFaceI,
                const sweepData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const label neighbourCellI,
                const sweepData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const sweepData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const sweepData&, TrackingData& td) const;


        // Member Operators

            inline void operator=(const scalar value);

            // Needed for List IO
            inline bool operator==(const sweepData&) const;

            inline bool operator!=(const sweepData&) const;


        // IOstream Operators

            friend Ostream& operator<<
            (
                Ostream& os,
                const sweepData& svf
            )
            {
                return os << svf.value_ << svf.origin_;
            }

            friend Istream& operator>>(Istream& is, sweepData& svf)
            {
                return is >> svf.value_ >> svf.origin_;
            }
};


//- Data associated with sweepData type are contiguous
template<>
inline bool contiguous<sweepData>()
{
    return true;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sweepDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
