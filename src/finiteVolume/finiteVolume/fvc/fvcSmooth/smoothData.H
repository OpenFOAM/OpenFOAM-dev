/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::smoothData

Description
    Helper class used by the fvc::smooth and fvc::spread functions.

SourceFiles
    smoothData.H
    smoothDataI.H

\*---------------------------------------------------------------------------*/

#ifndef smoothData_H
#define smoothData_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class smoothData Declaration
\*---------------------------------------------------------------------------*/

class smoothData
{

public:

    //- Class used to pass additional data in
    class trackData
    {
    public:

        //- Cut off distance
        scalar maxRatio;
    };


private:

    scalar value_;

    // Private Member Functions

        //- Update - gets information from neighbouring face/cell and
        //  uses this to update itself (if necessary) and return true
        template<class TrackingData>
        inline bool update
        (
            const smoothData& svf,
            const scalar scale,
            const scalar tol,
            TrackingData& td
        );


public:


    // Constructors

        //- Construct null
        inline smoothData();

        //- Construct from inverse field value
        inline smoothData(const scalar value);


    // Member Functions

        // Access

            //- Return value
            scalar value() const
            {
                return value_;
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
                const smoothData&,
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
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
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
                const label thisCelli,
                const label neighbourFacei,
                const smoothData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const smoothData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const smoothData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const smoothData&, TrackingData& td) const;


        // Member Operators

            inline void operator=(const scalar value);

            // Needed for List IO
            inline bool operator==(const smoothData&) const;

            inline bool operator!=(const smoothData&) const;


        // IOstream Operators

            friend Ostream& operator<<
            (
                Ostream& os,
                const smoothData& svf
            )
            {
                return os  << svf.value_;
            }

            friend Istream& operator>>(Istream& is, smoothData& svf)
            {
                return is  >> svf.value_;
            }
};


//- Data associated with smoothData type are contiguous
template<>
inline bool contiguous<smoothData>()
{
    return true;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "smoothDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
