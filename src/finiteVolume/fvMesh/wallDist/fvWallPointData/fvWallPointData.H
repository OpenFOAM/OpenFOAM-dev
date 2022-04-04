/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::fvWallPointData

Description
    Holds information (coordinate and normal) regarding nearest wall point.

    Is like fvWallPoint but transfer extra (passive) data.
    Used e.g. in wall distance calculation with wall reflection vectors.

SourceFiles
    fvWallPointDataI.H
    fvWallPointData.C

\*---------------------------------------------------------------------------*/

#ifndef fvWallPointData_H
#define fvWallPointData_H

#include "fvWallPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvWallPointData;

// Forward declaration of friend functions and operators

template<class Type>
Istream& operator>>(Istream&, fvWallPointData<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const fvWallPointData<Type>&);


/*---------------------------------------------------------------------------*\
                        Class fvWallPointData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvWallPointData
:
    public fvWallPoint
{
    // Private Data

        //- Data at nearest wall center
        Type data_;


    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const fvWallPointData<Type>& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    typedef Type dataType;


    // Constructors

        //- Construct null
        inline fvWallPointData();

        //- Construct from origin, normal, distance
        inline fvWallPointData
        (
            const point& origin,
            const Type& data,
            const scalar distSqr
        );


    // Member Functions

        // Access

            inline const Type& data() const;

            inline Type& data();

            template<class TrackingData>
            inline const Type& data(TrackingData& td) const;


        // Needed by meshWave

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const fvPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of celli
            template<class TrackingData>
            inline bool updateCell
            (
                const fvMesh& mesh,
                const label thisCelli,
                const labelPair& neighbourPatchAndFacei,
                const fvWallPointData<Type>& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of facei
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const label neighbourCelli,
                const fvWallPointData<Type>& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of facei
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const fvWallPointData<Type>& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const fvWallPointData<Type>&
        );

        friend Istream& operator>> <Type>
        (
            Istream&,
            fvWallPointData<Type>&
        );
};


//- Data associated with fvWallPointData type are contiguous. List the usual
//  ones.

template<>
inline bool contiguous<fvWallPointData<bool>>()
{
    return contiguous<fvWallPoint>();
}
template<>
inline bool contiguous<fvWallPointData<label>>()
{
    return contiguous<fvWallPoint>();
}
template<>
inline bool contiguous<fvWallPointData<scalar>>()
{
    return contiguous<fvWallPoint>();
}
template<>
inline bool contiguous<fvWallPointData<vector>>()
{
    return contiguous<fvWallPoint>();
}
template<>
inline bool contiguous<fvWallPointData<sphericalTensor>>()
{
    return contiguous<fvWallPoint>();
}
template<>
inline bool contiguous<fvWallPointData<symmTensor>>()
{
    return contiguous<fvWallPoint>();
}
template<>
inline bool contiguous<fvWallPointData<tensor>>()
{
    return contiguous<fvWallPoint>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvWallPointData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvWallPointDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
