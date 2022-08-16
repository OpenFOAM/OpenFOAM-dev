/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::DeltaInfoData

Description
    Class to be used with FaceCellWave which calculates the delta (thickness)
    across a layered mesh

SourceFiles
    DeltaInfoDataI.H

\*---------------------------------------------------------------------------*/

#ifndef DeltaInfoData_H
#define DeltaInfoData_H

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
template<class Type>
class DeltaInfoData;
template<class Type>
Ostream& operator<<(Ostream&, const DeltaInfoData<Type>&);
template<class Type>
Istream& operator>>(Istream&, DeltaInfoData<Type>&);

/*---------------------------------------------------------------------------*\
                          Class DeltaInfoData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DeltaInfoData
{
    // Private Data

        //- Direction of propagation (+1 or -1) relative to the face normal.
        //  Only valid when on a face. Takes a value of 0 when in a cell.
        label direction_;

        //- The previous face index. Only valid when in a cell. Takes a value
        //  of -1 if on a face.
        label prevFace_;

        //- Data
        Type data_;


public:

    // Constructors

        //- Construct null
        inline DeltaInfoData();

        //- Construct for a face given a direction and data
        inline DeltaInfoData(const label direction, const Type& data);


    // Member Functions

        // Access

            //- Return the data
            inline const Type& data() const;


        // Needed by FaceCellWave

            //- Check whether the DeltaInfoData has been changed at all or still
            //  contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for checking
            //  consistency across cyclics.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const DeltaInfoData<Type>&,
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
                const label neighbourFacei,
                const DeltaInfoData<Type>& neighbourInfo,
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
                const DeltaInfoData<Type>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const DeltaInfoData<Type>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test equality
            template<class TrackingData>
            inline bool equal
            (
                const DeltaInfoData<Type>&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const DeltaInfoData<Type>&) const;
        inline bool operator!=(const DeltaInfoData<Type>&) const;


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const DeltaInfoData<Type>&
        );
        friend Istream& operator>> <Type>
        (
            Istream&,
            DeltaInfoData<Type>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DeltaInfoDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
