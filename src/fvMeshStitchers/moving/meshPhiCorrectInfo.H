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
    Foam::meshPhiCorrectInfo

SourceFiles
    meshPhiCorrectInfoI.H

\*---------------------------------------------------------------------------*/

#ifndef meshPhiCorrectInfo_H
#define meshPhiCorrectInfo_H

#include "meshPhiPreCorrectInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvPatch;
class fvMesh;
class transformer;

// Forward declaration of friend functions and operators
class meshPhiCorrectInfo;
Ostream& operator<<(Ostream&, const meshPhiCorrectInfo&);
Istream& operator>>(Istream&, meshPhiCorrectInfo&);

/*---------------------------------------------------------------------------*\
                          Class meshPhiCorrectInfo Declaration
\*---------------------------------------------------------------------------*/

class meshPhiCorrectInfo
{
public:

    // Public Classes

        //- Enumeration to define the mesh shape the info applies to
        enum class shape
        {
            invalid,
            face,
            cell
        };


        //- Tracking data. Mostly just a collection of references to the
        //  pre-correct info.
        class trackData
        {
            public:

                const List<meshPhiPreCorrectInfo>& internalFacePci_;
                const List<List<meshPhiPreCorrectInfo>>& patchFacePci_;
                const List<meshPhiPreCorrectInfo>& cellPci_;
                const DimensionedField<scalar, volMesh>& dVdtError_;

                trackData
                (
                    const List<meshPhiPreCorrectInfo>& internalFacePci,
                    const List<List<meshPhiPreCorrectInfo>>& patchFacePci,
                    const List<meshPhiPreCorrectInfo>& cellPci,
                    const DimensionedField<scalar, volMesh>& dVdtError
                )
                :
                    internalFacePci_(internalFacePci),
                    patchFacePci_(patchFacePci),
                    cellPci_(cellPci),
                    dVdtError_(dVdtError)
                {}
        };


private:

    // Private Data

        //- The primitive shape that this info applies to
        shape shape_;

        //- Flux correction (if on a face) or volume conservation correction
        //  (if in a cell)
        scalar delta_;


    // Private Member Functions

        //- Throw an error if this is not a face
        void validateFace() const;

        //- Throw an error if this is not a cell
        void validateCell() const;


public:

    // Constructors

        //- Construct null
        inline meshPhiCorrectInfo();

        //- Construct for given shape
        inline meshPhiCorrectInfo(const shape& s);


    // Member Functions

        // Access

            //- Return the flux correction
            inline scalar deltaPhi() const;

            //- Access the flux correction
            inline scalar& deltaPhi();

            //- Return the volume change rate error correction
            inline scalar deltaDVdtError() const;

            //- Access the volume change rate error correction
            inline scalar& deltaDVdtError();


        // Needed by FaceCellWave

            //- Check whether the meshPhiCorrectInfo has been changed at all
            //  or still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for checking
            //  consistency across cyclics.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const fvMesh& mesh,
                const meshPhiCorrectInfo&,
                const scalar,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const fvPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face
            template<class TrackingData>
            inline bool updateCell
            (
                const fvMesh& mesh,
                const label thisCelli,
                const labelPair& neighbourPatchAndFacei,
                const meshPhiCorrectInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const label neighbourCelli,
                const meshPhiCorrectInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const meshPhiCorrectInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test equality
            template<class TrackingData>
            inline bool equal
            (
                const meshPhiCorrectInfo&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const meshPhiCorrectInfo&) const;
        inline bool operator!=(const meshPhiCorrectInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const meshPhiCorrectInfo&);
        friend Istream& operator>>(Istream&, meshPhiCorrectInfo&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshPhiCorrectInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
