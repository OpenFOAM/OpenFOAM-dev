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
    Foam::pointMVCWeight

Description
    Container to calculate weights for interpolating directly from vertices
    of cell using Mean Value Coordinates.

    Based on (VTK's vtkMeanValueCoordinatesInterpolator's) implementation
    of "Spherical Barycentric Coordinates"
    2006 paper Eurographics Symposium on Geometry Processing
    by Torsten Langer, Alexander Belyaev and Hans-Peter Seide



SourceFiles
    pointMVCWeight.C

\*---------------------------------------------------------------------------*/

#ifndef pointMVCWeight_H
#define pointMVCWeight_H

#include "scalarField.H"
#include "vectorField.H"
#include "Map.H"
#include "DynamicList.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class pointMesh;
template<class T> class pointPatchField;
template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;
class face;

/*---------------------------------------------------------------------------*\
                           Class pointMVCWeight Declaration
\*---------------------------------------------------------------------------*/

class pointMVCWeight
{
protected:

    // Protected data

        //- Cell index
        const label cellIndex_;

        //- Weights applied to cell vertices
        scalarField weights_;


    // Protected Member Functions

        //- Calculate weights from single face's vertices only
        void calcWeights
        (
            const Map<label>& toLocal,
            const face& f,
            const DynamicList<point>& u,
            const scalarField& dist,
            scalarField& weights
        ) const;

        //- Calculate weights from all cell's vertices
        void calcWeights
        (
            const polyMesh& mesh,
            const labelList& toGlobal,
            const Map<label>& toLocal,
            const vector& position,
            const vectorField& uVec,
            const scalarField& dist,
            scalarField& weights
        ) const;

public:

    //- Debug switch
    static int debug;

    //- Tolerance used in calculating barycentric co-ordinates
    //  (applied to normalised values)
    static scalar tol;


    // Constructors

        //- Construct from components
        pointMVCWeight
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label facei = -1
        );


    // Member Functions

        //- Cell index
        inline label cell() const
        {
            return cellIndex_;
        }

        //- Interpolation weights (in order of cellPoints)
        inline const scalarField& weights() const
        {
            return weights_;
        }

        //- Interpolate field
        template<class Type>
        inline Type interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>& psip
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointMVCWeightI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
