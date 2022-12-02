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
    Foam::vtkWriteFieldOps

Description
    VTK ASCII and binary write functions

SourceFiles
    vtkWriteFieldOps.C

\*---------------------------------------------------------------------------*/

#ifndef vtkWriteFieldOps_H
#define vtkWriteFieldOps_H

#include "vtkWriteOps.H"
#include "pointFieldsFwd.H"
#include "vtkMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Namespace vtkWriteOps Declaration
\*---------------------------------------------------------------------------*/

namespace vtkWriteOps
{
    //- Write volField with cell values (including decomposed cells)
    template<class Type>
    void write
    (
        std::ostream&,
        const bool binary,
        const DimensionedField<Type, volMesh>&,
        const vtkMesh&
    );

    //- Write pointField on all mesh points. Interpolate to cell centre
    //  for decomposed cell centres.
    template<class Type>
    void write
    (
        std::ostream&,
        const bool binary,
        const PointField<Type>&,
        const vtkMesh&
    );

    //- Write interpolated field on points and original cell values on
    //  decomposed cell centres.
    template<class Type>
    void write
    (
        std::ostream&,
        const bool binary,
        const VolField<Type>&,
        const PointField<Type>&,
        const vtkMesh&
    );

    //- Write generic GeometricFields
    template<class Type, template<class> class PatchField, class GeoMesh>
    void write
    (
        std::ostream&,
        const bool binary,
        const PtrList<GeometricField<Type, PatchField, GeoMesh>>&,
        const vtkMesh&
    );

    //- Interpolate and write volFields
    template<class Type>
    void write
    (
        std::ostream&,
        const bool binary,
        const volPointInterpolation&,
        const PtrList<VolField<Type>>&,
        const vtkMesh&
    );

} // End namespace vtkWriteOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vtkWriteFieldOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
