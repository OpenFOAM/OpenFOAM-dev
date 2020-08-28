/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::fv::steadyStateD2dt2Scheme

Description
    SteadyState d2dt2 which returns 0.

SourceFiles
    steadyStateD2dt2Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef steadyStateD2dt2Scheme_H
#define steadyStateD2dt2Scheme_H

#include "d2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class steadyStateD2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class steadyStateD2dt2Scheme
:
    public fv::d2dt2Scheme<Type>
{
public:

    //- Runtime type information
    TypeName("steadyState");


    // Constructors

        //- Construct from mesh
        steadyStateD2dt2Scheme(const fvMesh& mesh)
        :
            d2dt2Scheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        steadyStateD2dt2Scheme(const fvMesh& mesh, Istream& is)
        :
            d2dt2Scheme<Type>(mesh, is)
        {}

        //- Disallow default bitwise copy construction
        steadyStateD2dt2Scheme(const steadyStateD2dt2Scheme&) = delete;


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::d2dt2Scheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmD2dt2
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmD2dt2
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmD2dt2
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const steadyStateD2dt2Scheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "steadyStateD2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
