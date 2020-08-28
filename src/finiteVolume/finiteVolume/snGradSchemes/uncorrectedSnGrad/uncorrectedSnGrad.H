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
    Foam::fv::uncorrectedSnGrad

Description
    Simple central-difference snGrad scheme using the non-orthogonal mesh
    delta-coefficients but without non-orthogonal correction.

SourceFiles
    uncorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef uncorrectedSnGrad_H
#define uncorrectedSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class uncorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uncorrectedSnGrad
:
    public snGradScheme<Type>
{
public:

    //- Runtime type information
    TypeName("uncorrected");


    // Constructors

        //- Construct from mesh
        uncorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}


        //- Construct from mesh and data stream
        uncorrectedSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~uncorrectedSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the uncorrectedSnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uncorrectedSnGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uncorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
