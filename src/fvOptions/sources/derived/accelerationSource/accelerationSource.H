/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::fv::accelerationSource

Description
    This fvOption applies an explicit acceleration force to components of the
    velocity field.

Usage
    Example usage:
    \verbatim
    accelerationSource
    {
        type        accelerationSource;
        active      on;
        selectionMode all;
        U           U;
        velocity    scale;
        value       (-2.572 0 0);
        scale
        {
            type        halfCosineRamp;
            start       0;
            duration    10;
        }
    }
    \endverbatim

SourceFiles
    accelerationSource.C

\*---------------------------------------------------------------------------*/

#ifndef accelerationSource_H
#define accelerationSource_H

#include "cellSetOption.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class accelerationSource Declaration
\*---------------------------------------------------------------------------*/

class accelerationSource
:
    public cellSetOption
{
    // Private Data

        //- Time-varying velocity
        autoPtr<Function1<vector>> velocity_;


    // Private Member Functions

        //- Source term to momentum equation
        template<class AlphaRhoFieldType>
        void add
        (
            const AlphaRhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


public:

    //- Runtime type information
    TypeName("accelerationSource");


    // Constructors

        //- Construct from components
        accelerationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~accelerationSource()
    {}


    // Member Functions

         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "accelerationSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
