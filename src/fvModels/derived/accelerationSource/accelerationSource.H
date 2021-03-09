/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
    This fvModel applies an explicit acceleration force to components of the
    velocity field.

Usage
    Example usage:
    \verbatim
    accelerationSource1
    {
        type            accelerationSource;

        selectionMode   all;

        U               U;

        velocity        scale;
        value           (-2.572 0 0);
        scale
        {
            type            halfCosineRamp;
            start           0;
            duration        10;
        }
    }
    \endverbatim

SourceFiles
    accelerationSource.C

\*---------------------------------------------------------------------------*/

#ifndef accelerationSource_H
#define accelerationSource_H

#include "fvModel.H"
#include "fvCellSet.H"
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
    public fvModel
{
    // Private Data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Name of the velocity field
        word UName_;

        //- Time-varying velocity
        autoPtr<Function1<vector>> velocity_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Source term to momentum equation
        template<class AlphaRhoFieldType>
        void add
        (
            const AlphaRhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;


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

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

            //- Source term to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


        // Mesh motion

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);


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
