/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::fv::interRegionExplicitPorositySource

Description
    Inter-region explicit porosity source.

Usage
    Example usage, here employing the Darcy-Forchheimer model:
    \verbatim
    interRegionExplicitPorositySourceCoeffs
    {
        type            DarcyForchheimer;

        DarcyForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (5e7 -1000 -1000);
            f   f [0 -1 0 0 0 0 0] (0 0 0);

            coordinateSystem
            {
                type    cartesian;
                origin  (0 0 0);
                coordinateRotation
                {
                    type    axesRotation;
                    e1  (0.70710678 0.70710678 0);
                    e2  (0 0 1);
                }
            }
        }
    }
    \endverbatim

SourceFiles
    interRegionExplicitPorositySource.C

\*---------------------------------------------------------------------------*/

#ifndef interRegionExplicitPorositySource_H
#define interRegionExplicitPorositySource_H

#include "interRegionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
              Class interRegionExplicitPorositySource Declaration
\*---------------------------------------------------------------------------*/

class interRegionExplicitPorositySource
:
    public interRegionModel
{
    // Private data

        //- Velocity field name, default = U
        word UName_;

        //- Porosity mapping filter
        volScalarField::Internal filter_;

        //- Run-time selectable porosity model
        mutable autoPtr<porosityModel> porosityPtr_;


    // Protected Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("interRegionExplicitPorositySource");


    // Constructors

        //- Construct from components
        interRegionExplicitPorositySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        interRegionExplicitPorositySource
        (
            const interRegionExplicitPorositySource&
        );


    //- Destructor
    virtual ~interRegionExplicitPorositySource()
    {}


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Add explicit and implicit contributions to compressible equation

            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interRegionExplicitPorositySource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
