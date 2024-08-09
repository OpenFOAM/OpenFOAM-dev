/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::porosityModel

Description
    Top level model for porosity models

SourceFiles
    porosityModel.C
    porosityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef porosityModel_H
#define porosityModel_H

#include "fvMesh.H"
#include "fvMatricesFwd.H"
#include "coordinateSystem.H"
#include "dimensionedVector.H"
#include "wordRe.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class porosityModel Declaration
\*---------------------------------------------------------------------------*/

class porosityModel
:
    public regIOobject
{
    // Private Member Functions

        //- Disallow default bitwise copy construction
        porosityModel(const porosityModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const porosityModel&) = delete;


protected:

    // Protected data

        //- Porosity name
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Name of cellZone
        word zoneName_;

        //- Local co-ordinate system
        coordinateSystem coordSys_;


    // Protected Member Functions

        //- Transform the model data wrt mesh changes
        virtual void calcTransformModelData() = 0;

        //- Adjust negative resistance values to be multiplier of max value
        void adjustNegativeResistance(dimensionedVector& resist);

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const = 0;

        virtual void correct(fvVectorMatrix& UEqn) const = 0;

        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const = 0;

        //- Return label index
        label fieldIndex(const label index) const;


public:

    //- Runtime type information
    TypeName("porosityModel");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        porosityModel,
        mesh,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& coeffDict,
            const word& cellZoneName
        ),
        (name, mesh, dict, coeffDict, cellZoneName)
    );

    //- Constructor
    porosityModel
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict,
        const dictionary& coeffDict,
        const word& cellZoneName = word::null
    );


    //- Selector
    static autoPtr<porosityModel> New
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Destructor
    virtual ~porosityModel();


    // Member Functions

        //- Return const access to the porosity model name
        inline const word& name() const;

        //- Return const access to the cell zone name
        inline const word& zoneName() const;

        //- Return the force over the cell zone(s)
        virtual tmp<vectorField> force
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add resistance
        virtual void addResistance(fvVectorMatrix& UEqn);

        //- Add resistance
        virtual void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC
        );


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


    // I-O

        //- Write
        virtual bool writeData(Ostream& os) const
        {
            return true;
        }

        //- Inherit read from regIOobject
        using regIOobject::read;

        //- Read porosity dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "porosityModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
