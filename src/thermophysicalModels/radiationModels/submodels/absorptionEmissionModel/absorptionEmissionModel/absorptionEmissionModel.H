/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::radiation::absorptionEmissionModel

Description
    Model to supply absorption and emission coefficients for radiation
    modelling

\*---------------------------------------------------------------------------*/

#ifndef radiationAbsorptionEmissionModel_H
#define radiationAbsorptionEmissionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                   Class absorptionEmissionModel Declaration
\*---------------------------------------------------------------------------*/

class absorptionEmissionModel
{

protected:

    // Protected data

        //- Radiation model dictionary
        const dictionary dict_;

        //- Reference to the fvMesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("absorptionEmissionModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            absorptionEmissionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        absorptionEmissionModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Selector
    static autoPtr<absorptionEmissionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~absorptionEmissionModel();


    // Member Functions

        // Access

            //- Reference to the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }


            // Absorption coefficient

                //- Absorption coefficient (net)
                virtual tmp<volScalarField> a(const label bandI = 0) const;

                //- Absorption coefficient for continuous phase
                virtual tmp<volScalarField> aCont(const label bandI = 0) const;

                //- Absorption coefficient for dispersed phase
                virtual tmp<volScalarField> aDisp(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient (net)
                virtual tmp<volScalarField> e(const label bandI = 0) const;

                //- Return emission coefficient for continuous phase
                virtual tmp<volScalarField> eCont(const label bandI = 0) const;

                //- Return emission coefficient for dispersed phase
                virtual tmp<volScalarField> eDisp(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution (net)
                virtual tmp<volScalarField> E(const label bandI = 0) const;

                //- Emission contribution for continuous phase
                virtual tmp<volScalarField> ECont(const label bandI = 0) const;

                //- Emission contribution for dispersed phase
                virtual tmp<volScalarField> EDisp(const label bandI = 0) const;


            //- Const access to the number of bands - defaults to 1 for grey
            //  absorption/emission
            virtual label nBands() const;

            //- Const access to the bands - defaults to Vector2D::one for grey
            //  absorption/emission
            virtual const Vector2D<scalar>& bands(const label n) const;

            //- Flag for whether the absorption/emission is for a grey gas
            virtual bool isGrey() const;

            //- Correct absorption coefficients
            virtual void correct
            (
                volScalarField& a,
                PtrList<volScalarField>& aj
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

