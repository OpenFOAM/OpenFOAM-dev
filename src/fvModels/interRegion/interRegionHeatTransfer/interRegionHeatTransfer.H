/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::fv::interRegionHeatTransfer

Description
    Model for inter-region heat exchange. Requires specification of a model for
    the heat transfer coefficient (htc) and the area per unit volume (Av).
    These are then used to apply the following source to the energy equation:

    \f[
        -htc*Av*(T_{nbr,mapped} - T)
    \f]

    If the semiImplicit option is set, then this becomes:

    \f[
        -htc*Av*(T_{nbr,mapped} - T) + htc*Av/Cp*h - Sp(htc*Av/Cp, h);
    \f]

Usage
    Example usage:
    \verbatim
    interRegionHeatTransfer
    {
        type            interRegionHeatTransfer;

        libs            ("libinterRegionFvModels.so");

        nbrRegion       other;

        interpolationMethod intersection;
        master          true;

        semiImplicit    no;

        Av              200;

        heatTransferCoefficientModel constant;

        htc             10;
    }
    \endverbatim

See also
    Foam::fv::heatTransferCoefficientModel

SourceFiles
    interRegionHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef interRegionHeatTransfer_H
#define interRegionHeatTransfer_H

#include "interRegionModel.H"
#include "heatTransferCoefficientModel.H"
#include "heatTransferAv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class interRegionHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class interRegionHeatTransfer
:
    public interRegionModel
{
    // Private data

        //- Flag to activate semi-implicit coupling
        bool semiImplicit_;

        //- Name of temperature field; default = "T"
        word TName_;

        //- Name of neighbour temperature field; default = "T"
        word TNbrName_;

        //- The heat transfer area per unit volume
        autoPtr<heatTransferAv> heatTransferAv_;

        //- The heat transfer model
        autoPtr<heatTransferCoefficientModel> heatTransferCoefficientModel_;


    // Private member functions

        //- Non-virtual readalpha*
        void readCoeffs(const dictionary& dict);

        //- Get the neighbour heat transfer
        const interRegionHeatTransfer& nbrHeatTransfer() const;


public:

    //- Runtime type information
    TypeName("interRegionHeatTransfer");


    // Constructors

        //- Construct from dictionary
        interRegionHeatTransfer
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~interRegionHeatTransfer();


    // Member Functions

        // Access

            //- Return the heat transfer coefficient
            inline const volScalarField& htc() const;


        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Source term to energy equation
            virtual void addSup
            (
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Source term to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;


        // Correction

            //- Correct the model
            virtual void correct();


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
