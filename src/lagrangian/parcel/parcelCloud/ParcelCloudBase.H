/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::ParcelCloudBase

Description
    Base template for parcel clouds. Inserts the parcelCloudBase virtualisation
    layer into the class. Also defines default zero-return source methods to
    enable the less functional clouds to be used in more complex situations.

\*---------------------------------------------------------------------------*/

#ifndef ParcelCloudBase_H
#define ParcelCloudBase_H

#include "Cloud.H"
#include "parcelCloudBase.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidThermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ParticleType>
class ParcelCloudBase
:
    public Cloud<ParticleType>,
    virtual public parcelCloudBase
{
protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;


public:

    // Constructors

        //- Construct given carrier fields
        ParcelCloudBase
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            const bool readFields = true
        )
        :
            Cloud<ParticleType>(rho.mesh(), cloudName, false),
            mesh_(rho.mesh())
        {}

        //- Construct given carrier fields and thermo
        ParcelCloudBase
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const fluidThermo& carrierThermo,
            const bool readFields = true
        )
        :
            Cloud<ParticleType>(rho.mesh(), cloudName, false),
            mesh_(rho.mesh())
        {}

        //- Copy constructor with new name
        ParcelCloudBase
        (
            ParcelCloudBase<ParticleType>& c,
            const word& name
        )
        :
            Cloud<ParticleType>(c.mesh_, name, c),
            mesh_(c.mesh_)
        {}

        //- Copy constructor with new name - creates bare cloud
        ParcelCloudBase
        (
            const fvMesh& mesh,
            const word& name,
            const ParcelCloudBase<ParticleType>& c
        )
        :
            Cloud<ParticleType>(mesh, name, IDLList<ParticleType>()),
            mesh_(mesh)
        {}


    // Member Functions

        // Access

            //- Return reference to the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }


        // Sources

            // Momentum

                //- Return momentum source term [kg m/s^2]
                tmp<fvVectorMatrix> SU(const volVectorField& U) const
                {
                    return tmp<fvVectorMatrix>
                    (
                        new fvVectorMatrix(U, dimMass*dimAcceleration)
                    );
                }

                //- Momentum transfer [kg m/s]
                tmp<volVectorField::Internal> UTrans() const
                {
                    return volVectorField::Internal::New
                    (
                        this->name() + ":UTrans",
                        this->mesh(),
                        dimensionedVector(dimMass*dimVelocity, Zero)
                    );
                }

                //- Momentum transfer coefficient [kg]
                tmp<volScalarField::Internal> UCoeff() const
                {
                    return volScalarField::Internal::New
                    (
                        this->name() + ":UCoeffs",
                        this->mesh(),
                        dimensionedScalar(dimMass, Zero)
                    );
                }


            // Energy

                //- Return sensible enthalpy source term [J/s]
                tmp<fvScalarMatrix> Sh(const volScalarField& hs) const
                {
                    return tmp<fvScalarMatrix>
                    (
                        new fvScalarMatrix(hs, dimEnergy/dimTime)
                    );
                }

                //- Sensible enthalpy transfer [J]
                tmp<volScalarField::Internal> hsTrans() const
                {
                    return volScalarField::Internal::New
                    (
                        this->name() + ":hsTrans",
                        this->mesh(),
                        dimensionedScalar(dimEnergy, Zero)
                    );
                }

                //- Sensible enthalpy transfer coefficient [J/T]
                tmp<volScalarField::Internal> hsCoeff() const
                {
                    return volScalarField::Internal::New
                    (
                        this->name() + ":hsCoeffs",
                        this->mesh(),
                        dimensionedScalar(dimEnergy/dimTemperature, Zero)
                    );
                }

                //- Return equivalent particulate emission [kg/m/s^3]
                tmp<volScalarField> Ep() const
                {
                    return volScalarField::New
                    (
                        this->name() + ":radiation:Ep",
                        this->mesh(),
                        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
                    );
                }

                //- Return equivalent particulate absorption [1/m]
                tmp<volScalarField> ap() const
                {
                    return volScalarField::New
                    (
                        this->name() + ":radiation:ap",
                        this->mesh(),
                        dimensionedScalar(dimless/dimLength, 0)
                    );
                }

                //- Return equivalent particulate scattering factor [1/m]
                tmp<volScalarField> sigmap() const
                {
                    return volScalarField::New
                    (
                        this->name() + ":radiation:sigmap",
                        this->mesh(),
                        dimensionedScalar(dimless/dimLength, 0)
                    );
                }


            // Mass

                //- Return mass source term for specie [kg/s]
                tmp<fvScalarMatrix> SYi
                (
                    const label i,
                    const volScalarField& Yi
                ) const
                {
                    return tmp<fvScalarMatrix>
                    (
                        new fvScalarMatrix(Yi, dimMass/dimTime)
                    );
                }

                //- Return total mass source term [kg/s]
                tmp<fvScalarMatrix> Srho(const volScalarField& rho) const
                {
                    return tmp<fvScalarMatrix>
                    (
                        new fvScalarMatrix(rho, dimMass/dimTime)
                    );
                }

                //- Return total mass source [kg/m^3/s]
                tmp<volScalarField::Internal> Srho() const
                {
                    return tmp<volScalarField::Internal>
                    (
                        volScalarField::Internal::New
                        (
                            this->name() + ":Srho",
                            this->mesh(),
                            dimensionedScalar(dimDensity/dimTime, Zero)
                        )
                    );
                }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
