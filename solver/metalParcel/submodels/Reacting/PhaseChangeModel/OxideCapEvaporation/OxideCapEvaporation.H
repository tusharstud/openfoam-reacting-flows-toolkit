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
    Foam::OxideCapEvaporation

Description
    Aluminum droplet evaporation model with oxide cap hindrance and Clausius-Clapeyron
    vapor pressure (Zhang et al., 2022, Eqs. 16-20, 26-32).

\*---------------------------------------------------------------------------*/

#ifndef OxideCapEvaporation_H
#define OxideCapEvaporation_H

#include "PhaseChangeModel.H"
#include "liquidMixtureProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class OxideCapEvaporation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class OxideCapEvaporation
:
    public PhaseChangeModel<CloudType>
{
protected:

    // Protected data

        //- Global liquid properties data
        const liquidMixtureProperties& liquids_;

        //- Whether or not to permit condensation (default false)
        Switch condensation_;

        //- List of active liquid names
        List<word> activeLiquids_;

        //- Mapping between liquid and carrier species
        List<label> liqToCarrierMap_;

        //- Mapping between local and global liquid species
        List<label> liqToLiqMap_;

        //- Vapor pressure Clausius-Clapeyron constants: ln(Psat/Pref) = A - B/T
        scalar A_clausius_;
        scalar B_clausius_;
        scalar Pref_;

        //- Alumina melting temperature [K]
        scalar TmeltOxide_;

        //- Latent heat of evaporation [J/kg]
        scalar L_evap_;


    // Protected Member Functions

        //- Sherwood number as a function of Reynolds and Schmidt numbers
        scalar Sh(const scalar Re, const scalar Sc) const;

        //- Calculate the carrier phase component volume fractions at celli
        tmp<scalarField> calcXc(const label celli) const;

        //- Dynamic oxide cap blocking fraction beta (Eqs. 26-32)
        scalar calcBeta(const scalar Tp) const;


public:

    //- Runtime type information
    TypeName("OxideCapEvaporation");


    // Constructors

        //- Construct from dictionary
        OxideCapEvaporation(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        OxideCapEvaporation(const OxideCapEvaporation<CloudType>& pcm);

        //- Construct and return a clone
        virtual autoPtr<PhaseChangeModel<CloudType>> clone() const
        {
            return autoPtr<PhaseChangeModel<CloudType>>
            (
                new OxideCapEvaporation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~OxideCapEvaporation();


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const scalar dt,
            const scalar Re,
            const scalar Pr,
            const scalar d,
            const scalar nu,
            const scalar T,
            const scalar Ts,
            const scalar pc,
            const scalar Tc,
            const scalarField& X,
            scalarField& dMassPC
        ) const;

        //- Return the enthalpy per unit mass
        virtual scalar dh
        (
            const label idc,
            const label idl,
            const scalar p,
            const scalar T
        ) const;

        //- Return vapourisation temperature
        virtual scalar Tvap(const scalarField& X) const;

        //- Return maximum/limiting temperature
        virtual scalar TMax(const scalar p, const scalarField& X) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "OxideCapEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
