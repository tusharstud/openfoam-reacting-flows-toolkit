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

\*---------------------------------------------------------------------------*/

#include "Kavanau.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Kavanau<CloudType>::Kavanau
(
    const dictionary& dict,
    CloudType& cloud
)
:
    HeatTransferModel<CloudType>(dict, cloud, typeName)
{}


template<class CloudType>
Foam::Kavanau<CloudType>::Kavanau(const Kavanau<CloudType>& htm)
:
    HeatTransferModel<CloudType>(htm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Kavanau<CloudType>::~Kavanau()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::Kavanau<CloudType>::Nu
(
    const scalar Re,
    const scalar Pr
) const
{
    // 1. Standard Ranz-Marshall Nusselt Number
    const scalar NuRanz = 2.0 + 0.6 * Foam::sqrt(max(Re, scalar(0.0))) * Foam::cbrt(max(Pr, scalar(1e-4)));

    // 2. Kavanau transition regime correction: Nu_Kava = Nu_Ranz / (1 + 3.42 * Nu_Ranz * Ma / (Re * Pr))
    scalar corr = 0.0;
    if (this->coeffDict().found("transitionFactor"))
    {
        const scalar tf = this->coeffDict().template lookup<scalar>("transitionFactor");
        corr = 3.42 * NuRanz * tf / max(Pr, scalar(1e-4));
    }
    else if (this->coeffDict().found("Ma"))
    {
        const scalar Ma = this->coeffDict().template lookup<scalar>("Ma");
        corr = 3.42 * NuRanz * Ma / (max(Re, scalar(1e-4)) * max(Pr, scalar(1e-4)));
    }

    return NuRanz / (1.0 + max(corr, scalar(0.0)));
}


template<class CloudType>
Foam::scalar Foam::Kavanau<CloudType>::htc
(
    const scalar dp,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW
) const
{
    const scalar NuRanz = 2.0 + 0.6 * Foam::sqrt(max(Re, scalar(0.0))) * Foam::cbrt(max(Pr, scalar(1e-4)));

    scalar NuKava = NuRanz;
    if (this->coeffDict().found("transitionFactor"))
    {
        const scalar tf = this->coeffDict().template lookup<scalar>("transitionFactor");
        NuKava = NuRanz / (1.0 + 3.42 * NuRanz * tf / max(Pr, scalar(1e-4)));
    }
    else if (this->coeffDict().found("Ma"))
    {
        const scalar Ma = this->coeffDict().template lookup<scalar>("Ma");
        NuKava = NuRanz / (1.0 + 3.42 * NuRanz * Ma / (max(Re, scalar(1e-4)) * max(Pr, scalar(1e-4))));
    }
    else
    {
        // Dynamic formulation: Ma / Re = mu / (rho * c_sound * dp)
        const scalar Tref = max(this->owner().constProps().T0(), scalar(300.0));
        const scalar cSound = 340.0 * Foam::sqrt(Tref / 300.0);
        const scalar rhoGas = 1.18 * (300.0 / Tref);
        const scalar muGas = 1.8e-5 * Foam::pow(Tref / 300.0, 0.7);
        const scalar MaOverRe = muGas / max(rhoGas * cSound * max(dp, scalar(1e-8)), scalar(1e-12));
        NuKava = NuRanz / (1.0 + 3.42 * NuRanz * MaOverRe / max(Pr, scalar(1e-4)));
    }

    scalar htc = NuKava * kappa / max(dp, scalar(1e-9));

    if (this->BirdCorrection() && (mag(htc) > rootVSmall) && (mag(NCpW) > rootVSmall))
    {
        const scalar phit = min(NCpW / htc, scalar(50.0));
        if (phit > 0.001)
        {
            htc *= phit / (exp(phit) - 1.0);
        }
    }

    return htc;
}

// ************************************************************************* //
