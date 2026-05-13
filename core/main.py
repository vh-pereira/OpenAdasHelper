# Copyright 2016-2018 Euratom
# Copyright 2016-2018 United Kingdom Atomic Energy Authority
# Copyright 2016-2018 Centro de Investigaciones Energéticas, Medioambientales y Tecnológicas
#
# Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the
# European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/software/page/eupl5
#
# Unless required by applicable law or agreed to in writing, software distributed
# under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied.
#
# See the Licence for the specific language governing permissions and limitations
# under the Licence.

from .openadas.utility.utils import DEFAULT_REPOSITORY_PATH
from .openadas import utility as repository
from .utility.rates import *

class OpenADAS:
    """
    OpenADAS atomic data source.

    :param str data_path: OpenADAS local repository path.
    :param bool permit_extrapolation: If true, informs interpolation objects to allow extrapolation
                                      beyond the limits of the tabulated data. Default is False.
    :param bool missing_rates_return_null: If true, allows Null rate objects to be returned when
                                           the requested atomic data is missing. Default is False.
    :param bool wavelength_element_fallback: If true, allows to use the element's wavelength when
                                             the isotope's wavelength is not available.
                                             Default is False.
    """

    def __init__(self, populate=True, data_path=None, permit_extrapolation=False, missing_rates_return_null=False,
                 wavelength_element_fallback=False):

        super().__init__()
        self._data_path = data_path or DEFAULT_REPOSITORY_PATH

        self._permit_extrapolation = permit_extrapolation

        self._missing_rates_return_null = missing_rates_return_null

        self._wavelength_element_fallback = wavelength_element_fallback

        self._populate() if populate else None

    @property
    def data_path(self):
        return self._data_path

    def _populate(self, download=True, repository_path=None, adas_path=None):
        """
        Populates the OpenADAS repository with a typical set of rates and wavelengths.

        If an ADAS file is not note found an attempt will be made to download the
        file from the OpenADAS website. This behaviour can be disabled by setting
        the download argument to False.

        :param download: Attempt to download the ADAS files if missing (default=True).
        :param repository_path: Alternate path for the OpenADAS repository (default=None).
        :param adas_path: Alternate path in which to search for ADAS files (default=None) .
        """

        # install a common selection of open adas files
        rates = {
            'adf15': (
                ('H', 0, 'adf15/pec12#h/pec12#h_pju#h0.dat'),
                ('He', 0, 'adf15/pec96#he/pec96#he_pju#he0.dat'),
                ('He', 1, 'adf15/pec96#he/pec96#he_pju#he1.dat'),
                ('Be', 0, 'adf15/pec96#be/pec96#be_pju#be0.dat'),
                ('Be', 1, 'adf15/pec96#be/pec96#be_pju#be1.dat'),
                ('Be', 2, 'adf15/pec96#be/pec96#be_pju#be2.dat'),
                ('Be', 3, 'adf15/pec96#be/pec96#be_pju#be3.dat'),
                ('C', 0, 'adf15/pec96#c/pec96#c_vsu#c0.dat'),
                ('C', 1, 'adf15/pec96#c/pec96#c_vsu#c1.dat'),
                ('C', 2, 'adf15/pec96#c/pec96#c_vsu#c2.dat'),
                ('C', 5, 'adf15/pec96#c/pec96#c_pju#c5.dat'),
                # (neon,      0, 'adf15/pec96#ne/pec96#ne_pju#ne0.dat'),     #TODO: OPENADAS DATA CORRUPT
                # (neon,      1, 'adf15/pec96#ne/pec96#ne_pju#ne1.dat'),     #TODO: OPENADAS DATA CORRUPT
                ('N', 0, 'adf15/pec96#n/pec96#n_vsu#n0.dat'),
                ('N', 1, 'adf15/pec96#n/pec96#n_vsu#n1.dat'),
                # ('N',  2, 'adf15/pec96#n/pec96#n_vsu#n2.dat'),    #TODO: OPENADAS DATA CORRUPT
            )
        }
        repository.install_files(rates, download=download, repository_path=repository_path, adas_path=adas_path)

    def wavelength(self, ion, charge, transition):
        """
        Spectral line wavelength for a given transition.

        :param ion: Element object defining the ion type.
        :param charge: Charge state of the ion.
        :param transition: Tuple containing (initial level, final level)
        :return: Wavelength in nanometers.
        """

        return repository.get_wavelength(ion, charge, transition, repository_path=self._data_path)

    def impact_excitation_pec(self, ion, charge, transition):
        """
        Electron impact excitation photon emission coefficient for a given species.

        Open-ADAS data is interpolated with cubic spline in log-log space.
        Nearest neighbour extrapolation is used when permit_extrapolation is True.

        :param ion: Element object defining the ion type.
        :param charge: Charge state of the ion.
        :param transition: Tuple containing (initial level, final level).
        :return: Impact excitation photon emission coefficient in ph.m^3.s^-1 as a
                 function of electron density and temperature.
        """

        # extract element from isotope because there are no isotope rates in ADAS
        # keep the isotope for the wavelength
        ion_element = ion

        try:
            # read electron impact excitation PEC from json file in the repository
            data = repository.get_pec_excitation_rate(ion_element, charge, transition, repository_path=self._data_path)

        except RuntimeError:
            if self._missing_rates_return_null:
                return NullImpactExcitationPEC()
            raise

        return ImpactExcitationPEC(data, extrapolate=self._permit_extrapolation)

    def recombination_pec(self, ion, charge, transition):
        """
        Recombination photon emission coefficient for a given species.

        Open-ADAS data is interpolated with cubic spline in log-log space.
        Nearest neighbour extrapolation is used when permit_extrapolation is True.

        :param ion: Element object defining the ion type.
        :param charge: Charge state of the ion after recombination.
        :param transition: Tuple containing (initial level, final level).
        :return: Recombination photon emission coefficient in ph.m^3.s^-1 as a function of electron
                 density and temperature.
        """

        try:
            # read free electron recombination PEC from json file in the repository
            data = repository.get_pec_recombination_rate(ion, charge, transition,
                                                         repository_path=self._data_path)


        except (FileNotFoundError, KeyError):
            if self._missing_rates_return_null:
                return NullRecombinationPEC()
            raise

        # obtain isotope's rest wavelength for a given transition
        # the wavelength is used ot convert the PEC from photons/s/m3 to W/m3
        wavelength = self.wavelength(ion, charge, transition)

        return RecombinationPEC(data, extrapolate=self._permit_extrapolation)

    def thermal_cx_pec(self, donor_element, donor_charge, receiver_element, receiver_charge, transition):
        """
        Thermal CX photon emission coefficient for a given species.

        Open-ADAS data is interpolated with cubic spline in log-log space.
        Nearest neighbour extrapolation is used when permit_extrapolation is True.

        :param donor_element: Element object defining the donor ion type.
        :param donor_charge: Charge state of the donor ion.
        :param receiver_element: Element object defining the receiver ion type.
        :param receiver_charge: Charge state of the receiver ion.
        :param transition: Tuple containing (initial level, final level) of the receiver
                           in charge state receiver_charge - 1.
        :return: Thermal charge exchange photon emission coefficient in ph.m^3.s^-1
                 as a function of electron density, electron temperature and donor temperature.
        """

        try:
            # read thermal CX rate from json file in the repository
            data = repository.get_pec_thermal_cx_rate(donor_element, donor_charge,
                                                      receiver_element, receiver_charge,
                                                      transition,
                                                      repository_path=self._data_path)

        except RuntimeError:
            if self._missing_rates_return_null:
                return NullThermalCXPEC()
            raise

        # obtain isotope's rest wavelength for a given transition
        # the wavelength is used ot convert the PEC from photons/s/m3 to W/m3
        wavelength = self.wavelength(receiver_element, receiver_charge - 1, transition)

        return ThermalCXPEC(data, extrapolate=self._permit_extrapolation)
