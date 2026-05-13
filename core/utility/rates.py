
# Copyright 2016-2024 Euratom
# Copyright 2016-2024 United Kingdom Atomic Energy Authority
# Copyright 2016-2024 Centro de Investigaciones Energéticas, Medioambientales y Tecnológicas
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

import numpy as np
from scipy.interpolate import RegularGridInterpolator



class IonisationPEC:
    """
    Ionisation PEC.

    Data is interpolated with cubic spline in log-log space.
    Nearest neighbour extrapolation is used if extrapolate is True.

    :param dict data: Ionisation PEC dictionary containing the following entries:

    |   'ne': 1D array of size (N) with electron density in m^-3,
    |   'te': 1D array of size (M) with electron temperature in eV,
    |   'PEC': 2D array of size (N, M) with ionisation PEC in m^3.s^-1.

    :param bint extrapolate: Enable extrapolation (default=False).

    :ivar tuple density_range: Electron density interpolation range.
    :ivar tuple temperature_range: Electron temperature interpolation range.
    :ivar dict raw_data: Dictionary containing the raw data.
    """

    def __init__(self, data:dict, extrapolate=False):

        self.raw_data = data

        # unpack
        ne = data['ne']
        te = data['te']
        pec = np.log10(data['rate'])

        # store limits of data
        self.density_range = ne.min(), ne.max()
        self.temperature_range = te.min(), te.max()

        # interpolate PEC
        # using nearest extrapolation to avoid infinite values at 0 for some PECs
        extrapolation_type = np.nan if not extrapolate else 'none'
        self._pec = RegularGridInterpolator(
            (np.log10(ne), np.log10(te)),
            pec,
            'cubic',
            bounds_error=False,
            fill_value=extrapolation_type
            )

    def evaluate(self, density, temperature):

        # need to handle zeros, also density and temperature can become negative due to cubic interpolation
        if density.all() <= 0 or temperature.all() <= 0:
            return 0

        pts = np.column_stack((np.log10(density), np.log10(temperature)))
        # calculate PEC and convert from log10 space to linear space
        return 10 ** self._pec(pts)

class NullIonisationPEC:
    """
    An ionisation PEC that always returns zero.
    Needed for use cases where the required atomic data is missing.
    """
    def __init__(self):
        pass

    def evaluate(self, density, temperature):
        return 0.0

class RecombinationPEC:
    """
    Recombination PEC.

    Data is interpolated with cubic spline in log-log space.
    Nearest neighbour extrapolation is used if extrapolate is True.

    :param dict data: Recombination PEC dictionary containing the following entries:

    |       'ne': 1D array of size (N) with electron density in m^-3,
    |       'te': 1D array of size (M) with electron temperature in eV,
    |       'PEC': 2D array of size (N, M) with recombination PEC in m^3.s^-1.

    :param bint extrapolate: Enable extrapolation (default=False).

    :ivar tuple density_range: Electron density interpolation range.
    :ivar tuple temperature_range: Electron temperature interpolation range.
    :ivar dict raw_data: Dictionary containing the raw data.
    """

    def __init__(self, data:dict, extrapolate=False):

        self.raw_data = data

        # unpack
        ne = data['ne']
        te = data['te']
        pec = np.log10(data['rate'])

        # store limits of data
        self.density_range = ne.min(), ne.max()
        self.temperature_range = te.min(), te.max()

        # interpolate PEC
        # using nearest extrapolation to avoid infinite values at 0 for some PECs
        extrapolation_type = np.nan if not extrapolate else 'none'
        self._pec = RegularGridInterpolator(
            (np.log10(ne), np.log10(te)),
            pec,
            'cubic',
            bounds_error=False,
            fill_value=extrapolation_type
        )

    def evaluate(self, density, temperature):
        # need to handle zeros, also density and temperature can become negative due to cubic interpolation
        if density.all() <= 0 or temperature.all() <= 0:
            return 0

        pts = np.column_stack((np.log10(density), np.log10(temperature)))
        # calculate PEC and convert from log10 space to linear space
        return 10 ** self._pec(pts)

class NullRecombinationPEC:
    """
    A recombination PEC that always returns zero.
    Needed for use cases where the required atomic data is missing.
    """
    def __init__(self):
        pass

    def evaluate(self, density, temperature):
        return 0.0

class ThermalCXPEC:
    """
    Thermal charge exchange PEC.

    Data is interpolated with cubic spline in log-log space.
    Linear extrapolation is used if extrapolate is True.

    :param dict data: CX PEC dictionary containing the following entries:

    |       'ne': 1D array of size (N) with electron density in m^-3,
    |       'te': 1D array of size (M) with electron temperature in eV,
    |       'PEC': 2D array of size (N, M) with thermal CX PEC in m^3.s^-1.

    :param bint extrapolate: Enable extrapolation (default=False).

    :ivar tuple density_range: Electron density interpolation range.
    :ivar tuple temperature_range: Electron temperature interpolation range.
    :ivar dict raw_data: Dictionary containing the raw data.
    """

    def __init__(self, data:dict, extrapolate=False):

        self.raw_data = data

        # unpack
        ne = data['ne']
        te = data['te']
        pec = np.log10(data['rate'])

        # store limits of data
        self.density_range = ne.min(), ne.max()
        self.temperature_range = te.min(), te.max()

        # interpolate PEC
        extrapolation_type = np.nan if not extrapolate else 'none'
        self._pec = RegularGridInterpolator(
            (np.log10(ne), np.log10(te)),
            pec,
            'cubic',
            bounds_error=False,
            fill_value=extrapolation_type
        )

    def evaluate(self, density, temperature):

        # need to handle zeros, also density and temperature can become negative due to cubic interpolation
        if density.all() <= 0 or temperature.all() <= 0:
            return 0

        pts = np.column_stack((np.log10(density), np.log10(temperature)))
        # calculate PEC and convert from log10 space to linear space
        return 10 ** self._pec(pts)

class NullThermalCXPEC:
    """
    A thermal CX PEC that always returns zero.
    Needed for use cases where the required atomic data is missing.
    """
    def __init__(self):
        pass

    def evaluate(self, density, temperature):
        return 0.0

class ImpactExcitationPEC:
    """
    Impact excitation PEC.

    Data is interpolated with cubic spline in log-log space.
    Linear extrapolation is used if extrapolate is True.

    :param dict data: Exc PEC dictionary containing the following entries:

    |       'ne': 1D array of size (N) with electron density in m^-3,
    |       'te': 1D array of size (M) with electron temperature in eV,
    |       'PEC': 2D array of size (N, M) with impact exc PEC in m^3.s^-1.

    :param bint extrapolate: Enable extrapolation (default=False).

    :ivar tuple density_range: Electron density interpolation range.
    :ivar tuple temperature_range: Electron temperature interpolation range.
    :ivar dict raw_data: Dictionary containing the raw data.
    """

    def __init__(self, data:dict, extrapolate=False):

        self.raw_data = data

        # unpack
        ne = data['ne']
        te = data['te']
        pec = np.log10(data['rate'])

        # store limits of data
        self.density_range = ne.min(), ne.max()
        self.temperature_range = te.min(), te.max()

        # interpolate PEC
        extrapolation_type = np.nan if not extrapolate else 'none'
        self._pec = RegularGridInterpolator(
            (np.log10(ne), np.log10(te)),
            pec,
            'cubic',
            bounds_error=False,
            fill_value=extrapolation_type
        )

    def evaluate(self, density, temperature):

        # need to handle zeros, also density and temperature can become negative due to cubic interpolation
        if density.all() <= 0 or temperature.all() <= 0:
            return 0

        pts = np.column_stack((np.log10(density), np.log10(temperature)))
        # calculate PEC and convert from log10 space to linear space
        return 10 ** self._pec(pts)

class NullImpactExcitationPEC:
    """
        A impact excitation PEC that always returns zero.
        Needed for use cases where the required atomic data is missing.
        """

    def __init__(self):
        pass

    def evaluate(self, density, temperature):
        return 0.0