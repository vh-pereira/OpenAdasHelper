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

import urllib.parse
import urllib.request
import numpy as np
import os
import json

from . import pec
from ...utility.recursivedict import RecursiveDict
from ...utility.conversion import PerCm3ToPerM3, Cm3ToM3
from .. import parse_files as parser

OPENADAS_FILE_URL = 'http://open.adas.ac.uk/download/'

DEFAULT_REPOSITORY_PATH = os.path.join(os.getcwd(), ".openadas/data")

"""
Utilities for managing the local rate repository.
"""


def encode_transition(transition):
    """
    Generate a key string from a transition.

    Both integer and string transition descriptions are handled.
    """

    upper, lower = transition

    upper = str(upper).lower()
    lower = str(lower).lower()

    return '{} -> {}'.format(upper, lower)

def update_wavelengths(wavelengths, repository_path=None):
    """
    Updates the wavelength files `/wavelength/<species>/<charge>.json`
    in atomic data repository.

    File contains multiple rates, indexed by the transitions.

    :param wavelengths: Dictionary in the form:

    |  { <species>: { <charge>: { <transition>: <wavelength> } } }, where
    |      <species> is the plasma species (Element/Isotope),
    |      <charge> is the charge of the plasma species,
    |      <transition> is the tuple containing (initial level, final level),
    |      <wavelength> is the transition's wavelength in nm.

    :param repository_path: Path to the atomic data repository.
    """

    repository_path = repository_path or DEFAULT_REPOSITORY_PATH

    for element, charge_states in wavelengths.items():
        for charge, transitions in charge_states.items():

            # sanitise and validate

            path = os.path.join(repository_path, 'wavelength/{}/{}.json'.format(element, charge))

            # read in any existing wavelengths
            try:
                with open(path, 'r') as f:
                    content = RecursiveDict.from_dict(json.load(f))
            except FileNotFoundError:
                content = RecursiveDict()

            # add/replace data for a transition
            for transition in transitions:
                key = encode_transition(transition)
                content[key] = float(wavelengths[element][charge][transition])

            # create directory structure if missing
            directory = os.path.dirname(path)
            if not os.path.isdir(directory):
                os.makedirs(directory)

            # write new data
            with open(path, 'w') as f:
                json.dump(content, f, indent=2, sort_keys=True)

def get_wavelength(element, charge, transition, repository_path=None):
    """
    Reads the wavelength for the given species, charge and transition from the repository.

    :param element: Plasma species (Element/Isotope).
    :param charge: Charge of the plasma species.
    :param transition: Tuple containing (initial level, final level).
    :param repository_path: Path to the atomic data repository.

    :return wavelength: Wavelength in nm.
    """

    repository_path = repository_path or DEFAULT_REPOSITORY_PATH
    path = os.path.join(repository_path, 'wavelength/{}/{}.json'.format(element, charge))
    try:
        with open(path, 'r') as f:
            content = json.load(f)
        return content[encode_transition(transition)]
    except (FileNotFoundError, KeyError):
        raise RuntimeError('Requested wavelength (element={}, charge={}, transition={})'
                           ' is not available.'.format(element, charge, transition))

def install_files(configuration, download=False, repository_path=None, adas_path=None):

    for adf in configuration:
        if adf.lower() == 'adf11scd':
            for args in configuration[adf]:
                install_adf11scd(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf11acd':
            for args in configuration[adf]:
                install_adf11acd(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf11ccd':
            for args in configuration[adf]:
                install_adf11ccd(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf11plt':
            for args in configuration[adf]:
                install_adf11plt(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf11prb':
            for args in configuration[adf]:
                install_adf11prb(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf11prc':
            for args in configuration[adf]:
                install_adf11prc(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf12':
            for args in configuration[adf]:
                install_adf12(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf15':
            for args in configuration[adf]:
                install_adf15(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf21':
            for args in configuration[adf]:
                install_adf21(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf22bmp':
            for args in configuration[adf]:
                install_adf22bmp(*args, download=download, repository_path=repository_path, adas_path=adas_path)
        if adf.lower() == 'adf22bme':
            for args in configuration[adf]:
                install_adf22bme(*args, download=download, repository_path=repository_path, adas_path=adas_path)


# todo: move print calls to logging

def install_adf11scd(element, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the ionisation rate defined in an ADF11 file to the repository.

    :param element: The element described by the rate file.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # decode file and write out rates
    rate_adas = parse_adf11(element, path)

    rate_cherab = _notation_adf11_adas2cherab(rate_adas, "scd")  # convert from adas to cherab notation

    repository.update_ionisation_rates(rate_cherab, repository_path)

def install_adf11acd(element, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the recombination rate defined in an ADF11 file to the repository.

    :param element: The element described by the rate file.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # decode file and write out rates
    rate_adas = parse_adf11(element, path)

    rate_cherab = _notation_adf11_adas2cherab(rate_adas, "acd")  # convert from adas to cherab notation

    repository.update_recombination_rates(rate_cherab, repository_path)

def install_adf11ccd(donor_element, donor_charge, receiver_element, file_path, download=False,
                     repository_path=None, adas_path=None):
    """
    Adds the thermal charge exchange rate defined in an ADF11 file to the repository.

    :param donor_element: Element donating the electron, for the case of ADF11 files it is
      neutral hydrogen.
    :param donor_charge: Charge of the donor atom/ion.
    :param receiver_element: Element receiving the electron.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # decode file and write out rates
    rate_adas = parse_adf11(receiver_element, path)
    rate_cherab = _notation_adf11_adas2cherab(rate_adas, "ccd")  # convert from adas to cherab notation

    # reshape rate dictionary to match cherab convention
    rate_cherab_ccd = RecursiveDict()
    rate_cherab_ccd[donor_element][donor_charge] = rate_cherab

    repository.update_thermal_cx_rates(rate_cherab_ccd, repository_path)

def install_adf11plt(element, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the line radiated power rates defined in an ADF11 file to the repository.

    :param element: The element described by the rate file.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # decode file and write out rates
    rate_adas = parse_adf11(element, path)

    rate_cherab = _notation_adf11_adas2cherab(rate_adas, "plt")  # convert from adas to cherab notation

    repository.update_line_power_rates(rate_cherab, repository_path)

def install_adf11prb(element, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the continuum radiated power rates defined in an ADF11 file to the repository.

    :param element: The element described by the rate file.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # decode file and write out rates
    rate_adas = parse_adf11(element, path)

    rate_cherab = _notation_adf11_adas2cherab(rate_adas, "prb")  # convert from adas to cherab notation

    repository.update_continuum_power_rates(rate_cherab, repository_path)

def install_adf11prc(element, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the CX radiated power rates defined in an ADF11 file to the repository.

    :param element: The element described by the rate file.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # decode file and write out rates
    rate_adas = parse_adf11(element, path)

    rate_cherab = _notation_adf11_adas2cherab(rate_adas, "prc")  # convert from adas to cherab notation

    repository.update_cx_power_rates(rate_cherab, repository_path)

def install_adf12(donor_ion, donor_metastable, receiver_ion, receiver_charge, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the rates in the ADF12 file to the repository.

    :param donor_ion: The donor ion element described by the rate file.
    :param donor_metastable: The donor ion metastable level.
    :param receiver_ion: The receiver ion element described by the rate file.
    :param receiver_charge: The receiver ion ionisation level described by the rate file.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # decode file and write out rates
    rates = parse_adf12(donor_ion, donor_metastable, receiver_ion, receiver_charge, path)
    repository.update_beam_cx_rates(rates, repository_path)

def install_adf15(element, ionisation, file_path, download=False, repository_path=None, adas_path=None, header_format=None):
    """
    Adds the rates in the ADF15 file to the repository.

    :param element: The element described by the rate file.
    :param ionisation: The ionisation level described by the rate file.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # decode file and write out rates
    rates, wavelengths = parser.parse_adf15(element, ionisation, path, header_format=header_format)

    if (rates, wavelengths) == (None, None):
        return None

    if 'thermalcx' in rates:
        cx_rates = rates.pop('thermalcx')
        # CX rates for Tdon = Trec (2D function of Ne, Te)
        # converting to 3D function of Ne, Te, Tdon
        pec.update_pec_thermal_cx_rates(_thermalcx_adf15_2dto3d_converter(cx_rates))

    pec.update_pec_rates(rates, repository_path)
    update_wavelengths(wavelengths, repository_path)

def install_adf21(beam_species, target_ion, target_charge, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the beam stopping rate defined in an ADF21 file to the repository.

    :param beam_species: Beam neutral atom (Element/Isotope).
    :param target_ion: Target species (Element/Isotope).
    :param target_charge: Charge of the target species.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # # decode file and write out rates
    rate = parse_adf21(beam_species, target_ion, target_charge, path)
    repository.update_beam_stopping_rates(rate, repository_path)

def install_adf22bmp(beam_species, beam_metastable, target_ion, target_charge, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the beam population rate defined in an ADF22 BMP file to the repository.

    :param beam_species: Beam neutral atom (Element/Isotope).
    :param beam_metastable: Metastable/excitation level of beam neutral atom.
    :param target_ion: Target species (Element/Isotope).
    :param target_charge: Charge of the target species.
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # # decode file and write out rates
    rate = parse_adf22bmp(beam_species, beam_metastable, target_ion, target_charge, path)
    repository.update_beam_population_rates(rate, repository_path)

def install_adf22bme(beam_species, target_ion, target_charge, transition, file_path, download=False, repository_path=None, adas_path=None):
    """
    Adds the beam emission rate defined in an ADF22 BME file to the repository.

    :param beam_species: Beam neutral atom (Element/Isotope).
    :param target_ion: Target species (Element/Isotope).
    :param target_charge: Charge of the target species.
    :param transition: Tuple containing (initial level, final level).
    :param file_path: Path relative to ADAS root.
    :param download: Attempt to download file if not present (Default=True).
    :param repository_path: Path to the repository in which to install the rates (optional).
    :param adas_path: Path to ADAS files repository (optional).
    """

    print('Installing {}...'.format(file_path))
    path = _locate_adas_file(file_path, download, adas_path, repository_path)
    if not path:
        raise ValueError('Could not locate the specified ADAS file.')

    # # decode file and write out rates
    rate = parse_adf22bme(beam_species, target_ion, target_charge, transition, path)
    repository.update_beam_emission_rates(rate, repository_path)

def _locate_adas_file(file_path, download=False, adas_path=None, repository_path=None):

    path = None
    repository_path = repository_path or DEFAULT_REPOSITORY_PATH

    # is file in adas path?
    if adas_path:
        target = os.path.join(adas_path, file_path)
        if os.path.isfile(target):
            path = target

    # download file?
    if not path and download:
        target = os.path.join(repository_path, "_download_cache", file_path)

        # is file in cache? if not download...
        if os.path.isfile(target):
            print(' - ADF file already in cache')
            path = target
        else:

            # create directory structure if missing
            directory = os.path.dirname(target)
            if not os.path.isdir(directory):
                os.makedirs(directory)

            print(" - downloading ADF file '{}' to '{}'".format(file_path, target))

            url = urllib.parse.urljoin(OPENADAS_FILE_URL, file_path.replace('#', '][').lstrip('/'))
            urllib.request.urlretrieve(url, target)
            path = target

    return path

def _notation_adf11_adas2cherab(rate_adas, filetype):
    """
    Converts adas unit, charge and numeric notation to cherab notation

    :param rate_adas: Nested dictionary of shape rate_adas[element][charge][te, ne, rates]
    :param filetype: string denoting adas adf11 file type to decide whether charge conversion is to be applied.
      Will be applied for file types: "scd", "ccd", "plt", "pls"
    :return: nested dictionary with cherab rates and units notation
    """

    # Charge correction will be applied if there is difference between adas and cherab charge notation
    if filetype in ["scd", "plt", "pls"]:
        charge_correction = int(-1)
    else:
        charge_correction = int(0)

    # adas units, charge and number notation to be changed to cherab notation
    rate_cherab = RecursiveDict()
    for i in rate_adas.keys():
        for j in rate_adas[i].keys():
            # convert from adas log10 in [cm**-3] notation to cherab [m**-3] electron density notation
            rate_cherab[i][j + charge_correction]["ne"] = PerCm3ToPerM3.to(10**rate_adas[i][j]["ne"])
            # convert from adas log10 to cherab electron temperature notation
            rate_cherab[i][j + charge_correction]["te"] = 10**rate_adas[i][j]["te"]
            rate_cherab[i][j + charge_correction]["rates"] = Cm3ToM3.to(10**rate_adas[i][j]["rates"])

    return rate_cherab

def _thermalcx_adf15_2dto3d_converter(rates):
    """
    Converts thermal CX PEC rates parsed from a standard ADF 15 file
    to the format supported by the repository.

    In the standard ADF 15 file, it is assumed that the donor is H0 and Tdon = Trec.
    """
    new_rates = RecursiveDict()
    for element, charge_states in rates.items():
        for charge, transitions in charge_states.items():
            for transition, rate in transitions.items():
                data = np.empty((len(rate['ne']), len(rate['te']), 2))
                data[:, :, :] = rate['rate'][:, :, None]
                new_rate = {'ne': rate['ne'],
                            'te': rate['te'],
                            'td': np.array([0.01, 10000]),
                            'rate': data}
                new_rates['H'][0][element][charge + 1][transition] = new_rate

    return new_rates