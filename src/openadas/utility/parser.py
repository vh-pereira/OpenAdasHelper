import numpy as np
import re

from src.utility.conversion import PerCm3ToPerM3, Cm3ToM3


def _group_by_block(source_file, match_string):
    """
    Generator the splits the ADF15 file into blocks.

    Groups lines of file into blocks based on precursor '  6561.9A   24...'

    Note: comment section not filtered out of last block, don't over-read!
    """

    buffer = []
    for line in source_file:
        if re.match(match_string, line, re.IGNORECASE):
            if buffer:
                yield buffer
            buffer = [line]
        else:
            buffer.append(line)
    yield buffer

def extract_rate(file, block_num):
    """
    Reads and converts the rate data for the specified block.
    """

    # search from start of file
    file.seek(0)

    wavelength_match = r"^\s*[0-9]*\.[0-9]* ?a? +.*$"
    block_id_match = r"^\s*[0-9]*\.[0-9]* ?a?\s*([0-9]*)\s*([0-9]*).*/type *= *([a-zA-Z]*).*/isel *= * ([0-9]*)$"

    for block in _group_by_block(file, wavelength_match):
        match = re.match(block_id_match, block[0], re.IGNORECASE)

        if not match:
            continue

        if int(match.groups()[3]) == block_num:
            # get number of n, T and rate data points:
            num_n = int(match.groups()[0])
            num_t = int(match.groups()[1])
            num_r = num_n * num_t

            block.pop(0)

            # Load density values
            nn = 0
            density = []
            while nn != num_n:
                next_line = block.pop(0)
                components = next_line.split()
                for value in components:
                    nn += 1
                    density.append(float(value))

            # Load temperature values
            nt = 0
            temperature = []
            while nt != num_t:
                next_line = block.pop(0)
                components = next_line.split()
                for value in components:
                    nt += 1
                    temperature.append(float(value))

            # Load rate values
            nr = 0
            rates = []
            while nr != num_r:
                next_line = block.pop(0)
                components = next_line.split()
                for value in components:
                    nr += 1
                    rates.append(float(value))

            density = np.array(density)
            temperature = np.array(temperature)
            rates = np.array(rates)
            rates = rates.reshape((num_n, num_t))

            # convert units from cm^-3 to m^-3
            density = PerCm3ToPerM3.to(density)
            rates = Cm3ToM3.to(rates)

            return {'ne': density, 'te': temperature, 'rate': rates}

    # If code gets to here, block wasn't found.
    raise RuntimeError('Block number {} was not found in the ADF15 file.'.format(block_num))