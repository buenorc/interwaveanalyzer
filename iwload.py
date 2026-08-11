# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno
Interwave Analyzer - Load variables module

Interwave Analyzer - Version 2 (2026)
Load variables module version: 2.260810

-------------------------------------------------------------------------------

de Carvalho Bueno, R; Bleninger, T. B.; Lorke, A.
Internal wave analyzer for thermally stratified lakes
Environmental Modelling & Software, Elsevier, 2020


Developed by Rafael de Carvalho Bueno
https://buenorc.github.io/

Improvements and betterments by
Andreas Lorke & Tobias Bleninger

Report problems and improvements to email adresss below
decarvalhobueno@gmail.com

for more information, see:
https://buenorc.github.io/pages/interwave.html

-------------------------------------------------------------------------------

PURPOSE
-------

This module reads every input file of Interwave Analyzer.  Compared with
version 1.x it introduces three changes:

1.  *Format detection.*  The column separator, the presence of a header and
    the layout of the date columns are detected from the file itself, so the
    same analysis accepts tab, comma, semicolon or whitespace separated files
    and either the five-column calendar layout or a single ISO 8601 timestamp.

2.  *Vectorised reading.*  Files are parsed with a single NumPy call instead of
    a per-cell Python loop, and the file is no longer read twice (once to count
    the lines, once to read them).

3.  *Physically meaningful time matching.*  Auxiliary records (meteorology,
    water level) are mapped onto the analysis grid by the operators of
    :mod:`iwtime` rather than by nearest-neighbour search on a pseudo-serial
    number.  Wind direction is combined as a vector and the wind speed is also
    returned in its stress-equivalent form.
"""

import os
import re
import io

import numpy as np

import iwpars as pars   # package of additional parameters
import iwmod as mod     # package of functions
import iwtime as tmb    # time base and resampling operators


# Tokens accepted as "no data" in any numeric column
_NA_TOKENS = ('', 'na', 'n/a', 'nan', 'null', 'none', '-', '--', '#n/a',
              '#value!', '#div/0!', '?')

# Candidate column separators, in the order they are tested
_DELIMITERS = ('\t', ';', ',', '|')

_ISO_RE = re.compile(
    r'^\s*\d{4}[-/]\d{1,2}[-/]\d{1,2}([ T]\d{1,2}:\d{2}(:\d{2})?)?\s*$')
_DATE_ONLY_RE = re.compile(r'^\s*\d{4}[-/]\d{1,2}[-/]\d{1,2}\s*$')
_TIME_RE = re.compile(r'^\s*\d{1,2}:\d{2}(:\d{2})?\s*$')


# =============================================================================
#  Low level helpers
# =============================================================================

def _clean_path(path):
    """Strip stray newlines and surrounding quotes from a user supplied path."""
    if path is None:
        return None
    path = str(path).strip().strip('"').strip("'")
    return path.replace('\n', '').replace('\r', '')


def _is_number(token):
    """True when ``token`` parses as a finite or non-finite float."""
    try:
        float(token)
        return True
    except (TypeError, ValueError):
        return False


def _to_float(token):
    """Parse a token to float, mapping recognised no-data markers to NaN."""
    if token is None:
        return np.nan
    text = token.strip()
    if text.lower() in _NA_TOKENS:
        return np.nan
    # Accept the decimal comma used by many European loggers, but only when
    # the comma is not the column separator (handled by the caller).
    try:
        return float(text)
    except ValueError:
        try:
            return float(text.replace(',', '.'))
        except ValueError:
            return np.nan


def _read_sample(path, n_lines=60):
    """Return the first ``n_lines`` non-empty lines of a text file."""
    lines = []
    with open(path, 'r', errors='replace') as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if line.strip():
                lines.append(line)
            if len(lines) >= n_lines:
                break
    return lines


def detect_delimiter(lines):
    """
    Detect the column separator of a tabular text file.

    Inputs:
        lines : list[str] - sample lines taken from the file

    Outputs:
        str or None - the separator; ``None`` means "any run of whitespace",
        which is the convention used by NumPy readers.

    Notes:
        A separator is accepted when it splits every sample line into the same
        number of fields, and that number is larger than one.  Candidates are
        tested in a fixed order so that a tab-separated file containing decimal
        commas is not mistaken for a comma-separated file.
    """
    if not lines:
        return None

    body = lines[1:] if len(lines) > 1 else lines

    for delimiter in _DELIMITERS:
        counts = {len(line.split(delimiter)) for line in body}
        if len(counts) == 1 and counts.pop() > 1:
            return delimiter

    counts = {len(line.split()) for line in body}
    if len(counts) == 1 and counts.pop() > 1:
        return None

    # Inconsistent line lengths: fall back to the separator that appears most
    # often, so that the reader can still emit a useful diagnostic.
    scores = {d: sum(line.count(d) for line in body) for d in _DELIMITERS}
    best = max(scores, key=scores.get)
    return best if scores[best] else None


def _split(line, delimiter):
    """Split a line according to a detected separator."""
    return line.split(delimiter) if delimiter else line.split()


def detect_header(lines, delimiter):
    """
    Decide whether the first line of the file is a header.

    Inputs:
        lines     : list[str] - sample lines
        delimiter : str or None - column separator

    Outputs:
        bool - True when the first line must be skipped.

    Notes:
        The first line is treated as a header when at least one of its fields
        cannot be parsed as a number.  A file whose very first record happens
        to be numeric is therefore read in full, which is the desired behaviour
        for headerless files.
    """
    if not lines:
        return False
    fields = _split(lines[0], delimiter)
    return not all(_is_number(f) for f in fields if f.strip() != '')


def detect_time_layout(rows, header=None):
    """
    Identify how the timestamp is encoded in the first columns of a record.

    Inputs:
        rows   : list[list[str]] - fields of several sample data lines
        header : list[str] or None - header fields, when the file has one

    Outputs:
        (layout, n_time_columns)

            layout in {'ymdhm', 'ymdhms', 'iso', 'date_time'}

            'ymdhm'     : Year, Month, Day, Hour, Minute  (5 columns, legacy)
            'ymdhms'    : Year, Month, Day, Hour, Minute, Second (6 columns)
            'iso'       : one ISO 8601 timestamp, e.g. 2015-06-17T07:30
            'date_time' : an ISO date and a clock time in two columns

    Notes:
        The legacy five-column layout remains fully supported and is still the
        format written by most thermistor loggers; the additional layouts exist
        so that files exported from spreadsheets or databases can be used
        without pre-processing.

        Telling a six-column calendar apart from a five-column calendar
        followed by a data column cannot be done from a single record: a water
        level of exactly 22.0 m looks like a valid seconds field.  Two
        independent criteria are therefore combined.  When the file has a
        header, the name of the sixth column decides.  Otherwise the sixth
        field must be an integer in [0, 60) in *every* sampled record, which no
        physical measurement satisfies in practice.
    """
    if not rows or not rows[0]:
        return 'ymdhm', 5

    fields = rows[0]
    first = fields[0].strip()

    if _DATE_ONLY_RE.match(first) and len(fields) > 1 \
            and _TIME_RE.match(fields[1].strip()):
        return 'date_time', 2

    if _ISO_RE.match(first):
        return 'iso', 1

    if len(fields) < 6:
        return 'ymdhm', 5

    # ---- decisive criterion: the header names the sixth column --------------
    if header and len(header) >= 6:
        name = header[5].strip().strip('"').lower()
        return ('ymdhms', 6) if name.startswith('sec') or name in ('s', 'ss') \
            else ('ymdhm', 5)

    # ---- headerless file: every sampled record must look like a seconds field
    def _is_second(token):
        if not _is_number(token):
            return False
        value = float(token)
        return 0 <= value < 60 and float(value).is_integer()

    def _is_calendar(row):
        try:
            year, month, day = float(row[0]), float(row[1]), float(row[2])
            hour, minute = float(row[3]), float(row[4])
        except (ValueError, IndexError):
            return False
        return (1900 <= year <= 2200 and 1 <= month <= 12 and 1 <= day <= 31
                and 0 <= hour <= 24 and 0 <= minute < 60)

    usable = [r for r in rows if len(r) >= 6]
    if usable and all(_is_calendar(r) and _is_second(r[5]) for r in usable):
        return 'ymdhms', 6

    return 'ymdhm', 5


class TableFormat(object):
    """
    Description of the layout of an input file.

    Attributes:
        path        : str  - file inspected
        delimiter   : str or None - column separator (None = whitespace)
        has_header  : bool - whether the first line must be skipped
        header      : list[str] - header fields (empty when absent)
        n_columns   : int  - number of columns of a data record
        layout      : str  - timestamp layout, see :func:`detect_time_layout`
        n_time      : int  - number of columns occupied by the timestamp
        n_data      : int  - number of remaining (value) columns
    """

    def __init__(self, path, delimiter, has_header, header, n_columns,
                 layout, n_time, timed=True):
        self.path = path
        self.delimiter = delimiter
        self.has_header = has_header
        self.header = header
        self.n_columns = n_columns
        self.layout = layout
        self.n_time = n_time
        self.timed = timed
        self.n_data = max(0, n_columns - n_time)

    def describe(self):
        """One-line human readable summary, written to the diagnosis file."""
        names = {'\t': 'tab', ';': 'semicolon', ',': 'comma', '|': 'pipe'}
        sep = names.get(self.delimiter, 'whitespace')
        text = (f"{os.path.basename(self.path)}: {sep}-separated, "
                f"{'with' if self.has_header else 'no'} header, "
                f"{self.n_columns} columns")
        if self.timed:
            text += (f", time layout '{self.layout}', "
                     f"{self.n_data} data column(s)")
        return text


def inspect(path, timed=True):
    """
    Detect the format of a tabular input file.

    Inputs:
        path  : str  - file to inspect
        timed : bool - False for files without a timestamp (sensor and
                transect files), so that no time layout is inferred

    Outputs:
        TableFormat

    Raises:
        IOError - the file does not exist, is empty or holds no data record.
    """
    path = _clean_path(path)
    if not path or not os.path.isfile(path):
        raise IOError(f"input file not found: {path}")

    lines = _read_sample(path)
    if not lines:
        raise IOError(f"input file is empty: {path}")

    delimiter = detect_delimiter(lines)
    has_header = detect_header(lines, delimiter)

    header = [h.strip() for h in _split(lines[0], delimiter)] \
        if has_header else []
    body = lines[1:] if has_header else lines
    if not body:
        raise IOError(f"input file contains no data records: {path}")

    rows = [_split(line, delimiter) for line in body]

    if timed:
        layout, n_time = detect_time_layout(rows, header)
    else:
        layout, n_time = 'none', 0

    return TableFormat(path, delimiter, has_header, header, len(rows[0]),
                       layout, n_time, timed=timed)


# =============================================================================
#  Generic timed table reader
# =============================================================================

def read_timed_table(path, fmt=None, max_columns=None):
    """
    Read a tabular file whose leading columns encode a timestamp.

    Inputs:
        path        : str - file to read
        fmt         : TableFormat or None - reuse a previous inspection
        max_columns : int or None - keep only the first ``max_columns`` value
                      columns (the remaining ones are ignored)

    Outputs:
        dict with

            epoch  : np.ndarray[float] - POSIX seconds, strictly increasing
            values : np.ndarray[float] - shape (n_records, n_data), NaN where
                     the file contained a no-data marker
            fmt    : TableFormat
            report : dict - sampling diagnostics, see
                     :func:`iwtime.sampling_report`
            dropped: int  - records removed because their timestamp could not
                     be parsed or was duplicated

    Notes:
        Records are sorted by time and duplicated timestamps are discarded,
        which guarantees the strictly increasing time vector required by the
        interpolation operators.
    """
    fmt = fmt or inspect(path)
    text = _load_text(fmt.path)

    if fmt.layout in ('iso', 'date_time'):
        stamps, values = _parse_text_time(text, fmt)
    else:
        stamps, values = _parse_numeric_time(text, fmt)

    if max_columns is not None and values.shape[1] > max_columns:
        values = values[:, :max_columns]

    epoch = tmb.to_epoch(stamps)
    valid = np.isfinite(epoch)
    dropped = int(np.count_nonzero(~valid))

    epoch = epoch[valid]
    values = values[valid]

    before = epoch.size
    epoch, values = tmb.sort_unique(epoch, values)
    dropped += before - epoch.size

    return {'epoch': epoch, 'values': values, 'fmt': fmt,
            'report': tmb.sampling_report(epoch, os.path.basename(fmt.path)),
            'dropped': dropped}


def _load_text(path):
    """Read the whole file into memory as text (files are at most tens of MB)."""
    with open(path, 'r', errors='replace') as handle:
        return handle.read()


def _parse_numeric_time(text, fmt):
    """
    Parse a fully numeric file (calendar columns followed by values).

    A single ``np.genfromtxt`` call replaces the per-cell Python loop used in
    version 1.x, which dominated the loading time of long records.
    """
    stream = io.StringIO(text)
    table = np.genfromtxt(stream,
                          delimiter=fmt.delimiter,
                          skip_header=1 if fmt.has_header else 0,
                          dtype=float,
                          missing_values=_NA_TOKENS,
                          filling_values=np.nan,
                          invalid_raise=False,
                          usemask=False,
                          encoding='utf-8',
                          comments='#')

    table = np.atleast_2d(table)
    if table.shape[1] < fmt.n_time:
        raise IOError(f"{os.path.basename(fmt.path)}: expected at least "
                      f"{fmt.n_time} columns, found {table.shape[1]}")

    second = table[:, 5] if fmt.layout == 'ymdhms' else None
    stamps = tmb.components_to_datetime64(table[:, 0], table[:, 1],
                                          table[:, 2], table[:, 3],
                                          table[:, 4], second)

    # A NaN in any calendar column makes the timestamp meaningless
    bad = ~np.all(np.isfinite(table[:, :fmt.n_time]), axis=1)
    stamps = stamps.astype('datetime64[s]')
    values = table[:, fmt.n_time:]
    if np.any(bad):
        values = values.copy()
        values[bad] = np.nan
        # mark the timestamp itself as unusable
        stamps = np.where(bad, np.datetime64('NaT', 's'), stamps)

    return stamps, np.atleast_2d(values)


def _parse_text_time(text, fmt):
    """Parse a file whose timestamp is an ISO string (one or two columns)."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if fmt.has_header:
        lines = lines[1:]

    stamps = np.empty(len(lines), dtype='datetime64[s]')
    values = np.full((len(lines), fmt.n_data), np.nan, dtype=float)

    for i, line in enumerate(lines):
        fields = _split(line, fmt.delimiter)
        try:
            if fmt.layout == 'iso':
                token = fields[0].strip().replace('/', '-').replace(' ', 'T')
            else:
                token = (fields[0].strip().replace('/', '-') + 'T'
                         + fields[1].strip())
            stamps[i] = np.datetime64(token, 's')
        except (ValueError, IndexError):
            stamps[i] = np.datetime64('NaT', 's')

        for j in range(fmt.n_time, min(len(fields), fmt.n_columns)):
            values[i, j - fmt.n_time] = _to_float(fields[j])

    return stamps, values


# =============================================================================
#  Input files of Interwave Analyzer
# =============================================================================

def read_temperature(path):
    """
    Read the thermistor-chain record (``.tem``).

    Inputs:
        path : str - temperature file

    Outputs:
        dict with

            epoch  : np.ndarray[float] - POSIX seconds
            temp   : np.ndarray[float] - shape (n_times, n_sensors), degC,
                     ordered from the sensor closest to the surface to the
                     deepest one
            n_sensors : int
            fmt, report, dropped : see :func:`read_timed_table`

    Notes:
        The number of sensors is taken from the *data* records instead of from
        the header line.  Version 1.x counted the words of the header, so a
        header such as ``Wind Speed`` (two words for one column) produced a
        wrong sensor count.
    """
    table = read_timed_table(path)
    table['temp'] = table.pop('values')
    table['n_sensors'] = int(table['temp'].shape[1])

    if table['n_sensors'] < 2:
        raise IOError(f"{os.path.basename(table['fmt'].path)}: at least two "
                      f"temperature columns are required, found "
                      f"{table['n_sensors']}")
    return table


def read_meteorology(path, radiation=0):
    """
    Read the meteorological record (``.met``).

    Inputs:
        path      : str - meteorological file
        radiation : int - 1 when a short-wave radiation column is present

    Outputs:
        dict with

            epoch     : np.ndarray[float]
            speed     : np.ndarray[float] - wind speed (m/s)
            direction : np.ndarray[float] - wind direction (deg)
            radiation : np.ndarray[float] - short-wave radiation (W/m2) or an
                        array of NaN when not supplied
            fmt, report, dropped

    Notes:
        The radiation column is read whenever it exists in the file, even if
        the user did not request the radiation analysis, so that the diagnosis
        file can report the discrepancy.
    """
    table = read_timed_table(path)
    values = table.pop('values')

    if values.shape[1] < 2:
        raise IOError(f"{os.path.basename(table['fmt'].path)}: wind speed and "
                      f"wind direction columns are required")

    table['speed'] = values[:, 0]
    table['direction'] = values[:, 1] % 360.0
    table['has_radiation'] = bool(values.shape[1] >= 3)
    table['radiation'] = (values[:, 2] if table['has_radiation']
                          else np.full(values.shape[0], np.nan))
    table['radiation_requested'] = bool(radiation)

    return table


def read_level(path):
    """
    Read the water-level record (``.niv``).

    Inputs:
        path : str - water-level file

    Outputs:
        dict with ``epoch``, ``level`` (m) and the usual diagnostics.
    """
    table = read_timed_table(path)
    values = table.pop('values')
    if values.shape[1] < 1:
        raise IOError(f"{os.path.basename(table['fmt'].path)}: a water-level "
                      f"column is required")
    table['level'] = values[:, 0]
    return table


def read_sensors(path, n_expected=None):
    """
    Read the thermistor mounting file (``.sen``).

    Inputs:
        path       : str - sensor file
        n_expected : int or None - number of temperature columns found in the
                     temperature file; used for consistency checking

    Outputs:
        dict with

            position      : np.ndarray[float] - distance of every sensor to its
                            reference (m)
            specification : np.ndarray[int]   - 1 = distance measured downwards
                            from the water surface (floating chain),
                            2 = distance measured upwards from the bed
                            (fixed chain)
            messages      : list[str]

    Raises:
        IOError - the file cannot be read, or the sensor count disagrees with
        the temperature file.
    """
    fmt = inspect(path, timed=False)
    stream = io.StringIO(_load_text(fmt.path))
    table = np.genfromtxt(stream, delimiter=fmt.delimiter,
                          skip_header=1 if fmt.has_header else 0,
                          dtype=float, filling_values=np.nan,
                          invalid_raise=False, usemask=False,
                          encoding='utf-8', comments='#')
    table = np.atleast_2d(table)

    if table.shape[1] < 3:
        raise IOError(f"{os.path.basename(fmt.path)}: expected three columns "
                      f"(sensor, position, specification), found "
                      f"{table.shape[1]}")

    position = table[:, 1].astype(float)
    specification = np.rint(table[:, 2]).astype(int)

    messages = []
    unknown = ~np.isin(specification, (1, 2))
    if np.any(unknown):
        messages.append(f"{os.path.basename(fmt.path)}: "
                        f"{int(np.count_nonzero(unknown))} sensor(s) with an "
                        f"unknown specification; assumed to be 1 (measured "
                        f"from the water surface)")
        specification[unknown] = 1

    if n_expected is not None and position.size != n_expected:
        raise IOError(f"{os.path.basename(fmt.path)}: the file describes "
                      f"{position.size} sensors but the temperature file "
                      f"contains {n_expected} temperature columns")

    if not np.all(np.isfinite(position)):
        raise IOError(f"{os.path.basename(fmt.path)}: non-numeric sensor "
                      f"position")

    return {'position': position, 'specification': specification,
            'fmt': fmt, 'messages': messages}


# =============================================================================
#  Basin geometry
# =============================================================================

def loadLen(filepath):
    """
    Read a transect file (``.len``) with the columns {H(m), Ls(m), Xr(m)}.

    Inputs:
        filepath : str - transect file

    Outputs:
        dict with ``depths``, ``dists`` and ``refs`` as NumPy arrays, sorted by
        increasing depth.

    Notes:
        Sorting is applied here so that every downstream consumer (slope
        computation, hypsography reconstruction, interpolation) can rely on a
        monotonic depth axis.
    """
    fmt = inspect(filepath, timed=False)
    stream = io.StringIO(_load_text(fmt.path))
    table = np.genfromtxt(stream, delimiter=fmt.delimiter,
                          skip_header=1 if fmt.has_header else 0,
                          dtype=float, filling_values=np.nan,
                          invalid_raise=False, usemask=False,
                          encoding='utf-8', comments='#')
    table = np.atleast_2d(table)

    if table.shape[1] < 2:
        raise IOError(f"{os.path.basename(fmt.path)}: expected at least the "
                      f"columns H(m) and Ls(m)")

    depths = table[:, 0].astype(float)
    dists = table[:, 1].astype(float)
    refs = (table[:, 2].astype(float) if table.shape[1] >= 3
            else np.zeros_like(depths))

    keep = np.isfinite(depths) & np.isfinite(dists)
    depths, dists, refs = depths[keep], dists[keep], np.nan_to_num(refs[keep])

    order = np.argsort(depths, kind='stable')

    return {"depths": depths[order], "dists": dists[order],
            "refs": refs[order]}


def loadData(type_length, fna, fna_trans, len_basin, change_basin, maxDepth):
    """
    Build the basin geometry according to the selected configuration.

    Inputs:
        type_length  : int   - 1 uniform basin, 2 single transect,
                       3 two orthogonal transects
        fna          : str   - main (longitudinal) transect file
        fna_trans    : str   - transverse transect file
        len_basin    : float - basin length for the uniform case (m)
        change_basin : float - azimuth of the main transect (deg)
        maxDepth     : float - mean water depth, used for the uniform case (m)

    Outputs:
        (longData, transData) - geometry dictionaries, ``transData`` is None
        unless ``type_length`` is 3.
    """
    longData = None
    transData = None

    if type_length == 1:
        # Uniform basin: a vertical wall of constant length
        dz = 0.2
        depths = np.linspace(dz, maxDepth + dz, 4)
        longData = {
            "depths": depths,
            "dists": np.full_like(depths, len_basin, dtype=float),
            "refs": np.zeros_like(depths, dtype=float),
            "orientation": 270,
        }

    elif type_length == 2:
        longData = loadLen(fna)
        longData["orientation"] = 270

    elif type_length == 3:
        longData = loadLen(fna)
        orientation_main = change_basin if change_basin else 270
        longData = mod.basinOrientation(longData, orientation_main)

        transData = loadLen(fna_trans)
        transData = mod.basinOrientation(transData,
                                         (orientation_main + 90) % 360)

    return longData, transData


def lengthType(type_length, fna, additional_params):
    """
    Resolve the final basin configuration from the additional parameters.

    Inputs:
        type_length       : int  - configuration selected in the GUI
        fna               : str  - main transect file
        additional_params : dict - additional-parameter block

    Outputs:
        (type_length, fna, fna_trans, change_basin, warnTransv, missFile)

    Notes:
        A transverse transect supplied through ``pathBathy`` promotes a
        type 2 configuration to type 3.  The path is resolved relative to the
        main transect file, and also accepted as an absolute path.
    """
    fna_trans = None
    warnTransv = False
    missFile = None
    change_basin = 270

    if type_length == 2:
        path_bathy = pars.extractPathBathy(additional_params)

        if path_bathy:
            candidates = [path_bathy,
                          os.path.join(os.path.dirname(_clean_path(fna) or ''),
                                       path_bathy)]
            found = next((c for c in candidates if os.path.isfile(c)), None)

            if found:
                fna_trans = found
                type_length = 3
                change_basin = pars.extractChangeBasin(additional_params)
            else:
                warnTransv = True
                missFile = path_bathy

    return type_length, fna, fna_trans, change_basin, warnTransv, missFile


# =============================================================================
#  Assembly of the analysis data set
# =============================================================================

def sensor_depths(level, position, specification, z0):
    """
    Convert sensor mounting distances into height above the vertical datum.

    Inputs:
        level         : np.ndarray[float] - water-surface elevation at every
                        analysis time (m), shape (n_times,)
        position      : np.ndarray[float] - mounting distance of every sensor
                        (m), shape (n_sensors,)
        specification : np.ndarray[int]   - 1 measured from the surface,
                        2 measured from the bed
        z0            : float - elevation of the vertical datum (m)

    Outputs:
        np.ndarray[float] - shape (n_times, n_sensors), height of every sensor
        above the datum.

    Notes:
        For a floating chain (specification 1) the sensor follows the water
        surface, so its height above the datum changes with the water level;
        for a fixed chain (specification 2) it does not.  The computation is
        fully vectorised, replacing the double Python loop of version 1.x.
    """
    level = np.asarray(level, dtype=float)[:, None]
    position = np.asarray(position, dtype=float)[None, :]
    specification = np.asarray(specification)[None, :]

    return np.where(specification == 1, level - position, z0 + position)


def order_profile(height, temp):
    """
    Sort each temperature profile by decreasing sensor height.

    Inputs:
        height : np.ndarray[float] - shape (n_times, n_sensors), height above
                 datum
        temp   : np.ndarray[float] - shape (n_times, n_sensors)

    Outputs:
        (height_sorted, temp_sorted, n_reordered)
            n_reordered : int - number of profiles whose sensor order had to be
            changed.

    Notes:
        Version 1.x attempted this with a single pass of adjacent swaps that
        (i) indexed ``h[t][-1]`` when the loop variable was zero, thereby
        comparing against the *last* sensor of the profile, and (ii) modified
        the temperature array in place through an alias, so that the "raw" and
        "reordered" temperature fields were in fact the same object.  Here the
        profile is sorted properly and the inputs are never modified.
    """
    height = np.asarray(height, dtype=float)
    temp = np.asarray(temp, dtype=float)

    order = np.argsort(-height, axis=1, kind='stable')
    rows = np.arange(height.shape[0])[:, None]

    height_sorted = height[rows, order]
    temp_sorted = temp[rows, order]

    changed = int(np.count_nonzero(
        np.any(order != np.arange(height.shape[1])[None, :], axis=1)))

    return height_sorted, temp_sorted, changed


def build_dataset(temp_path, meteo_path, sensor_path, level_path,
                  radiation=0, level_mode=1, level_value=0.0, z0=0.0,
                  options=None):
    """
    Read every input file and assemble the data set on a common time base.

    Inputs:
        temp_path   : str   - temperature file (``.tem``)
        meteo_path  : str   - meteorological file (``.met``)
        sensor_path : str   - sensor file (``.sen``)
        level_path  : str   - water-level file (``.niv``), used when
                      ``level_mode`` is 2
        radiation   : int   - 1 when short-wave radiation is analysed
        level_mode  : int   - 1 constant water level, 2 water-level record
        level_value : float - constant water-surface elevation (m)
        z0          : float - elevation of the vertical datum (m)
        options     : dict or None - analysis-window options, see
                      :func:`iwtime.build_timebase`.  Recognised keys:
                      ``start``, ``end``, ``step``, ``trim_to_meteo``,
                      ``max_gap``, ``wind_align``, ``fill_gap``

    Outputs:
        dict with

            epoch      : np.ndarray[float] - analysis grid (POSIX seconds)
            datetime   : list[datetime]    - the same grid as Python datetimes
            dt          : float - analysis step (h), the unit used everywhere
                          else in the program
            temp        : np.ndarray - raw temperature on the analysis grid
            temp_sorted : np.ndarray - profile sorted from surface to bed
            height      : np.ndarray - sensor height above datum, sorted
            level       : np.ndarray - water-surface elevation (m)
            wind, wind_stress_speed, direction, radiation : np.ndarray
            n_sensors   : int
            messages    : list[str] - every diagnostic produced while loading
            quality     : dict      - coverage statistics per variable

    Notes:
        This function is the single entry point used by the backend; it
        replaces the former sequence ``file_len`` / ``temper_read`` /
        ``serial_cota`` / ``serial_wind``.
    """
    options = dict(options or {})
    messages = []

    # ---- temperature: defines the master clock ------------------------------
    temperature = read_temperature(temp_path)
    messages.append(temperature['fmt'].describe())
    messages.extend(temperature['report']['messages'])
    if temperature['dropped']:
        messages.append(f"temperature: {temperature['dropped']} record(s) "
                        f"discarded (unparsable or duplicated timestamp)")

    sensors = read_sensors(sensor_path, temperature['n_sensors'])
    messages.append(sensors['fmt'].describe())
    messages.extend(sensors['messages'])

    # ---- meteorology --------------------------------------------------------
    meteo = read_meteorology(meteo_path, radiation)
    messages.append(meteo['fmt'].describe())
    messages.extend(meteo['report']['messages'])
    if meteo['dropped']:
        messages.append(f"meteorology: {meteo['dropped']} record(s) discarded")
    if radiation and not meteo['has_radiation']:
        messages.append("radiation analysis requested but the meteorological "
                        "file has no radiation column; radiation disabled")
    if meteo['has_radiation'] and not radiation:
        messages.append("the meteorological file contains a radiation column "
                        "that was not selected for analysis")

    # ---- water level --------------------------------------------------------
    level_record = None
    level_span = None
    if level_mode == 2:
        level_record = read_level(level_path)
        messages.append(level_record['fmt'].describe())
        messages.extend(level_record['report']['messages'])
        level_span = (level_record['epoch'][0], level_record['epoch'][-1])

    # ---- analysis time base -------------------------------------------------
    timebase = tmb.build_timebase(
        temperature['epoch'],
        meteo_span=(meteo['epoch'][0], meteo['epoch'][-1]),
        level_span=level_span,
        start=options.get('start'),
        end=options.get('end'),
        step=options.get('step'),
        trim_to_meteo=bool(options.get('trim_to_meteo', False)))
    messages.extend(timebase['messages'])

    epoch = timebase['epoch']
    dt_seconds = timebase['dt']

    if timebase['index'] is not None:
        temp = temperature['temp'][timebase['index']]
    else:
        temp = tmb.resample_matrix(temperature['epoch'], temperature['temp'],
                                   epoch)

    # ---- water level on the analysis grid -----------------------------------
    max_gap = options.get('max_gap')
    if level_mode == 2:
        level, level_mode_used = tmb.align_scalar(
            level_record['epoch'], level_record['level'], epoch,
            mode='auto', max_gap=max_gap, dt_dst=dt_seconds)
        messages.append(f"water level combined by '{level_mode_used}'")
        n_missing = int(np.count_nonzero(~np.isfinite(level)))
        if n_missing:
            messages.append(f"water level: {n_missing} sample(s) outside the "
                            f"record; filled with the record mean")
            level = np.where(np.isfinite(level), level,
                             np.nanmean(level_record['level']))
    else:
        level = np.full(epoch.size, float(level_value))

    # ---- wind on the analysis grid ------------------------------------------
    wind = tmb.align_wind(meteo['epoch'], meteo['speed'], meteo['direction'],
                          epoch, mode=options.get('wind_align', 'auto'),
                          max_gap=max_gap, dt_dst=dt_seconds)
    messages.append(f"wind combined by '{wind['mode']}' "
                    f"(meteorological step "
                    f"{meteo['report']['dt'] / 60.0:.1f} min, analysis step "
                    f"{dt_seconds / 60.0:.1f} min)")

    fill_gap = int(options.get('fill_gap', 0))
    if fill_gap > 0:
        for key in ('speed', 'speed_rms'):
            wind[key], filled, left = tmb.fill_short_gaps(wind[key], fill_gap)
            if filled:
                messages.append(f"wind {key}: {filled} sample(s) interpolated "
                                f"across short gaps, {left} left undefined")

    if radiation and meteo['has_radiation']:
        solar, _ = tmb.align_scalar(meteo['epoch'], meteo['radiation'], epoch,
                                    mode=options.get('wind_align', 'auto'),
                                    max_gap=max_gap, dt_dst=dt_seconds)
    else:
        solar = np.zeros(epoch.size, dtype=float)

    # ---- sensor geometry ----------------------------------------------------
    height = sensor_depths(level, sensors['position'],
                           sensors['specification'], z0)
    height_sorted, temp_sorted, reordered = order_profile(height, temp)
    if reordered:
        messages.append(f"{reordered} profile(s) had to be reordered so that "
                        f"sensor height decreases with column index")

    # ---- quality report -----------------------------------------------------
    quality = {
        'n_times': int(epoch.size),
        'n_sensors': int(temperature['n_sensors']),
        'temperature_missing': int(np.count_nonzero(~np.isfinite(temp))),
        'wind_missing': int(np.count_nonzero(~np.isfinite(wind['speed']))),
        'direction_missing': int(np.count_nonzero(
            ~np.isfinite(wind['direction']))),
        'wind_mode': wind['mode'],
        'dt_temperature': float(temperature['report']['dt']),
        'dt_meteorology': float(meteo['report']['dt']),
        'dt_analysis': float(dt_seconds),
        'resampled': bool(timebase['resampled']),
    }

    return {
        'epoch': epoch,
        'datetime': tmb.to_pydatetime(tmb.from_epoch(epoch)),
        'dt': dt_seconds / 3600.0,          # hours, as used by the backend
        'dt_seconds': dt_seconds,
        'temp': temp,
        'temp_sorted': temp_sorted,
        'height': height_sorted,
        'height_raw': height,
        'level': level,
        'wind': wind['speed'],
        'wind_stress_speed': wind['speed_rms'],
        'direction': wind['direction'],
        'radiation': solar,
        'has_radiation': bool(meteo['has_radiation'] and radiation),
        'n_sensors': int(temperature['n_sensors']),
        'sensor_position': sensors['position'],
        'sensor_specification': sensors['specification'],
        'messages': messages,
        'quality': quality,
    }


# =============================================================================
#  Backwards compatible helpers
# =============================================================================

def file_len(fname, col='off'):
    """
    Number of data records (and optionally of value columns) of a file.

    Retained for backwards compatibility; new code should use
    :func:`inspect` together with :func:`read_timed_table`.
    """
    fmt = inspect(fname)
    with open(fmt.path, 'r', errors='replace') as handle:
        count = sum(1 for line in handle if line.strip())
    if fmt.has_header:
        count -= 1
    return (count, fmt.n_data) if col == 'on' else count


def find_nearest(array, value):
    """Index of the element of ``array`` closest to ``value``."""
    return int(np.abs(np.asarray(array) - value).argmin())


def level(sen, qt):
    """
    Legacy reader of the sensor file, returning ``(position, specification)``.
    """
    record = read_sensors(sen, None)
    return record['position'][:qt], record['specification'][:qt]
