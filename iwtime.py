# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno
Interwave Analyzer - Time-base and resampling module

Interwave Analyzer - Version 2 (2026)
Time-base module version: 2.260810

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

Up to version 1.x every time series was mapped onto the temperature record
through an integer "serial number" built as

    serial = minute + 60*hour + 1440*day + 43200*month + 518400*year

which assumes 30-day months and 12*30-day years.  That mapping is *not*
monotonic in physical time: 31 January and 1 February share exactly the same
serial number, so meteorological and water-level records were mismatched by up
to one day around every long month.  This module replaces that construction by
true calendar arithmetic (``numpy.datetime64`` -> POSIX seconds) and adds the
resampling operators required to combine records that were not sampled on the
same clock.

The module is deliberately free of any dependency other than NumPy so that it
can be imported by the loader before SciPy is available.
"""

import numpy as np

# Reference epoch used to convert datetime64 to float seconds
_EPOCH = np.datetime64('1970-01-01T00:00:00', 's')

# Sentinel used throughout Interwave Analyzer for "not provided"
NODATA = -999


# =============================================================================
#  Calendar conversion
# =============================================================================

def components_to_datetime64(year, month, day, hour, minute, second=None):
    """
    Build a ``datetime64[s]`` vector from separated calendar components.

    Inputs:
        year, month, day, hour, minute : array_like
            Calendar components.  They may be float arrays (as produced by a
            plain text reader); they are rounded to the nearest integer.
        second : array_like or None
            Optional seconds column.

    Outputs:
        np.ndarray[datetime64[s]] - absolute time of every record.

    Notes:
        Calendar arithmetic is delegated to NumPy, therefore leap years and
        variable month lengths are handled exactly.  Out-of-range components
        (e.g. hour = 24) are normalised by the timedelta additions below, which
        makes the reader tolerant to loggers that write 24:00 instead of 00:00
        of the following day.
    """
    year = np.rint(np.asarray(year, dtype=float)).astype('int64')
    month = np.rint(np.asarray(month, dtype=float)).astype('int64')
    day = np.rint(np.asarray(day, dtype=float)).astype('int64')
    hour = np.rint(np.asarray(hour, dtype=float)).astype('int64')
    minute = np.rint(np.asarray(minute, dtype=float)).astype('int64')

    if second is None:
        second = np.zeros_like(minute)
    else:
        second = np.rint(np.asarray(second, dtype=float)).astype('int64')

    # Build the first day of the month, then add the remaining offsets.  Using
    # timedelta additions (instead of composing a string) keeps the operation
    # vectorised and normalises out-of-range components.  Casting an integer
    # array to 'datetime64[Y]' interprets it as an offset from 1970, hence the
    # explicit subtraction.
    base = (year - 1970).astype('datetime64[Y]')
    base = base.astype('datetime64[M]') + (month - 1).astype('timedelta64[M]')

    stamp = (base.astype('datetime64[D]')
             + (day - 1).astype('timedelta64[D]'))
    stamp = (stamp.astype('datetime64[s]')
             + hour.astype('timedelta64[h]')
             + minute.astype('timedelta64[m]')
             + second.astype('timedelta64[s]'))

    return stamp


def to_epoch(stamp):
    """
    Convert ``datetime64`` to POSIX seconds (float64).

    Inputs:
        stamp : np.ndarray[datetime64] - absolute times.

    Outputs:
        np.ndarray[float] - seconds since 1970-01-01T00:00:00.
    """
    stamp = np.asarray(stamp, dtype='datetime64[s]')
    return (stamp - _EPOCH) / np.timedelta64(1, 's')


def from_epoch(seconds):
    """
    Convert POSIX seconds back to ``datetime64[s]``.

    Inputs:
        seconds : array_like - seconds since the epoch.

    Outputs:
        np.ndarray[datetime64[s]]
    """
    seconds = np.rint(np.asarray(seconds, dtype=float)).astype('int64')
    return _EPOCH + seconds.astype('timedelta64[s]')


def to_pydatetime(stamp):
    """
    Convert ``datetime64`` to a list of ``datetime.datetime`` objects.

    Matplotlib date axes and the dashboard both expect native datetimes, so the
    conversion is centralised here.

    Inputs:
        stamp : np.ndarray[datetime64]

    Outputs:
        list[datetime.datetime]
    """
    stamp = np.asarray(stamp, dtype='datetime64[us]')
    return stamp.astype('O').tolist()


# =============================================================================
#  Sampling diagnostics
# =============================================================================

def sampling_report(epoch, name='series'):
    """
    Characterise the sampling of a time vector.

    Inputs:
        epoch : np.ndarray[float] - POSIX seconds, expected to be increasing.
        name  : str - label used in the diagnostic messages.

    Outputs:
        dict with the keys

            name        : str
            n           : int    - number of records
            dt          : float  - representative time step (s), the median of
                                   the positive forward differences
            dt_min/dt_max : float - extreme steps (s)
            regular     : bool   - True when every step matches ``dt`` to
                                   within 1 % (or 1 s, whichever is larger)
            duplicated  : int    - number of non-increasing steps
            gaps        : int    - number of steps larger than 1.5*dt
            largest_gap : float  - largest step (s)
            start, end  : float  - first and last epoch second
            messages    : list[str] - human readable diagnostics

    Notes:
        The median (rather than the difference between the first two records,
        as used up to version 1.x) is robust to a single corrupt timestamp at
        the beginning of the file.
    """
    epoch = np.asarray(epoch, dtype=float)
    out = {'name': name, 'n': int(epoch.size), 'messages': []}

    if epoch.size < 2:
        out.update(dt=np.nan, dt_min=np.nan, dt_max=np.nan, regular=False,
                   duplicated=0, gaps=0, largest_gap=np.nan,
                   start=float(epoch[0]) if epoch.size else np.nan,
                   end=float(epoch[-1]) if epoch.size else np.nan)
        out['messages'].append(f"{name}: fewer than two records")
        return out

    diff = np.diff(epoch)
    positive = diff[diff > 0]

    dt = float(np.median(positive)) if positive.size else np.nan
    duplicated = int(np.count_nonzero(diff <= 0))
    gaps = int(np.count_nonzero(diff > 1.5 * dt)) if np.isfinite(dt) else 0

    tolerance = max(0.01 * dt, 1.0) if np.isfinite(dt) else np.inf
    regular = bool(positive.size == diff.size
                   and np.all(np.abs(diff - dt) <= tolerance))

    out.update(dt=dt,
               dt_min=float(diff.min()),
               dt_max=float(diff.max()),
               regular=regular,
               duplicated=duplicated,
               gaps=gaps,
               largest_gap=float(diff.max()),
               start=float(epoch[0]),
               end=float(epoch[-1]))

    if duplicated:
        out['messages'].append(
            f"{name}: {duplicated} non-increasing timestamp(s)")
    if gaps:
        out['messages'].append(
            f"{name}: {gaps} gap(s), largest {diff.max() / 60.0:.1f} min "
            f"against a nominal step of {dt / 60.0:.1f} min")
    if not regular and not duplicated and not gaps:
        out['messages'].append(
            f"{name}: irregular sampling (step {diff.min() / 60.0:.1f} to "
            f"{diff.max() / 60.0:.1f} min)")

    return out


def sort_unique(epoch, *arrays):
    """
    Sort a record by time and drop duplicated timestamps.

    Inputs:
        epoch     : np.ndarray[float] - POSIX seconds
        *arrays   : np.ndarray - companion arrays whose first axis matches
                    ``epoch``

    Outputs:
        (epoch, *arrays) - filtered copies, strictly increasing in time.

    Notes:
        Duplicated timestamps keep their *first* occurrence.  A strictly
        increasing time vector is a precondition for ``np.interp`` and
        ``np.searchsorted`` used by the alignment operators below.
    """
    epoch = np.asarray(epoch, dtype=float)
    order = np.argsort(epoch, kind='stable')
    epoch = epoch[order]
    arrays = [np.asarray(a)[order] for a in arrays]

    keep = np.ones(epoch.size, dtype=bool)
    if epoch.size > 1:
        keep[1:] = np.diff(epoch) > 0

    epoch = epoch[keep]
    arrays = [a[keep] for a in arrays]

    return (epoch, *arrays)


# =============================================================================
#  Alignment operators
# =============================================================================

def _edges_from_centres(t, dt):
    """Interval edges of a target grid whose samples are interval centres."""
    edges = np.empty(t.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (t[:-1] + t[1:])
    edges[0] = t[0] - 0.5 * dt
    edges[-1] = t[-1] + 0.5 * dt
    return edges


def _block_reduce(t_src, values, edges, weights=None):
    """
    Sum ``values`` (optionally weighted) inside the bins defined by ``edges``.

    Returns (sum, count_or_weight_sum).  Implemented with ``np.add.at`` on the
    bin indices, i.e. a single O(N) pass instead of one search per target
    sample.
    """
    idx = np.searchsorted(edges, t_src, side='right') - 1
    inside = (idx >= 0) & (idx < edges.size - 1)
    idx = idx[inside]

    values = np.asarray(values, dtype=float)[inside]
    finite = np.isfinite(values)
    idx = idx[finite]
    values = values[finite]

    total = np.zeros(edges.size - 1, dtype=float)
    norm = np.zeros(edges.size - 1, dtype=float)

    if weights is None:
        np.add.at(total, idx, values)
        np.add.at(norm, idx, 1.0)
    else:
        w = np.asarray(weights, dtype=float)[inside][finite]
        np.add.at(total, idx, values * w)
        np.add.at(norm, idx, w)

    return total, norm


def align_scalar(t_src, v_src, t_dst, mode='auto', max_gap=None, dt_dst=None):
    """
    Map a scalar time series onto a target time grid.

    Inputs:
        t_src   : np.ndarray[float] - source POSIX seconds (increasing)
        v_src   : np.ndarray[float] - source values
        t_dst   : np.ndarray[float] - target POSIX seconds (increasing)
        mode    : {'auto', 'interp', 'nearest', 'mean'}
            'interp'  - linear interpolation in time
            'nearest' - nearest neighbour (legacy behaviour of version 1.x)
            'mean'    - arithmetic mean of every source sample falling inside
                        the target interval, falling back to interpolation for
                        empty intervals
            'auto'    - 'mean' when the source is sampled at least 1.5 times
                        faster than the target, otherwise 'interp'
        max_gap : float or None
            Target samples whose distance to the closest source sample exceeds
            ``max_gap`` seconds are returned as NaN.  ``None`` disables the
            check.  This prevents a single record from being propagated across
            a multi-day instrument outage, which the legacy nearest-neighbour
            search did silently.
        dt_dst  : float or None - target step (s); inferred when omitted.

    Outputs:
        (values, mode_used)
            values    : np.ndarray[float] - resampled series, NaN where the
                        gap criterion is violated
            mode_used : str - the operator actually applied
    """
    t_src = np.asarray(t_src, dtype=float)
    v_src = np.asarray(v_src, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)

    if t_src.size == 0:
        return np.full(t_dst.size, np.nan), 'empty'

    if t_src.size == 1:
        out = np.full(t_dst.size, v_src[0], dtype=float)
        if max_gap is not None:
            out[np.abs(t_dst - t_src[0]) > max_gap] = np.nan
        return out, 'constant'

    if dt_dst is None:
        dt_dst = float(np.median(np.diff(t_dst))) if t_dst.size > 1 else np.nan
    dt_src = float(np.median(np.diff(t_src)))

    if mode == 'auto':
        mode = 'mean' if (np.isfinite(dt_dst) and dt_src * 1.5 <= dt_dst) \
            else 'interp'

    finite = np.isfinite(v_src)

    if mode == 'nearest':
        idx = np.searchsorted(t_src, t_dst)
        idx = np.clip(idx, 1, t_src.size - 1)
        left = t_dst - t_src[idx - 1]
        right = t_src[idx] - t_dst
        idx = np.where(left <= right, idx - 1, idx)
        out = v_src[idx]

    elif mode == 'mean':
        edges = _edges_from_centres(t_dst, dt_dst)
        total, count = _block_reduce(t_src, v_src, edges)
        out = np.divide(total, count, out=np.full_like(total, np.nan),
                        where=count > 0)
        # Intervals without any source sample fall back to interpolation
        empty = ~np.isfinite(out)
        if np.any(empty) and np.any(finite):
            out[empty] = np.interp(t_dst[empty], t_src[finite], v_src[finite])

    else:  # 'interp'
        mode = 'interp'
        if not np.any(finite):
            return np.full(t_dst.size, np.nan), 'interp'
        out = np.interp(t_dst, t_src[finite], v_src[finite])

    out = np.asarray(out, dtype=float)

    if max_gap is not None:
        out[_gap_mask(t_src[finite] if np.any(finite) else t_src,
                      t_dst, max_gap)] = np.nan

    return out, mode


def align_direction(t_src, dir_src, t_dst, speed_src=None, mode='auto',
                    max_gap=None, dt_dst=None):
    """
    Map a directional (circular) time series onto a target time grid.

    Inputs:
        t_src     : np.ndarray[float] - source POSIX seconds
        dir_src   : np.ndarray[float] - direction in degrees (meteorological
                    convention, 0 deg = North, increasing clockwise)
        t_dst     : np.ndarray[float] - target POSIX seconds
        speed_src : np.ndarray[float] or None - magnitude used to weight the
                    average.  When provided the result is the direction of the
                    vector-mean momentum instead of the direction of the
                    unit-vector mean.
        mode, max_gap, dt_dst : see :func:`align_scalar`

    Outputs:
        (direction, mode_used) - direction wrapped to [0, 360).

    Notes:
        Directions must never be interpolated as plain scalars: the linear
        average of 350 deg and 10 deg is 180 deg, i.e. exactly the opposite of
        the correct answer (0 deg).  The series is therefore decomposed into
        its Cartesian components,

            u = V sin(theta),   v = V cos(theta),

        which are averaged or interpolated independently, and recombined with
        ``atan2``.  With ``speed_src`` supplied this is the vector average of
        the wind, which is the quantity that matters for the momentum actually
        transferred to the basin.
    """
    t_src = np.asarray(t_src, dtype=float)
    dir_src = np.asarray(dir_src, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)

    magnitude = np.ones_like(dir_src) if speed_src is None \
        else np.asarray(speed_src, dtype=float)

    theta = np.deg2rad(dir_src)
    u = magnitude * np.sin(theta)
    v = magnitude * np.cos(theta)

    u_dst, used = align_scalar(t_src, u, t_dst, mode=mode, max_gap=max_gap,
                               dt_dst=dt_dst)
    v_dst, _ = align_scalar(t_src, v, t_dst, mode=used, max_gap=max_gap,
                            dt_dst=dt_dst)

    out = np.degrees(np.arctan2(u_dst, v_dst)) % 360.0
    # A tiny negative angle rounds to exactly 360.0 under the modulo; fold it
    # back so that the result always lies in [0, 360).
    out[out >= 360.0] = 0.0
    out[~np.isfinite(u_dst) | ~np.isfinite(v_dst)] = np.nan

    return out, used


def align_wind(t_src, speed_src, dir_src, t_dst, mode='auto', max_gap=None,
               dt_dst=None):
    """
    Map a wind record onto a target grid, preserving momentum flux.

    Inputs:
        t_src     : np.ndarray[float] - source POSIX seconds
        speed_src : np.ndarray[float] - wind speed (m/s)
        dir_src   : np.ndarray[float] - wind direction (degrees)
        t_dst     : np.ndarray[float] - target POSIX seconds
        mode, max_gap, dt_dst : see :func:`align_scalar`

    Outputs:
        dict with

            speed      : np.ndarray - arithmetic mean speed (m/s)
            speed_rms  : np.ndarray - stress-equivalent speed (m/s), see below
            direction  : np.ndarray - speed-weighted vector direction (deg)
            mode       : str        - operator applied
            covered    : np.ndarray[bool] - True where source data exist

    Notes:
        When the meteorological record is sampled faster than the temperature
        record (the usual situation), block averaging the *speed* biases the
        forcing low, because the surface stress scales with the square of the
        wind speed,

            tau = rho_air C_D U^2 .

        The mean stress over an averaging interval is therefore proportional to
        <U^2>, not to <U>^2.  ``speed_rms`` returns

            U_rms = sqrt( <U^2> ) ,

        i.e. the constant wind speed that would deliver the same momentum flux
        over the interval.  It is the quantity used for the wind stress,
        friction velocity, Richardson and Wedderburn numbers, while ``speed``
        is retained for reporting and for the wind spectra, where the linear
        signal is the quantity of interest.
    """
    t_src = np.asarray(t_src, dtype=float)
    speed_src = np.asarray(speed_src, dtype=float)
    dir_src = np.asarray(dir_src, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)

    speed, used = align_scalar(t_src, speed_src, t_dst, mode=mode,
                               max_gap=max_gap, dt_dst=dt_dst)

    speed_sq, _ = align_scalar(t_src, speed_src ** 2, t_dst, mode=used,
                               max_gap=max_gap, dt_dst=dt_dst)
    with np.errstate(invalid='ignore'):
        speed_rms = np.sqrt(np.clip(speed_sq, 0.0, None))

    direction, _ = align_direction(t_src, dir_src, t_dst, speed_src=speed_src,
                                   mode=used, max_gap=max_gap, dt_dst=dt_dst)

    covered = np.isfinite(speed)

    return {'speed': speed, 'speed_rms': speed_rms, 'direction': direction,
            'mode': used, 'covered': covered}


def _gap_mask(t_src, t_dst, max_gap):
    """
    Boolean mask of target samples farther than ``max_gap`` from any source
    sample.
    """
    t_src = np.asarray(t_src, dtype=float)
    if t_src.size == 0:
        return np.ones(np.size(t_dst), dtype=bool)

    idx = np.searchsorted(t_src, t_dst)
    idx_hi = np.clip(idx, 0, t_src.size - 1)
    idx_lo = np.clip(idx - 1, 0, t_src.size - 1)

    distance = np.minimum(np.abs(t_dst - t_src[idx_lo]),
                          np.abs(t_src[idx_hi] - t_dst))
    return distance > max_gap


def fill_short_gaps(values, max_run):
    """
    Interpolate across NaN runs no longer than ``max_run`` samples.

    Inputs:
        values  : np.ndarray[float] - series possibly containing NaN
        max_run : int - longest run of consecutive NaN that may be filled

    Outputs:
        (filled, n_filled, n_left)
            filled   : np.ndarray - copy with short gaps interpolated
            n_filled : int - number of samples filled
            n_left   : int - number of samples still NaN

    Notes:
        Long outages are deliberately left as NaN so that they propagate into
        the diagnostics rather than being disguised as data.
    """
    values = np.asarray(values, dtype=float).copy()
    bad = ~np.isfinite(values)

    if not np.any(bad) or np.all(bad):
        return values, 0, int(np.count_nonzero(bad))

    good_idx = np.flatnonzero(~bad)
    filled = 0

    # Identify contiguous runs of NaN
    edges = np.flatnonzero(np.diff(bad.astype(np.int8)))
    starts = (edges + 1)[bad[edges + 1]]
    if bad[0]:
        starts = np.concatenate(([0], starts))
    ends = (edges + 1)[~bad[edges + 1]]
    if bad[-1]:
        ends = np.concatenate((ends, [bad.size]))

    x = np.arange(values.size, dtype=float)
    for start, stop in zip(starts, ends):
        # Interior runs only: an unbounded run at either end cannot be
        # interpolated without extrapolating.
        if start == 0 or stop == values.size:
            continue
        if stop - start > max_run:
            continue
        values[start:stop] = np.interp(x[start:stop], x[good_idx],
                                       values[good_idx])
        filled += stop - start

    return values, int(filled), int(np.count_nonzero(~np.isfinite(values)))


# =============================================================================
#  Analysis time base
# =============================================================================

def build_timebase(epoch_temp, meteo_span=None, level_span=None,
                   start=None, end=None, step=None, trim_to_meteo=False):
    """
    Define the time grid on which the whole analysis is carried out.

    Inputs:
        epoch_temp    : np.ndarray[float] - POSIX seconds of the temperature
                        record (the master clock of the analysis)
        meteo_span    : (float, float) or None - first/last meteorological
                        record, used when ``trim_to_meteo`` is set
        level_span    : (float, float) or None - first/last water-level record
                        (reported only)
        start, end    : float or None - explicit POSIX bounds requested by the
                        user
        step          : float or None - requested step in seconds; when given
                        the analysis grid is rebuilt uniformly at this step
        trim_to_meteo : bool - restrict the analysis window to the interval
                        actually covered by the meteorological record

    Outputs:
        dict with

            epoch      : np.ndarray[float] - analysis grid (POSIX seconds)
            dt         : float  - analysis step (s)
            index      : np.ndarray[int] or None - indices of ``epoch_temp``
                         retained; ``None`` when the grid was rebuilt and the
                         temperature has to be resampled instead of subset
            resampled  : bool   - True when the grid is not a subset of the
                         temperature record
            messages   : list[str]

    Notes:
        Version 1.x hard-wired the analysis to *every* record of the
        temperature file, starting at its first date.  The window is now an
        explicit object, which makes it possible to (i) start the analysis on
        an arbitrary date, (ii) analyse a sub-period without editing the input
        files, (iii) work on a coarser grid for a quick look, and (iv) avoid
        extrapolating the meteorological forcing beyond its own record.
    """
    epoch_temp = np.asarray(epoch_temp, dtype=float)
    messages = []

    lo = epoch_temp[0]
    hi = epoch_temp[-1]

    if trim_to_meteo and meteo_span is not None:
        new_lo = max(lo, meteo_span[0])
        new_hi = min(hi, meteo_span[1])
        if new_hi > new_lo:
            if new_lo > lo or new_hi < hi:
                messages.append(
                    "analysis window trimmed to the meteorological record")
            lo, hi = new_lo, new_hi
        else:
            messages.append(
                "temperature and meteorological records do not overlap; "
                "trimToMeteo ignored")

    if start is not None:
        if start > hi:
            messages.append("timeStart is after the end of the temperature "
                            "record; ignored")
        else:
            lo = max(lo, start)
    if end is not None:
        if end < lo:
            messages.append("timeEnd is before the start of the analysis "
                            "window; ignored")
        else:
            hi = min(hi, end)

    if level_span is not None:
        if level_span[0] > lo or level_span[1] < hi:
            messages.append(
                "water-level record does not cover the whole analysis window")

    inside = (epoch_temp >= lo) & (epoch_temp <= hi)
    if np.count_nonzero(inside) < 3:
        messages.append("requested window keeps fewer than three temperature "
                        "records; the full record is used instead")
        inside = np.ones(epoch_temp.size, dtype=bool)
        lo, hi = epoch_temp[0], epoch_temp[-1]

    index = np.flatnonzero(inside)
    native = epoch_temp[index]
    dt_native = float(np.median(np.diff(native))) if native.size > 1 else np.nan

    report = sampling_report(native, 'temperature')

    if step is None and report['regular']:
        return {'epoch': native, 'dt': dt_native, 'index': index,
                'resampled': False, 'messages': messages + report['messages']}

    # Either the user asked for a different step, or the native record is not
    # evenly spaced.  In both cases a uniform grid is generated; the caller
    # resamples the temperature onto it.
    dt_target = float(step) if step else dt_native
    if not np.isfinite(dt_target) or dt_target <= 0:
        dt_target = dt_native

    if step and dt_target < dt_native:
        messages.append(
            f"requested timeStep ({dt_target / 60.0:.1f} min) is finer than "
            f"the temperature record ({dt_native / 60.0:.1f} min); "
            f"the record step is used instead")
        dt_target = dt_native

    if report['regular'] and abs(dt_target - dt_native) <= 1e-6:
        # Nothing left to do: the requested grid is the native one.
        return {'epoch': native, 'dt': dt_native, 'index': index,
                'resampled': False, 'messages': messages + report['messages']}

    n = int(np.floor((hi - lo) / dt_target)) + 1
    grid = lo + dt_target * np.arange(n, dtype=float)

    if not report['regular']:
        messages.append(
            f"temperature record resampled onto a uniform "
            f"{dt_target / 60.0:.1f} min grid")
    elif step:
        messages.append(
            f"temperature record resampled from {dt_native / 60.0:.1f} min to "
            f"{dt_target / 60.0:.1f} min as requested")

    return {'epoch': grid, 'dt': dt_target, 'index': None,
            'resampled': True, 'messages': messages + report['messages']}


def resample_matrix(t_src, matrix, t_dst):
    """
    Linearly interpolate every column of a 2-D field onto a new time grid.

    Inputs:
        t_src  : np.ndarray[float] - source POSIX seconds (increasing)
        matrix : np.ndarray - shape (len(t_src), n_columns)
        t_dst  : np.ndarray[float] - target POSIX seconds

    Outputs:
        np.ndarray - shape (len(t_dst), n_columns)
    """
    matrix = np.asarray(matrix, dtype=float)
    t_src = np.asarray(t_src, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)

    if matrix.ndim == 1:
        finite = np.isfinite(matrix)
        if not np.any(finite):
            return np.full(t_dst.size, np.nan)
        return np.interp(t_dst, t_src[finite], matrix[finite])

    out = np.empty((t_dst.size, matrix.shape[1]), dtype=float)
    for j in range(matrix.shape[1]):
        column = matrix[:, j]
        finite = np.isfinite(column)
        if not np.any(finite):
            out[:, j] = np.nan
        else:
            out[:, j] = np.interp(t_dst, t_src[finite], column[finite])
    return out
