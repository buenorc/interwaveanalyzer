# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Additional parameters module

Interwave Analyzer - Version 2 (2026) 
Additional parameters module version: 2.260305

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

"""

# ------------- Additional Parameters Configuration ---------------------------

"""
Generic parameter extractor.

    - Looks for 'key' in additional_params
    - Removes surrounding '#'
    - Casts to desired type
    - Returns default if missing or invalid
"""

def extractParam(additional_params, key, cast_type=str, default=None,
                 aliases=()):

    raw = additional_params.get(key)

    # Accept documented spelling variants of the same keyword
    for alias in aliases:
        if raw:
            break
        raw = additional_params.get(alias)

    # Keyword matching is case-insensitive: a parameter typed as "namebasin"
    # or "NameBasin" used to be ignored without any message.
    if not raw:
        wanted = {key.lower()} | {a.lower() for a in aliases}
        for name, value in additional_params.items():
            if name.strip().lower() in wanted:
                raw = value
                break

    if not raw:
        return default

    cleaned = raw.strip().strip('#').strip()

    if not cleaned:
        return default

    try:
        return cast_type(cleaned)
    except (ValueError, TypeError):
        return default


# ------------- Additional Parameters -----------------------------------------

def extractName(additional_params):
    return extractParam(
        additional_params,
        key="nameBasin",
        cast_type=str,
        default="No name provided"
    )

def extractPathBathy(additional_params):
    # "pathBathymetry" is accepted as well: it is the spelling used in some of
    # the distributed example data sets, and the mismatch previously caused the
    # transverse transect to be ignored without any message.
    return extractParam(
        additional_params,
        key="pathBathy",
        cast_type=str,
        default=None,
        aliases=("pathBathymetry", "pathTransect")
    )

def extractChangeBasin(additional_params):
    return extractParam(
        additional_params,
        key="changeBasin",
        cast_type=float,
        default=270
    )


# ------------- Analysis window, resampling and parametrisations --------------

"""
The parameters below were introduced in version 2.260810.  They control the
time base of the analysis and the way the auxiliary records are combined with
the temperature record.  All of them are optional: when none is supplied the
analysis behaves exactly as before, i.e. it runs over the whole temperature
record at the temperature sampling step.

    timeStart      ISO date, first instant of the analysis window
    timeEnd        ISO date, last instant of the analysis window
    timeStep       minutes, analysis step (coarser than the record only)
    trimToMeteo    1 to restrict the window to the meteorological record
    maxGap         minutes, longest interval over which an auxiliary record
                   may be extended before it is reported as missing
    fillGap        number of consecutive missing analysis steps that may be
                   interpolated in the wind record
    windAlign      auto | interp | nearest | mean, operator used to bring the
                   meteorological record onto the analysis grid
    dragLaw        smooth | legacy | largepond, surface drag parametrisation
"""


def _extractDate(additional_params, key):
    """Parse an ISO date parameter into POSIX seconds, or None."""
    import numpy as np

    raw = extractParam(additional_params, key=key, cast_type=str, default=None)
    if not raw:
        return None

    token = raw.strip().replace('/', '-').replace(' ', 'T')
    try:
        return float((np.datetime64(token, 's')
                      - np.datetime64('1970-01-01T00:00:00', 's'))
                     / np.timedelta64(1, 's'))
    except ValueError:
        return None


def extractTimeOptions(additional_params):
    """
    Collect every option that controls the analysis time base.

    Inputs:
        additional_params : dict - additional-parameter block from the GUI

    Outputs:
        dict accepted by ``iwload.build_dataset`` as its ``options`` argument.
    """
    step = extractParam(additional_params, key="timeStep", cast_type=float,
                        default=None)
    max_gap = extractParam(additional_params, key="maxGap", cast_type=float,
                           default=None)

    align = extractParam(additional_params, key="windAlign", cast_type=str,
                         default="auto").strip().lower()
    if align not in ("auto", "interp", "nearest", "mean"):
        align = "auto"

    return {
        'start': _extractDate(additional_params, "timeStart"),
        'end': _extractDate(additional_params, "timeEnd"),
        'step': step * 60.0 if step else None,
        'trim_to_meteo': bool(extractParam(additional_params,
                                           key="trimToMeteo",
                                           cast_type=int, default=0)),
        'max_gap': max_gap * 60.0 if max_gap else None,
        'fill_gap': extractParam(additional_params, key="fillGap",
                                 cast_type=int, default=0),
        'wind_align': align,
    }


def extractDragLaw(additional_params):
    """
    Surface drag parametrisation selected by the user.

    Outputs:
        str - 'smooth' (default), 'legacy' or 'largepond'
    """
    law = extractParam(additional_params, key="dragLaw", cast_type=str,
                       default="smooth").strip().lower()
    return law if law in ("smooth", "legacy", "largepond") else "smooth"


def extractMixingEfficiency(additional_params):
    """
    Fraction of the wind power that reaches the stratified interior.

    Outputs:
        float - default 0.0016
    """
    value = extractParam(additional_params, key="mixingEfficiency",
                         cast_type=float, default=0.0016)
    return value if (value and 0 < value < 1) else 0.0016