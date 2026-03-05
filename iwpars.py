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

def extractParam(additional_params, key, cast_type=str, default=None):

    raw = additional_params.get(key)

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
    return extractParam(
        additional_params,
        key="pathBathy",
        cast_type=str,
        default=None
    )

def extractChangeBasin(additional_params):
    return extractParam(
        additional_params,
        key="changeBasin",
        cast_type=float,
        default=270
    )