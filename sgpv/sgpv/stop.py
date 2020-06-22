# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 03:20:06 2020

@author: Kristjan

A helper module to emulate some functions found in R used to display warning
and error messages
"""


def stop(text: str):
    
    """
    Emulates the behavior of R's stop-function
    -> not test outside of IPython-console -> color only working when runnning
    without colorama

    Parameters
    ----------
    text : str
        DESCRIPTION.

    Returns
    -------
    None.

    """
    #from colorama import init
    from termcolor import colored
    import sys
    #init()
    print(colored(text, 'red'))
    sys.exit(1)


def warning(text: str):
    """
    Emulate the behavior of the R-function with the same name

    Parameters
    ----------
    text : str
        DESCRIPTION.

    Returns
    -------
    None.
    """
    from termcolor import colored
    print(colored(text, 'red'))
    
     