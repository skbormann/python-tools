# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 03:20:06 2020

@author: Kristjan
"""

def stop(text):
    """
    Emulates the behavior of R's stop-function
    -> not test outside of IPython-console -> color only working when runnning
    without colorama

    Parameters
    ----------
    text : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    from colorama import init
    from termcolor import colored
    import sys
    init()
    print(colored(text,'red'))
    sys.exit(1)
    