#!/usr/bin/env python

from distutils.core import setup

setup(name = "fenics-ascot",
      version = "1.2.0",
      description = "(A)utomated (S)tability (CO)ndition (T)ester",
      author = "Marie E. Rognes",
      author_email = "meg@simula.no",
      packages = ["ascot"],
      package_dir={"ascot": "ascot"})
