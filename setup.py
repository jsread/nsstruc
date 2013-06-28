import os
from setuptools import setup

# Everything I know about documentation I learned from
# http://pythonhosted.org/an_example_pypi_project/
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
        return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
            name = "nsstruc",
            version = "0.0.1",
            author = "Jocelyn Read",
            author_email = "jread@fullerton.edu",
            description = ("Some python code for neutron stars and EOS "
                           ),
            license = "BSD",
            keywords = "scientific scipy",
            url = "",
            packages=['eos','tov', 'tests'],
            long_description=read('README'),
            classifiers=[
                        "Development Status :: 1 - Planning",
                        "License :: OSI Approved :: BSD License",
                    ],
    )
