"""Setuptools config for access-atmosphere."""

from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("rtm", ["lib/rtm.pyx"], libraries=["gfortran", "m"])
]

setup(
    name="access-atmosphere",
    version="0.0.1",
    packages=["access-atmosphere"],
    python_requires=">=3.9",
    install_requires=[
        "cdsapi",
    ],
    author="Richard Lindsley",
    author_email="lindsley@remss.com",
    url="http://gitlab.remss.com/access/atmospheric-rtm/",
    include_package_data=True,
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
    zip_safe=False,
)
