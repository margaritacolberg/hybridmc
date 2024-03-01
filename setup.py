import os
from setuptools import Extension, setup

# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="hybridmc",
    version=os.environ.get('GIT_DESCRIBE_TAG', '0.1.0').strip('v'),
    author="Vigneshwar Rajesh",
    author_email="vignesh.rajesh@mail.utoronto.ca",
    description="A project using c++ for Monte Carlo sampling based simulations.",
    long_description="",
    zip_safe=False,
    python_requires=">=3.8",
)