from setuptools import setup, find_packages

with open("README.md","r") as rfile:
  readme = rfile.read()

requirements = []
setup(name = "liquidcosmo",
      version= "0.0.1",
      author = "Nils Schöneberg",
      description="Cosmological analyses made easy",
      long_description = readme,
      long_description_content_type="text/markdown",
      url="https://github.com/schoeneberg/liquidcosmo/",
      packages=find_packages(),
      install_requires=requirements,
      classifiers=[
          "Programming Language :: Python :: 3.8",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
      ],
)
